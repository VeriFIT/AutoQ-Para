#include <iostream>
#include "antlr4-runtime.h"
#include "qasm3Lexer.h"
#include "qasm3Parser.h"
#include "qasm3ParserVisitor.h"
#include "qasm3ParserListener.h"
#include "qasm3ParserBaseVisitor.h"
#include "qasm3ParserBaseListener.h"
#include "ASTBuilder.hpp"
#include "ASTNode.hpp"
#include <stdexcept>

using namespace std;
using namespace antlr4;

static ASTNodePtr gateOperandToAst(const GateOperand &op)
{
    ASTNodePtr base = id(op.first);
    if (!op.second.has_value() || op.second->empty())
        return base;
    return index_expr(base, *op.second);
}


any ASTBuilder::visitProgram(qasm3Parser::ProgramContext *ctx)
{
    auto root = program();
    for (auto child : ctx->children)
    {
        auto childAstAny = visit(child);
        if (childAstAny.has_value() && childAstAny.type() == typeid(ASTNodePtr))
            root->statements.push_back(any_cast<ASTNodePtr>(childAstAny));
        else if (dynamic_cast<antlr4::tree::TerminalNode *>(child))
            ; // pass // cout << "Terminal: " << term->getText() << endl;
        else
            cerr << "bad any_cast: type is " << childAstAny.type().name() << endl;
    }
    return (ASTNodePtr)root;
}

any ASTBuilder::visitVersion(qasm3Parser::VersionContext *ctx)
{
    return (ASTNodePtr)version(ctx->VersionSpecifier()->getText());
}

// any ASTBuilder::visitAnnotation(qasm3Parser::AnnotationContext *ctx) {}
// any ASTBuilder::visitStatement(qasm3Parser::StatementContext *ctx){ return visitChildren(ctx); }

any ASTBuilder::visitStatementOrScope(qasm3Parser::StatementOrScopeContext *ctx)
{
    if (ctx->statement())
        return toAst(ctx->statement());
    else if (ctx->scope())
        return toAst(ctx->scope());
    return nullptr;
}

any ASTBuilder::visitIncludeStatement(qasm3Parser::IncludeStatementContext *ctx)
{
    return (ASTNodePtr)include(ctx->StringLiteral()->getText());
}

any ASTBuilder::visitConstDeclarationStatement(qasm3Parser::ConstDeclarationStatementContext *ctx)
{
    ASTNodePtr expr = nullptr;
    if (ctx->declarationExpression())
        expr = toAst(ctx->declarationExpression()->expression());
    return (ASTNodePtr)const_decl(
        toAst(ctx->scalarType()),
        ctx->Identifier()->getText(),
        expr);
}

any ASTBuilder::visitQuantumDeclarationStatement(qasm3Parser::QuantumDeclarationStatementContext *ctx)
{
    ASTNodePtr size = nullptr;
    if (ctx->qubitType()->designator())
        size = toAst(ctx->qubitType()->designator()->expression());
    return (ASTNodePtr)qubit_decl(ctx->Identifier()->getText(), size);
}

// --------------------------------------------------------
// |   ccx  q, a[0][0]          -- Statement               |
// |      ↓  enter                                         |
// |  (ccx)(q  ,  a   [0][0] )  -- GateCallStatement       |
// |      ↓  enter                                         |
// |   ccx (q) , (a   [0][0] )  -- GateOperandList (V)isit |
// |      ↓  enter                                         |
// |   ccx  q(), (a) ([0][0] )  -- visitIndexedIdentifier  |
// |      ↓                                                |
// |   (a) (...)                -- Simplify                |
// --------------------------------------------------------
any ASTBuilder::visitGateCallStatement(qasm3Parser::GateCallStatementContext *ctx)
{
    // gateOperand: indexedIdentifier | HardwareQubit;
    // ! HardwareQubit statement is not yet implemented
    vector<GateOperand> operands;
    for (auto opCtx : ctx->gateOperandList()->gateOperand())
    {
        auto ret = visit(opCtx);
        assert(ret.type() == typeid(GateOperand));
        operands.push_back(any_cast<GateOperand>(ret));
    }
    return (ASTNodePtr)gatecall_vec(ctx->Identifier()->getText(), operands);
}

any ASTBuilder::visitMeasureArrowAssignmentStatement(qasm3Parser::MeasureArrowAssignmentStatementContext *ctx)
{
    auto measureOp = any_cast<GateOperand>(visit(ctx->measureExpression()->gateOperand()));
    ASTNodePtr qubit = gateOperandToAst(measureOp);

    ASTNodePtr target = nullptr;
    if (ctx->indexedIdentifier())
    {
        auto tgt = any_cast<GateOperand>(visit(ctx->indexedIdentifier()));
        target = gateOperandToAst(tgt);
    }

    return (ASTNodePtr)measure_stmt(qubit, target);
}

any ASTBuilder::visitGateStatement(qasm3Parser::GateStatementContext *ctx)
{
    ASTNodeList params;
    if (ctx->params)
        for (auto token : ctx->params->Identifier())
            params.push_back(id(token->getText()));

    ASTNodeList qubits;
    for (auto token : ctx->qubits->Identifier())
        qubits.push_back(id(token->getText()));

    ASTNodeList stmts;
    for (auto node : ctx->scope()->statementOrScope())
        stmts.push_back(toAst(node));

    return (ASTNodePtr)gate_stmt(ctx->Identifier()->getText(), params, qubits, block_vec(stmts));
}

any ASTBuilder::visitIndexedIdentifier(qasm3Parser::IndexedIdentifierContext *ctx)
{
    vector<ASTNodeList> indices_vec;
    for (auto opCtx : ctx->indexOperator())
    {
        auto any_val = visit(opCtx);
        if (!any_val.has_value())
        {
            throw runtime_error("visitIndexedIdentifier: empty index expression");
        }
        else if (any_val.type() == typeid(ASTNodeList))
        {
            indices_vec.push_back(any_cast<ASTNodeList>(any_val));
        }
        else if (any_val.type() == typeid(ASTNodePtr))
        {
            ASTNodeList tmp;
            tmp.push_back(any_cast<ASTNodePtr>(any_val));
            indices_vec.push_back(tmp);
        }
        else
        {
            throw runtime_error("visitIndexedIdentifier: unexpected return type from visit(indexOperator)");
        }
    }
    return gateoperand_vec(ctx->Identifier()->getText(), indices_vec);
}

any ASTBuilder::visitIndexExpression(qasm3Parser::IndexExpressionContext *ctx)
{
    auto base = toAst(ctx->expression());
    auto idxAny = visit(ctx->indexOperator());
    vector<ASTNodeList> indices;
    if (idxAny.type() == typeid(ASTNodeList))
    {
        indices.push_back(any_cast<ASTNodeList>(idxAny));
    }
    else if (idxAny.type() == typeid(ASTNodePtr))
    {
        ASTNodeList tmp;
        tmp.push_back(any_cast<ASTNodePtr>(idxAny));
        indices.push_back(tmp);
    }
    else
    {
        // unreachable in normal grammar
    }
    return (ASTNodePtr)index_expr(base, indices);
}

any ASTBuilder::visitIndexOperator(qasm3Parser::IndexOperatorContext *ctx)
{
    ASTNodeList stmts;
    for (auto &&exprCtx : ctx->expression())
        stmts.push_back(toAst(exprCtx));
    return stmts;
}

// any ASTBuilder::visitAssignmentStatement(qasm3Parser::AssignmentStatementContext *ctx){}

any ASTBuilder::visitForStatement(qasm3Parser::ForStatementContext *ctx)
{
    ASTNodeList stmts;

    for (auto node : ctx->statementOrScope()->scope()->statementOrScope())
    {
        auto res = toAst(node);
        stmts.push_back(res);
    }

    return (ASTNodePtr)fnode(
        toAst(ctx->scalarType()),
        id(ctx->Identifier()->getText()),
        toAst(ctx->rangeExpression()),
        block_vec(stmts));
}

any ASTBuilder::visitIfStatement(qasm3Parser::IfStatementContext *ctx)
{
    auto cond = toAst(ctx->expression());

    ASTNodePtr then_b = nullptr;
    ASTNodePtr else_b = nullptr;

    if (ctx->if_body->scope())
    {
        ASTNodeList stmts;
        for (auto node : ctx->if_body->scope()->statementOrScope())
            stmts.push_back(toAst(node));
        then_b = block_vec(stmts);
    }
    else
    {
        then_b = toAst(ctx->if_body);
    }

    if (ctx->else_body)
    {
        if (ctx->else_body->scope())
        {
            ASTNodeList stmts;
            for (auto node : ctx->else_body->scope()->statementOrScope())
                stmts.push_back(toAst(node));
            else_b = block_vec(stmts);
        }
        else
        {
            else_b = toAst(ctx->else_body);
        }
    }

    return (ASTNodePtr)ifnode(cond, then_b, else_b);
}

any ASTBuilder::visitRangeExpression(qasm3Parser::RangeExpressionContext *ctx)
{
    ASTNodePtr start = nullptr, end = nullptr, step = nullptr;

    assert(ctx->expression().size() >= 2);
    start = toAst(ctx->expression()[0]);
    end = toAst(ctx->expression()[1]);

    if (ctx->expression().size() == 3)
        step = toAst(ctx->expression()[2]);
    else if (ctx->expression().size() > 3)
        assert(false);
    else
        ; // pass

    return (ASTNodePtr)(step ? range(start, end, step) : range(start, end));
}
// 1+1;
any ASTBuilder::visitExpressionStatement(qasm3Parser::ExpressionStatementContext *ctx)
{
    return toAst(ctx->expression());
}
// (1+1) => 1+1
any ASTBuilder::visitParenthesisExpression(qasm3Parser::ParenthesisExpressionContext *ctx)
{
    return toAst(ctx->expression());
}

// TILDE | EXCLAMATION_POINT | MINUS
any ASTBuilder::visitUnaryExpression(qasm3Parser::UnaryExpressionContext *ctx)
{
    // TODO expr can be further enhanced for unary expressions by adding more structures.
    return (ASTNodePtr)expr("@" + ctx->op->getText(), id("0"), toAst(ctx->expression()));
}

#define DEFINE_EXPR_VISITOR(name)                                                           \
    any ASTBuilder::visit##name##Expression(qasm3Parser::name##ExpressionContext *ctx) \
    {                                                                                       \
        assert(ctx->expression().size() == 2);                                              \
        return (ASTNodePtr)expr(ctx->op->getText(),                                         \
                                toAst(ctx->expression()[0]),                       \
                                toAst(ctx->expression()[1]));                      \
    }

//  DOUBLE_ASTERISK: **
DEFINE_EXPR_VISITOR(Power)
// ASTERISK | SLASH | PERCENT
// 1*1*1 => ((1*1)*1)
DEFINE_EXPR_VISITOR(Multiplicative)
// PLUS | MINUS
// 1+1+1 => ((1+1)+1)
DEFINE_EXPR_VISITOR(Additive)
// >> | <<
DEFINE_EXPR_VISITOR(Bitshift)
// > | < | >= | <=
DEFINE_EXPR_VISITOR(Comparison)
// == | !=
DEFINE_EXPR_VISITOR(Equality)
// AMPERSAND: &
DEFINE_EXPR_VISITOR(BitwiseAnd)
// CARET: ^
DEFINE_EXPR_VISITOR(BitwiseXor)
// PIPE: |
DEFINE_EXPR_VISITOR(BitwiseOr)
// DOUBLE_AMPERSAND: &&
DEFINE_EXPR_VISITOR(LogicalAnd)
// DOUBLE_PIPE: ||
DEFINE_EXPR_VISITOR(LogicalOr)

any ASTBuilder::visitCallExpression(qasm3Parser::CallExpressionContext *ctx)
{
    ASTNodePtr callee = id(ctx->Identifier()->getText());
    ASTNodeList args;
    if (ctx->expressionList())
    {
        for (auto exprCtx : ctx->expressionList()->expression())
            args.push_back(toAst(exprCtx));
    }
    return (ASTNodePtr)call_expr(callee, args);
}

// (x + 1) => (id(x) + 1)
any ASTBuilder::visitLiteralExpression(qasm3Parser::LiteralExpressionContext *ctx)
{
    // TODO ... type
    return (ASTNodePtr)(id(ctx->getText()));
}

any ASTBuilder::visitScalarType(qasm3Parser::ScalarTypeContext *ctx)
{
    return (ASTNodePtr)type_node(ctx->getText());
}

any ASTBuilder::visitArrayType(qasm3Parser::ArrayTypeContext *ctx)
{
    return (ASTNodePtr)type_node(ctx->getText());
}

any ASTBuilder::visitClassicalDeclarationStatement(qasm3Parser::ClassicalDeclarationStatementContext *ctx)
{
    ASTNodePtr expr = nullptr;
    if (ctx->declarationExpression())
        expr = toAst(ctx->declarationExpression()->expression());
    ASTNodePtr type = ctx->scalarType() ? toAst(ctx->scalarType()) : toAst(ctx->arrayType());
    return (ASTNodePtr)classical_decl(type, ctx->Identifier()->getText(), expr);
}

