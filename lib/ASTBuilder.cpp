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

using namespace std;
using namespace antlr4;

#define visitAnd2ASTPN(expr) \
    any_cast<ASTNodePtr>(visit(expr))

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
        return visitAnd2ASTPN(ctx->statement());
    else if (ctx->scope())
        return visitAnd2ASTPN(ctx->scope());
    return nullptr;
}

any ASTBuilder::visitIncludeStatement(qasm3Parser::IncludeStatementContext *ctx)
{
    return (ASTNodePtr)include(ctx->StringLiteral()->getText());
}

any ASTBuilder::visitConstDeclarationStatement(qasm3Parser::ConstDeclarationStatementContext *ctx)
{
    return (ASTNodePtr)const_decl(
        ctx->scalarType()->getText(),
        ctx->Identifier()->getText(),
        ctx->declarationExpression()->getText());
}

any ASTBuilder::visitQuantumDeclarationStatement(qasm3Parser::QuantumDeclarationStatementContext *ctx)
{
    return (ASTNodePtr)qubit_decl(ctx->qubitType()->getText(), ctx->Identifier()->getText());
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
    vector<GateOperand> stmts;
    for (auto some : ctx->gateOperandList()->gateOperand())
    {
        auto ret = visit(some);
        assert(ret.type() == typeid(GateOperand));
        auto it = any_cast<GateOperand>(ret);
        stmts.push_back(it);
    }
    return (ASTNodePtr)gatecall_vec(ctx->Identifier()->getText(), stmts);
}

any ASTBuilder::visitIndexedIdentifier(qasm3Parser::IndexedIdentifierContext *ctx)
{
    vector<ASTNodeList> indices_vec;
    for (auto some : ctx->indexOperator())
    {
        auto any_val = visit(some);
        if (!any_val.has_value())
        {
            exit(30);
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
            exit(31);
        }
    }
    return gateoperand_vec(ctx->Identifier()->getText(), indices_vec);
}

any ASTBuilder::visitIndexExpression(qasm3Parser::IndexExpressionContext *ctx)
{
    return (ASTNodePtr)block(block(id("test:"), visitAnd2ASTPN(ctx->expression())), visitAnd2ASTPN(ctx->indexOperator()));
    // return visitAnd2ASTPN(ctx->indexOperator());
}

any ASTBuilder::visitIndexOperator(qasm3Parser::IndexOperatorContext *ctx)
{
    ASTNodeList stmts;
    for (auto &&i : ctx->expression())
        stmts.push_back(visitAnd2ASTPN(i));
    return stmts;
}

// any ASTBuilder::visitAssignmentStatement(qasm3Parser::AssignmentStatementContext *ctx){}

any ASTBuilder::visitForStatement(qasm3Parser::ForStatementContext *ctx)
{
    // TODO: ctx->scalarType();
    ASTNodeList stmts;

    for (auto some : ctx->statementOrScope()->scope()->statementOrScope())
    {
        auto res = visitAnd2ASTPN(some);
        stmts.push_back(res);
    }

    return (ASTNodePtr)fnode(
        id(ctx->Identifier()->getText()),
        visitAnd2ASTPN(ctx->rangeExpression()), // id(ctx->rangeExpression()->getText()),
        block_vec(stmts));
}

any ASTBuilder::visitRangeExpression(qasm3Parser::RangeExpressionContext *ctx)
{
    ASTNodePtr start = nullptr, end = nullptr, step = nullptr;

    assert(ctx->expression().size() >= 2);
    start = visitAnd2ASTPN(ctx->expression()[0]);
    end = visitAnd2ASTPN(ctx->expression()[1]);

    if (ctx->expression().size() == 3)
        step = visitAnd2ASTPN(ctx->expression()[2]);
    else if (ctx->expression().size() > 3)
        assert(false);
    else
        ; // pass

    return (ASTNodePtr)(step ? range(start, end, step) : range(start, end));
}
// 1+1;
any ASTBuilder::visitExpressionStatement(qasm3Parser::ExpressionStatementContext *ctx)
{
    return visitAnd2ASTPN(ctx->expression());
}
// (1+1) => 1+1
any ASTBuilder::visitParenthesisExpression(qasm3Parser::ParenthesisExpressionContext *ctx)
{
    return visitAnd2ASTPN(ctx->expression());
}

// TILDE | EXCLAMATION_POINT | MINUS
any ASTBuilder::visitUnaryExpression(qasm3Parser::UnaryExpressionContext *ctx)
{
    // TODO expr can be further enhanced for unary expressions by adding more structures.
    return (ASTNodePtr)expr("@" + ctx->op->getText(), id("0"), visitAnd2ASTPN(ctx->expression()));
}

#define DEFINE_EXPR_VISITOR(name)                                                           \
    any ASTBuilder::visit##name##Expression(qasm3Parser::name##ExpressionContext *ctx) \
    {                                                                                       \
        assert(ctx->expression().size() == 2);                                              \
        return (ASTNodePtr)expr(ctx->op->getText(),                                         \
                                visitAnd2ASTPN(ctx->expression()[0]),                       \
                                visitAnd2ASTPN(ctx->expression()[1]));                      \
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

// (x + 1) => (id(x) + 1)
any ASTBuilder::visitLiteralExpression(qasm3Parser::LiteralExpressionContext *ctx)
{
    // TODO ... type
    return (ASTNodePtr)(id(ctx->getText()));
}
