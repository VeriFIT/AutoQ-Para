#include "antlr4-runtime.h"
#include "g4/qasm3Lexer.h"
#include "g4/qasm3Parser.h"
#include "g4/qasm3ParserVisitor.h"
#include "g4/qasm3ParserListener.h"
#include "g4/qasm3ParserBaseVisitor.h"
#include "g4/qasm3ParserBaseListener.h"

// ANTLR4-generated API (qasm3ParserVisitor) for parsing definitions from https://openqasm.com/intro.html;
class ASTBuilder : public qasm3ParserBaseVisitor
{
public:
    std::any visitProgram(qasm3Parser::ProgramContext *ctx) override;
    std::any visitVersion(qasm3Parser::VersionContext *ctx) override;
    // ---- std::any visitStatement(qasm3Parser::StatementContext *ctx) override;
    /// std::any visitAnnotation(qasm3Parser::AnnotationContext *ctx) override;
    /// std::any visitScope(qasm3Parser::ScopeContext *ctx) override;
    // TODO std::any visitPragma(qasm3Parser::PragmaContext *ctx) override;
    std::any visitStatementOrScope(qasm3Parser::StatementOrScopeContext *ctx) override;
    // TODO std::any visitCalibrationGrammarStatement(qasm3Parser::CalibrationGrammarStatementContext *ctx) override;
    std::any visitIncludeStatement(qasm3Parser::IncludeStatementContext *ctx) override;
    // TODO std::any visitBreakStatement(qasm3Parser::BreakStatementContext *ctx) override;
    // TODO std::any visitContinueStatement(qasm3Parser::ContinueStatementContext *ctx) override;
    // TODO std::any visitEndStatement(qasm3Parser::EndStatementContext *ctx) override;
    std::any visitForStatement(qasm3Parser::ForStatementContext *ctx) override;
    // TODO std::any visitIfStatement(qasm3Parser::IfStatementContext *ctx) override;
    // TODO std::any visitReturnStatement(qasm3Parser::ReturnStatementContext *ctx) override;
    // TODO std::any visitWhileStatement(qasm3Parser::WhileStatementContext *ctx) override;
    // TODO std::any visitSwitchStatement(qasm3Parser::SwitchStatementContext *ctx) override;
    // TODO std::any visitSwitchCaseItem(qasm3Parser::SwitchCaseItemContext *ctx) override;
    // TODO std::any visitBarrierStatement(qasm3Parser::BarrierStatementContext *ctx) override;
    // TODO std::any visitBoxStatement(qasm3Parser::BoxStatementContext *ctx) override;
    // TODO std::any visitDelayStatement(qasm3Parser::DelayStatementContext *ctx) override;
    // TODO std::any visitNopStatement(qasm3Parser::NopStatementContext *ctx) override;
    std::any visitGateCallStatement(qasm3Parser::GateCallStatementContext *ctx) override;
    // TODO std::any visitMeasureArrowAssignmentStatement(qasm3Parser::MeasureArrowAssignmentStatementContext *ctx) override;
    // TODO std::any visitResetStatement(qasm3Parser::ResetStatementContext *ctx) override;
    // TODO std::any visitAliasDeclarationStatement(qasm3Parser::AliasDeclarationStatementContext *ctx) override;
    // TODO std::any visitClassicalDeclarationStatement(qasm3Parser::ClassicalDeclarationStatementContext *ctx) override;
    std::any visitConstDeclarationStatement(qasm3Parser::ConstDeclarationStatementContext *ctx) override;
    // TODO std::any visitIoDeclarationStatement(qasm3Parser::IoDeclarationStatementContext *ctx) override;
    // TODO std::any visitOldStyleDeclarationStatement(qasm3Parser::OldStyleDeclarationStatementContext *ctx) override;
    std::any visitQuantumDeclarationStatement(qasm3Parser::QuantumDeclarationStatementContext *ctx) override;
    // TODO std::any visitDefStatement(qasm3Parser::DefStatementContext *ctx) override;
    // TODO std::any visitExternStatement(qasm3Parser::ExternStatementContext *ctx) override;
    // TODO std::any visitGateStatement(qasm3Parser::GateStatementContext *ctx) override;
    /// std::any visitAssignmentStatement(qasm3Parser::AssignmentStatementContext *ctx) override;
    std::any visitExpressionStatement(qasm3Parser::ExpressionStatementContext *ctx) override;
    // TODO std::any visitCalStatement(qasm3Parser::CalStatementContext *ctx) override;
    // TODO std::any visitDefcalStatement(qasm3Parser::DefcalStatementContext *ctx) override;
    std::any visitBitwiseXorExpression(qasm3Parser::BitwiseXorExpressionContext *ctx) override;
    std::any visitAdditiveExpression(qasm3Parser::AdditiveExpressionContext *ctx) override;
    // TODO std::any visitDurationofExpression(qasm3Parser::DurationofExpressionContext *ctx) override;
    std::any visitParenthesisExpression(qasm3Parser::ParenthesisExpressionContext *ctx) override;
    std::any visitComparisonExpression(qasm3Parser::ComparisonExpressionContext *ctx) override;
    std::any visitMultiplicativeExpression(qasm3Parser::MultiplicativeExpressionContext *ctx) override;
    std::any visitLogicalOrExpression(qasm3Parser::LogicalOrExpressionContext *ctx) override;
    // TODO std::any visitCastExpression(qasm3Parser::CastExpressionContext *ctx) override;
    std::any visitPowerExpression(qasm3Parser::PowerExpressionContext *ctx) override;
    std::any visitBitwiseOrExpression(qasm3Parser::BitwiseOrExpressionContext *ctx) override;
    // TODO std::any visitCallExpression(qasm3Parser::CallExpressionContext *ctx) override;
    std::any visitBitshiftExpression(qasm3Parser::BitshiftExpressionContext *ctx) override;
    std::any visitBitwiseAndExpression(qasm3Parser::BitwiseAndExpressionContext *ctx) override;
    std::any visitEqualityExpression(qasm3Parser::EqualityExpressionContext *ctx) override;
    std::any visitLogicalAndExpression(qasm3Parser::LogicalAndExpressionContext *ctx) override;
    std::any visitIndexExpression(qasm3Parser::IndexExpressionContext *ctx) override;
    std::any visitUnaryExpression(qasm3Parser::UnaryExpressionContext *ctx) override;
    std::any visitLiteralExpression(qasm3Parser::LiteralExpressionContext *ctx) override;
    // TODO std::any visitAliasExpression(qasm3Parser::AliasExpressionContext *ctx) override;
    // TODO std::any visitDeclarationExpression(qasm3Parser::DeclarationExpressionContext *ctx) override;
    // TODO std::any visitMeasureExpression(qasm3Parser::MeasureExpressionContext *ctx) override;
    std::any visitRangeExpression(qasm3Parser::RangeExpressionContext *ctx) override;
    // TODO std::any visitSetExpression(qasm3Parser::SetExpressionContext *ctx) override;
    // TODO std::any visitArrayLiteral(qasm3Parser::ArrayLiteralContext *ctx) override;
    std::any visitIndexOperator(qasm3Parser::IndexOperatorContext *ctx) override;
    std::any visitIndexedIdentifier(qasm3Parser::IndexedIdentifierContext *ctx) override;
    // TODO std::any visitReturnSignature(qasm3Parser::ReturnSignatureContext *ctx) override;
    // TODO std::any visitGateModifier(qasm3Parser::GateModifierContext *ctx) override;
    // TODO std::any visitScalarType(qasm3Parser::ScalarTypeContext *ctx) override;
    // TODO std::any visitQubitType(qasm3Parser::QubitTypeContext *ctx) override;
    // TODO std::any visitArrayType(qasm3Parser::ArrayTypeContext *ctx) override;
    // TODO std::any visitArrayReferenceType(qasm3Parser::ArrayReferenceTypeContext *ctx) override;
    // TODO std::any visitDesignator(qasm3Parser::DesignatorContext *ctx) override;
    // TODO std::any visitDefcalTarget(qasm3Parser::DefcalTargetContext *ctx) override;
    // TODO std::any visitDefcalArgumentDefinition(qasm3Parser::DefcalArgumentDefinitionContext *ctx) override;
    // TODO std::any visitDefcalOperand(qasm3Parser::DefcalOperandContext *ctx) override;
    /// std::any visitGateOperand(qasm3Parser::GateOperandContext *ctx) override;
    // TODO std::any visitExternArgument(qasm3Parser::ExternArgumentContext *ctx) override;
    // TODO std::any visitArgumentDefinition(qasm3Parser::ArgumentDefinitionContext *ctx) override;
    // TODO std::any visitArgumentDefinitionList(qasm3Parser::ArgumentDefinitionListContext *ctx) override;
    // TODO std::any visitDefcalArgumentDefinitionList(qasm3Parser::DefcalArgumentDefinitionListContext *ctx) override;
    // TODO std::any visitDefcalOperandList(qasm3Parser::DefcalOperandListContext *ctx) override;
    // TODO std::any visitExpressionList(qasm3Parser::ExpressionListContext *ctx) override;
    // TODO std::any visitIdentifierList(qasm3Parser::IdentifierListContext *ctx) override;
    /// std::any visitGateOperandList(qasm3Parser::GateOperandListContext *ctx) override;
    // TODO std::any visitExternArgumentList(qasm3Parser::ExternArgumentListContext *ctx) override;
};