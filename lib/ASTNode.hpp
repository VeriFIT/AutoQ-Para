#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <string>
#include <optional>

using namespace std;

enum ASTKind
{
    Identifier,
    Int,
    Expr,
    Range,
    Program,
    Block,
    GateCall,
    GateStmt,
    For,
    If,
    Include,
    Version,
    ConstDecl,
    ClassicalDecl,
    QubitDecl,
    MeasureStmt,
    Index,
    CallExpr,
    Literal,
    Type,
    Unknown
};

struct ASTNode
{
    virtual ~ASTNode() = default;
    virtual void print(int indent = 0) const = 0;
    virtual ASTKind kind() const { return ASTKind::Unknown; }
};

struct CanShow
{
    virtual string toString() const = 0;
};

struct IdentifierNode : ASTNode, CanShow
{
    string name;
    IdentifierNode(string n) : name(n) {}
    ASTKind kind() const override { return ASTKind::Identifier; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Identifier(" << name << ")" << endl;
    }

    string toString() const override
    {
        return name;
    }
};

struct IntNode : ASTNode
{
    int value;
    IntNode(int v) : value(v) {}
    ASTKind kind() const override { return ASTKind::Int; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Int(" << value << ")" << endl;
    }
};

using ASTNodePtr = shared_ptr<ASTNode>;
using ASTNodeList = vector<ASTNodePtr>;

struct ExprNode : ASTNode, CanShow
{
    string op;
    ASTNodePtr expr1;
    ASTNodePtr expr2;
    ExprNode(string _op, ASTNodePtr e1, ASTNodePtr e2) : op(_op), expr1(e1), expr2(e2) {}
    ASTKind kind() const override { return ASTKind::Expr; }
    string toString() const override
    {
        auto e1 = dynamic_pointer_cast<CanShow>(expr1);
        auto s1 = e1 ? e1->toString() : "Expr( ?... )";
        auto e2 = dynamic_pointer_cast<CanShow>(expr2);
        auto s2 = e2 ? e2->toString() : "Expr( ?... )";
        // auto e2 = expr2->kind() == ASTKind::Expr ? dynamic_pointer_cast<CanShow>(expr2)->toString() : "Expr( ?... )";
        return "_(" + s1 + " " + op + " " + s2 + ")";
    }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Expr" << (*this).toString().substr(1) << endl;
    }
};

struct RangeNode : ASTNode
{
    ASTNodePtr start, end, step;
    RangeNode(ASTNodePtr s, ASTNodePtr e, ASTNodePtr st = nullptr)
        : start(s), end(e), step(st) {}
    ASTKind kind() const override { return ASTKind::Range; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Range start: " << dynamic_pointer_cast<CanShow>(start)->toString() << endl;
        cout << string(indent, ' ') << "        end: " << dynamic_pointer_cast<CanShow>(end)->toString() << endl;
        cout << string(indent, ' ') << "       step: " << (step ? dynamic_pointer_cast<CanShow>(step)->toString() : "#NULL") << endl;
    }
};

struct ProgramNode : ASTNode
{
    ASTNodeList statements;
    ProgramNode(ASTNodeList stmts)
        : statements(stmts) {}
    ASTKind kind() const override { return ASTKind::Program; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Program" << endl;
        for (const auto &stmt : statements)
            stmt->print(indent + 2);
    }
};

struct BlockNode : ASTNode
{
    ASTNodeList statements;
    BlockNode(ASTNodeList stmts)
        : statements(stmts) {}
    ASTKind kind() const override { return ASTKind::Block; }
    void print(int indent = 0) const override
    {
        if (statements.empty())
        {
            cout << string(indent, ' ') << "Block#NULL>" << endl;
            return;
        }
        else
        {
            cout << string(indent, ' ') << "Block {" << endl;
            for (const auto &stmt : statements)
                stmt->print(indent + 2);
            cout << string(indent, ' ') << "}Block end" << endl;
        }
    }
};

using GateOperand = pair<string, optional<vector<vector<shared_ptr<ASTNode>>>>>;

struct GateCallNode : ASTNode
{
    string name;
    ASTNodeList gate_args;
    vector<GateOperand> operand_list;
    GateCallNode(string g, ASTNodeList args, vector<GateOperand> operand_list)
        : name(g), gate_args(move(args)), operand_list(move(operand_list)) {}
    GateCallNode(string g, vector<GateOperand> operand_list)
        : name(g), operand_list(move(operand_list)) {}
    ASTKind kind() const override { return ASTKind::GateCall; }
    void print(int indent = 0) const override
    {
        auto operandToString = [](const GateOperand &arg) {
            string s = arg.first;
            if (!arg.second.has_value() || arg.second->empty() || arg.second->front().empty())
            {
                s += "#NULL";
            }
            else
            {
                for (const auto &exps : arg.second.value())
                {
                    s += "[";
                    for (size_t k = 0; k < exps.size(); k++)
                    {
                        auto e1 = dynamic_pointer_cast<CanShow>(exps[k]);
                        s += e1 ? e1->toString() : "Expr( ?... )";
                        if (k + 1 < exps.size())
                            s += ", ";
                    }
                    s += "]";
                }
            }
            return s;
        };

        cout << string(indent, ' ') << "GateCall(" << name;
        cout << ", gate_argv=(";
        for (size_t i = 0; i < gate_args.size(); ++i)
        {
            auto s = dynamic_pointer_cast<CanShow>(gate_args[i]);
            cout << (s ? s->toString() : "?");
            if (i + 1 < gate_args.size())
                cout << ",";
        }
        cout << ")";
        cout << ", argv=(";
        for (size_t i = 0; i < operand_list.size(); ++i)
        {
            cout << operandToString(operand_list[i]);
            if (i + 1 < operand_list.size())
                cout << ",";
        }
        cout << "))" << endl;
    }
};

struct GateStmtNode : ASTNode
{
    string name;
    ASTNodeList params;
    ASTNodeList qubits;
    ASTNodePtr body;
    GateStmtNode(string n, ASTNodeList p, ASTNodeList q, ASTNodePtr b)
        : name(move(n)), params(move(p)), qubits(move(q)), body(move(b)) {}
    ASTKind kind() const override { return ASTKind::GateStmt; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "GateStmt(" << name << ")" << endl;

        if (!params.empty())
        {
            cout << string(indent + 2, ' ') << "params:";
            for (const auto &p : params)
            {
                auto id = dynamic_pointer_cast<IdentifierNode>(p);
                cout << " " << (id ? id->toString() : "?");
            }
            cout << endl;
        }
        if (!qubits.empty())
        {
            cout << string(indent + 2, ' ') << "qubits:";
            for (const auto &q : qubits)
            {
                auto id = dynamic_pointer_cast<IdentifierNode>(q);
                cout << " " << (id ? id->toString() : "?");
            }
            cout << endl;
        }
        if (body)
        {
            cout << string(indent + 2, ' ') << "--- Body" << endl;
            body->print(indent + 2);
        }
    }
};

struct ForNode : ASTNode
{
    ASTNodePtr type;
    ASTNodePtr bind_target;
    ASTNodePtr range;
    ASTNodePtr body;
    ForNode(ASTNodePtr t, ASTNodePtr i, ASTNodePtr r, ASTNodePtr b)
        : type(t), bind_target(i), range(r), body(b) {}
    ASTKind kind() const override { return ASTKind::For; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "For" << endl;
        if (bind_target && type)
            cout << string(indent + 2, ' ')
                 << "--- Binding(" << dynamic_pointer_cast<CanShow>(bind_target)->toString() << ":" << dynamic_pointer_cast<CanShow>(type)->toString() << ")" << endl;
        if (range)
        {
            cout << string(indent + 2, ' ') << "--- Range" << endl;
            range->print(indent + 2);
        }
        if (body)
        {
            cout << string(indent + 2, ' ') << "--- Body" << endl;
            body->print(indent + 2);
        }
    }
};

struct IfNode : ASTNode
{
    ASTNodePtr cond;
    ASTNodePtr then_body;
    ASTNodePtr else_body;
    IfNode(ASTNodePtr c, ASTNodePtr t, ASTNodePtr e = nullptr)
        : cond(c), then_body(t), else_body(e) {}
    ASTKind kind() const override { return ASTKind::If; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "If" << endl;
        if (cond)
        {
            cout << string(indent + 2, ' ') << "--- Condition" << endl;
            cond->print(indent + 2);
        }
        if (then_body)
        {
            cout << string(indent + 2, ' ') << "--- Then" << endl;
            then_body->print(indent + 2);
        }
        if (else_body)
        {
            cout << string(indent + 2, ' ') << "--- Else" << endl;
            else_body->print(indent + 2);
        }
    }
};

struct IncludeNode : ASTNode
{
    string filename;
    IncludeNode(string f) : filename(f) {}
    ASTKind kind() const override { return ASTKind::Include; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Include(" << filename << ")" << endl;
    }
};

struct VersionNode : ASTNode
{
    string version;
    VersionNode(string v) : version(v) {}
    ASTKind kind() const override { return ASTKind::Version; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Version(" << version << ")\n";
    }
};

struct ConstDeclNode : ASTNode
{
    ASTNodePtr type;
    string name;
    ASTNodePtr expr;
    ConstDeclNode(ASTNodePtr t, string n, ASTNodePtr e)
        : type(move(t)), name(n), expr(e) {}
    ASTKind kind() const override { return ASTKind::ConstDecl; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "ConstDecl(";
        auto show = dynamic_pointer_cast<CanShow>(type);
        cout << (show ? show->toString() : "?") << " " << name;
        cout << " = " << dynamic_pointer_cast<CanShow>(expr)->toString();
        cout << ")" << endl;
    }
};

struct QubitDeclNode : ASTNode
{
    string name;
    ASTNodePtr size;
    QubitDeclNode(string n, ASTNodePtr s = nullptr) : name(n), size(s) {}
    ASTKind kind() const override { return ASTKind::QubitDecl; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ')
             << "QubitDecl" << (size ? "[ " + dynamic_pointer_cast<CanShow>(size)->toString() + " ]" : "[ 1 ]")
             << "( " << name << " )" << endl;
    }
};

struct MeasureStmtNode : ASTNode
{
    ASTNodePtr qubit;
    ASTNodePtr target;
    MeasureStmtNode(ASTNodePtr q, ASTNodePtr t = nullptr) : qubit(q), target(t) {}
    ASTKind kind() const override { return ASTKind::MeasureStmt; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Measure(";
        auto q = dynamic_pointer_cast<CanShow>(qubit);
        cout << (q ? q->toString() : "?");
        if (target)
        {
            auto t = dynamic_pointer_cast<CanShow>(target);
            cout << " -> " << (t ? t->toString() : "?");
        }
        cout << ")" << endl;
    }
};

struct IndexNode : ASTNode, CanShow
{
    ASTNodePtr base;
    vector<ASTNodeList> indices;
    IndexNode(ASTNodePtr b, vector<ASTNodeList> idx) : base(move(b)), indices(move(idx)) {}
    ASTKind kind() const override { return ASTKind::Index; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Index(";
        auto show = dynamic_pointer_cast<CanShow>(base);
        cout << (show ? show->toString() : "?");
        for (const auto &idxList : indices)
        {
            cout << "[";
            for (size_t i = 0; i < idxList.size(); ++i)
            {
                auto cs = dynamic_pointer_cast<CanShow>(idxList[i]);
                cout << (cs ? cs->toString() : "?");
                if (i + 1 < idxList.size())
                    cout << ", ";
            }
            cout << "]";
        }
        cout << ")" << endl;
    }
    string toString() const override
    {
        string s;
        auto show = dynamic_pointer_cast<CanShow>(base);
        s += show ? show->toString() : "?";
        for (const auto &idxList : indices)
        {
            s += "[";
            for (size_t i = 0; i < idxList.size(); ++i)
            {
                auto cs = dynamic_pointer_cast<CanShow>(idxList[i]);
                s += cs ? cs->toString() : "?";
                if (i + 1 < idxList.size())
                    s += ", ";
            }
            s += "]";
        }
        return s;
    }
};

struct CallExprNode : ASTNode, CanShow
{
    ASTNodePtr callee;
    ASTNodeList arguments;
    CallExprNode(ASTNodePtr c, ASTNodeList args)
        : callee(move(c)), arguments(move(args)) {}
    ASTKind kind() const override { return ASTKind::CallExpr; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "CallExpr(";
        auto show = dynamic_pointer_cast<CanShow>(callee);
        cout << (show ? show->toString() : "?");
        cout << ", args=(";
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            auto s = dynamic_pointer_cast<CanShow>(arguments[i]);
            cout << (s ? s->toString() : "?");
            if (i + 1 < arguments.size())
                cout << ", ";
        }
        cout << "))" << endl;
    }
    string toString() const override
    {
        string s;
        auto show = dynamic_pointer_cast<CanShow>(callee);
        s += show ? show->toString() : "?";
        s += "(";
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            auto cs = dynamic_pointer_cast<CanShow>(arguments[i]);
            s += cs ? cs->toString() : "?";
            if (i + 1 < arguments.size())
                s += ", ";
        }
        s += ")";
        return s;
    }
};

struct TypeNode : ASTNode, CanShow
{
    string name;
    TypeNode(string n) : name(n) {}
    ASTKind kind() const override { return ASTKind::Type; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "Type(" << name << ")" << endl;
    }
    string toString() const override { return "Type(" + name + ")" ; }
};

struct ClassicalDeclNode : ASTNode
{
    ASTNodePtr type;
    string name;
    ASTNodePtr expr;
    ClassicalDeclNode(ASTNodePtr t, string n, ASTNodePtr e = nullptr)
        : type(move(t)), name(n), expr(e) {}
    ASTKind kind() const override { return ASTKind::ClassicalDecl; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "ClassicalDecl(";
        auto show = dynamic_pointer_cast<CanShow>(type);
        cout << (show ? show->toString() : "?") << " " << name;
        if (expr)
            cout << " = " << dynamic_pointer_cast<CanShow>(expr)->toString();
        cout << ")" << endl;
    }
};

// __VA_ARGS__
#define MAKE_NODE_FACTORY(name, type, ...)                            \
    template <typename... Args>                                       \
    auto name(Args &&...args)                                         \
    {                                                                 \
        return make_shared<type>(make_nodes(forward<Args>(args)...)); \
    }
#define MAKE_NODE_FACTORY_INLINE(name, type, ...) \
    template <typename... Args>                   \
    inline auto name(Args &&...args)              \
    {                                             \
        return make_shared<type>(args...);        \
    }
// cout << (#__VA_ARGS__ << ...) << endl;

template <typename... Args>
ASTNodeList make_nodes(Args &&...args)
{
    ASTNodeList v;
    (v.push_back(args), ...);
    return v;
}

// ------------------- For expanding into the following
// template <typename... Args>
// auto program(Args &&...args)
// {
//     return make_shared<ProgramNode>(make_nodes(forward<Args>(args)...));
// }
MAKE_NODE_FACTORY(program, ProgramNode)
MAKE_NODE_FACTORY(block, BlockNode)

template <typename... Exps>
vector<shared_ptr<ASTNode>> exps(Exps &&...exps)
{
    vector<shared_ptr<ASTNode>> v;
    (v.push_back(forward<Exps>(exps)), ...);
    return v;
}

template <typename... Exprss>
GateOperand gateoperand(const string &key, Exprss &&...exprs)
{
    vector<vector<shared_ptr<ASTNode>>> v;
    (v.push_back(forward<Exprss>(exprs)), ...);
    return {key, v.empty() ? nullopt : make_optional(v)};
}

template <typename... Operands>
auto gatecall(const string &gate, const ASTNodeList &args, Operands &&...ops)
{
    vector<GateOperand> v;
    (v.push_back(forward<Operands>(ops)), ...);
    return make_shared<GateCallNode>(gate, args, move(v));
}

template <typename... Operands>
auto gatecall(const string &gate, Operands &&...ops)
{
    return gatecall(gate, {}, forward<Operands>(ops)...);
}

inline ASTNodeList astn_vec(const ASTNodeList &v)
{
    return v;
}

inline GateOperand gateoperand_vec(const string &key, const vector<vector<ASTNodePtr>> &v)
{
    return {key, v.empty() ? nullopt : make_optional(v)};
}

inline shared_ptr<GateCallNode> gatecall_vec(const string &gate, const ASTNodeList &args,
                                             const vector<GateOperand> &gq_list)
{
    return make_shared<GateCallNode>(gate, args, gq_list);
}

inline shared_ptr<GateCallNode> gatecall_vec(const string &gate, const vector<GateOperand> &gq_list)
{
    return make_shared<GateCallNode>(gate, ASTNodeList{}, gq_list);
}

inline shared_ptr<GateStmtNode> gate_stmt(const string &name, const ASTNodeList &params,
                                          const ASTNodeList &qubits, ASTNodePtr body)
{
    return make_shared<GateStmtNode>(name, params, qubits, body);
}

MAKE_NODE_FACTORY_INLINE(type_node, TypeNode, string)
MAKE_NODE_FACTORY_INLINE(classical_decl, ClassicalDeclNode, ASTNodePtr, string, ASTNodePtr)

template <typename T, typename I, typename R, typename B>
auto fnode(T &&t, I &&iter, R &&rng, B &&body)
{
    return make_shared<ForNode>(t, iter, rng, body);
}

template <typename A, typename B, typename C>
auto ifnode(A &&cond, B &&then_b, C &&else_b)
{
    return make_shared<IfNode>(cond, then_b, else_b);
}

// ------------------- For expanding into the following
// template <typename Arg>
// auto block_vec(Arg &&arg)
// {
//     return make_shared<BlockNode>(arg);
// }
MAKE_NODE_FACTORY_INLINE(block_vec, BlockNode, ASTNodePtr)
// inline auto id(string n)
// {
//      return make_shared<IdentifierNode>(n);
// }
MAKE_NODE_FACTORY_INLINE(id, IdentifierNode, string)
MAKE_NODE_FACTORY_INLINE(num, IntNode, int)
MAKE_NODE_FACTORY_INLINE(expr, ExprNode, string, ASTNodePtr, ASTNodePtr)
MAKE_NODE_FACTORY_INLINE(range, RangeNode, ASTNodePtr, ASTNodePtr)
MAKE_NODE_FACTORY_INLINE(include, IncludeNode, string)
MAKE_NODE_FACTORY_INLINE(version, VersionNode, string)
MAKE_NODE_FACTORY_INLINE(const_decl, ConstDeclNode, ASTNodePtr, string, ASTNodePtr)
MAKE_NODE_FACTORY_INLINE(qubit_decl, QubitDeclNode, string, ASTNodePtr)
MAKE_NODE_FACTORY_INLINE(index_expr, IndexNode, ASTNodePtr, vector<ASTNodeList>)
MAKE_NODE_FACTORY_INLINE(call_expr, CallExprNode, ASTNodePtr, ASTNodeList)
MAKE_NODE_FACTORY_INLINE(measure_stmt, MeasureStmtNode, ASTNodePtr, ASTNodePtr)
