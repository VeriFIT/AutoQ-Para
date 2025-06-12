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
    For,
    Include,
    Version,
    ConstDecl,
    QubitDecl,
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
        cout << string(indent, ' ') << "Range start:end:step" << endl;
        if (start)
            start->print(indent + 2);
        if (end)
            end->print(indent + 2);
        if (step)
            step->print(indent + 2);
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
        cout << string(indent, ' ') << "Block" << endl;
        for (const auto &stmt : statements)
            stmt->print(indent + 2);
    }
};

using GateOperand = pair<string, optional<vector<vector<shared_ptr<ASTNode>>>>>;

struct GateCallNode : ASTNode
{
    string name;
    vector<GateOperand> operand_list;
    GateCallNode(string g, vector<GateOperand> operand_list)
        : name(g), operand_list(operand_list) {}
    ASTKind kind() const override { return ASTKind::GateCall; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "GateCall(" << name << ", ";
        for (const auto &arg : operand_list)
        {
            cout << arg.first;
            if (!arg.second.has_value() || arg.second.value().empty() || arg.second.value().front().empty())
                cout << "#NULL";
            else
            {
                for (const auto &exps : arg.second.value())
                {
                    cout << "[";
                    for (size_t k = 0; k < exps.size(); k++)
                    {
                        auto e1 = dynamic_pointer_cast<CanShow>(exps[k]);
                        cout << (e1 ? (e1->toString() + (k + 1 < exps.size() ? ", " : "")) : "Expr( ?... )");
                    }
                    cout << "]";
                }
            }
            cout << " ";
        }
        cout << ")" << endl;
    }
};

struct ForNode : ASTNode
{
    ASTNodePtr iter;
    ASTNodePtr range;
    ASTNodePtr body;
    ForNode(ASTNodePtr i, ASTNodePtr r, ASTNodePtr b)
        : iter(i), range(r), body(b) {}
    ASTKind kind() const override { return ASTKind::For; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "For" << endl;
        if (iter)
            iter->print(indent + 2);
        if (range)
            range->print(indent + 2);
        if (body)
            body->print(indent + 2);
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
    string type;
    string name;
    string expr;
    ConstDeclNode(string t, string n, string e)
        : type(t), name(n), expr(e) {}
    ASTKind kind() const override { return ASTKind::ConstDecl; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "ConstDecl(" << type << " " << name << " = " << expr << ")" << endl;
    }
};

struct QubitDeclNode : ASTNode
{
    string name;
    string size;
    QubitDeclNode(string n, string s = "") : name(n), size(s) {}
    ASTKind kind() const override { return ASTKind::QubitDecl; }
    void print(int indent = 0) const override
    {
        cout << string(indent, ' ') << "QubitDecl(" << name;
        if (!size.empty())
            cout << " (" << size << ")";
        cout << ")\n";
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
auto gatecall(const string &gate, Operands &&...ops)
{
    vector<GateOperand> v;
    (v.push_back(forward<Operands>(ops)), ...);
    return make_shared<GateCallNode>(gate, move(v));
}

inline ASTNodeList astn_vec(const ASTNodeList &v)
{
    return v;
}

inline GateOperand gateoperand_vec(const string &key, const vector<vector<ASTNodePtr>> &v)
{
    return {key, v.empty() ? nullopt : make_optional(v)};
}

inline shared_ptr<GateCallNode> gatecall_vec(const string &gate, const vector<GateOperand> &gq_list)
{
    return make_shared<GateCallNode>(gate, gq_list);
}

template <typename _, typename __, typename ___>
auto fnode(_ &&arg1, __ &&arg2, ___ &&arg3)
{
    return make_shared<ForNode>(arg1, arg2, arg3);
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
MAKE_NODE_FACTORY_INLINE(const_decl, ConstDeclNode, string, string, string)
MAKE_NODE_FACTORY_INLINE(qubit_decl, QubitDeclNode, string, string)