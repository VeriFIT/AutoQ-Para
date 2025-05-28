// implements a parameterized version of the ripple-carry adder from https://arxiv.org/abs/quant-ph/0410184, using
// the circuit from https://github.com/openqasm/openqasm/blob/main/examples/adder.qasm for inspiration

OPENQASM 3.0;
include "stdgates.inc";

// constants declarations
const uint n=__nondet_uint();

// auxiliary gates
gate maj cin, b, a { // MAJority
    cx b, a;
    cx b, cin;
    ccx cin, b, a;
}

gate uma cin, b, a { // UnMajority and Add
    ccx cin, b, a;
    cx a, cin;
    cx cin, b;
}

qubit[n] a;     // first operand
qubit[n] b;     // second operand (and result in the end)
qubit c0;       // initial carry
qubit z;        // ancilla

// pre: a, b: operands, c0 = z = 0

maj c0, b[0], a[0];  // initialization

// first staircase (downwards)
for uint i in [0:n-2] {
	maj	a[i], b[i+1], a[i+1];
}

cnot a[n-1], z; // working with ancilla

// second staircase (upwards)
for uint i in [n-2:-1:0] {
	uma	a[i], b[i+1], a[i+1];
}

uma c0, b[0], a[0];  // finalization

// post: a_out = a_in, b_out = a_in + b_in, c0 = 0, z = z_in \xor b_out[n-1]
