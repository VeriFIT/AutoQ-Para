OPENQASM 3.0;
include "stdgates.inc";

const uint n = __nondet();

// qubits
qubit[n] q;
qubit[n] secret;
qubit ancilla;

// pre: q == 0...00, ancilla == 0, secret == ?...??

// state preparation
h q;
x ancilla;
h ancilla;

// oracle
ccx q, secret, ancilla;

// ALTERNATIVE ORACLE DEFINITION
//
// for uint i in [0:n-1] {
// 	ccx q[i], secret[i], ancilla;
// }


// post-processing
h q;
h ancilla;

// post: q == secret, ancilla == 1
