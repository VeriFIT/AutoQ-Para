OPENQASM 3.0;
include "stdgates.inc";
const uint n=*;
qubit[n] q;
qubit[n] secret;
qubit ancilla;

h q;
h ancilla;
ccx q, secret, ancilla;  // oracle(q, secret, ancilla);
h q;
h ancilla;
