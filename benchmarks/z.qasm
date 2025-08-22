// Z on all qubits
OPENQASM 3.0;
include "stdgates.inc";

const uint n = __nondet();

// qubits
qubit[n] q;

z q;
