// HXH on all qubits: should be equivalent to Z on all qubits
OPENQASM 3.0;
include "stdgates.inc";

const uint n = __nondet();

// qubits
qubit[n] q;

h q;
x q;
h q;
