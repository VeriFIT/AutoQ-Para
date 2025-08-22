OPENQASM 3.0;
include "stdgates.inc";

const uint n = __nondet();

// qubits
qubit[n] q;

// pre: q == 0...00

h q[0];

for uint i in [0:n-2] {
	cx q[i], q[i+1];
}

// post: q == 1/sqrt(2) * (|0...00> + |1...11>)
