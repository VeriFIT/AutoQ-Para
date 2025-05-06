OPENQASM 3.0;
include "stdgates.inc";

const uint n = __nondet();

// qubits
qubit[n] q;
qubit ancilla;

// pre: q == 0...00, ancilla == 0

// state preparation
h q;
x ancilla;
h ancilla;

for uint i in [0:n-1] {
	if (i % 2 = 0) {
		cx q[i], ancilla;
	}
	else { }
}

// post-processing
h q;
h ancilla;

// post: q == (10)*1?, ancilla == 1
