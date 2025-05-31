// implements a parameterized version of QFT (Quantum Fourier Transform), cf.
// https://www.ryanlarose.com/uploads/1/1/5/8/115879647/shor2-qft.pdf
// inspired by https://github.com/openqasm/openqasm/blob/main/examples/qft.qasm

OPENQASM 3.0;
include "stdgates.inc";

const uint n = __nondet();

// OL: not sure whether it is correct

// qubits
qubit[n] q;

// pre: ??????

// the core
for uint i in [0:n-1] {
	h q[i];
	for uint j in [i+1:n-1] {
		angle rot = pi/(2^{j-i});
		cphase(rot) q[i], q[j];
	}
}

// final swaps
for uint i in [0:n/2] {
	swap q[i], q[n-i-1];
}

// post: ?????????
