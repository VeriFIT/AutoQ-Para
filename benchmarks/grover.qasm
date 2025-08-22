// implements a parameterized version of Grover's algorithm with the oracle "all zeros"

OPENQASM 3.0;
include "stdgates.inc";

// constants declarations
// const uint n=__nondet_uint();
const uint n=5;

// static_assert: n > 2

// the number of iterations
const uint iter=floor(pi/4 * sqrt(2^n));

/////////////////////////////////////////////////////////////////////
// performs a multi-control Toffoli (CC...CX) gate using ancillas in the Oracle
def mct_oracle(qubit[n] control, qubit[n-1] anc, qubit target) {
	ccx control[0], control[1], anc[0];
	for uint i in [2:n-1] {
		ccx anc[i-2], control[i], anc[i-1];
	}

	cx anc[n-2], target;

	for uint i in [n-1:-1:2] {
		ccx anc[i-2], control[i], anc[i-1];
	}
	ccx control[0], control[1], anc[0];
}
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// performs a multi-control Toffoli (CC...CX) gate using ancillas in the Diffuser
def mct_diff(qubit[n] control, qubit[n-1] anc, qubit target) {
	ccx control[0], control[1], anc[0];
	for uint i in [2:n-2] {
		ccx anc[i-2], control[i], anc[i-1];
	}

	cx anc[n-3], target;

	for uint i in [n-2:-1:2] {
		ccx anc[i-2], control[i], anc[i-1];
	}
	ccx control[0], control[1], anc[0];
}
/////////////////////////////////////////////////////////////////////


qubit[n] w;          // working tape
qubit a;             // ancilla for phase oracle
qubit[n-1] mct_anc;  // ancillas for multi-control Toffolis

// start by Hadamarding all qubits
h w;

// put ancilla into |->
x a;
h a;

// the loop
for uint i in [0:n-2] {
	// phase oracle
	mct_diff(w, mct_anc, a);

	// Grover's diffuser
	h w;
	x w;

	// multi-control Z
	h w[n-1];
	mct_diff(w[0:n-2], mct_anc[0:n-3], w[n-1]);
	h w[n-1];

	x w;
	h w;
}

// post: TODO
