# Synchronized Weighted Tree Automata (SWTAs) implementation

The prototype written in Python in `python_prototype/` is capable 
of checking whether a given WA is zero, however, it does not support
algebraic representation of complex numbers. All further work is focused
on the C++ implementation located in `src/`.

## C++ implementation - tasks:

### Test suite
- [x] ACNs - test scaling during addition
- [ ] ACN Matrix - test row insertion into row-echelon matrix

### Features
- [x] Algebraic representation of complex numbers (ACN)
- [x] ACN Matrices
- [ ] Precise display of ACNs (for debugging)
- [ ] WA representation
- [ ] Checking whether a given WA is zero
- [ ] SWTA representation - probably using (some sort of) linear forms
- [ ] SWTA>WA conversion
- [ ] Transformer representation - using linear forms and transformers
- [ ] Applying transformers to automata
- [ ] Transformer composition
- [ ] BV (almost) end-to-end verification - with given transformers

# Compiling transformers from QASM
## List of patterns to support:
... 
