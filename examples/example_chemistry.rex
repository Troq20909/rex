// example_chemistry.rex
// Molecular energy, reaction rates, quantum chemistry Hamiltonians
// Run: python rex_transpiler.py examples/example_chemistry.rex

lab_ready;

molecule_energy benzene using huckel;
molecule_energy h2o using sto3g;
reaction_rate arrhenius;
quantum_chemistry_hamiltonian h2 basis sto3g;
