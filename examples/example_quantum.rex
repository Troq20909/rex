// example_quantum.rex
// Quantum chaos, magnetic walks, Loschmidt echo, spectral graph theory
// Run: python rex_transpiler.py examples/example_quantum.rex

lab_ready;

graph g = (undirected, [(0,1),(1,2),(2,3),(3,0),(0,2)]);

laplacian_spectrum g;
fiedler_value g;
eigenvalue_stats g;

magnetic_quantum_walk ctqw on g time 10 flux 0.5;
loschmidt_echo g;
participation_ratio g eigenstates;
spacing_analysis g;
chaos_stats g;

quantum_walk qw on g time 5;
