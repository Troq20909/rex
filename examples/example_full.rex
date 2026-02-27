// example_full.rex
// Full multi-domain demo – featured in README
// Covers every major Rex domain in one model
// Run: python rex_transpiler.py examples/example_full.rex

lab_ready;

// --- User-defined types and functions ---

struct Experiment {
    n_sites: f64;
    coupling: f64;
    temperature: f64;
}

fn compute_energy(m: f64, v: f64) -> f64 {
    let ke = 0.5 * m * v * v;
}

// --- Graph construction ---

graph lattice = (undirected, [(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)]);

// --- Spectral graph analysis ---

laplacian_spectrum lattice;
fiedler_value lattice;
eigenvalue_stats lattice;

// --- Quantum walks and chaos ---

magnetic_quantum_walk ctqw on lattice time 20 flux 0.314;
loschmidt_echo lattice;
participation_ratio lattice eigenstates;
spacing_analysis lattice;
chaos_stats lattice;

// --- Random matrix theory ---

rmt goe size 50;
spectrum sp of goe_rmt;

// --- Complex analysis ---

let pole = 2;
residue r at pole;
contour_integral ci radius 3;

// --- Chemistry ---

molecule_energy h2o using huckel;
reaction_rate arrhenius;
quantum_chemistry_hamiltonian h2 basis sto3g;

// --- Genomics ---

dna_sequence seq1 = "GCTAGCTAGCTAGCTAGCTA";
crispr_offtarget seq1 mismatches 1 pam "NGG";
gene_expression tumor;
population_genetics cohort size 500;
vector_titer aav9;
delivery_efficiency liver compartments 3;

// --- Thermal and fluid dynamics ---

heat_equation_1d bar length 2.0 dt 0.01 steps 50;
thermo_efficiency carnot hot 900 cold 300;
phase_transition water temp 373;
burgers_1d wave velocity 0.8 viscosity 0.02;
navier_stokes_2d_simple cavity viscosity 0.001 steps 200;
turbulence_stats cavity;

// --- Lab deployment ---

let confidence = 0.95;
validate confidence;
export_csv all_results to "full_results.csv";
lab_daq_output "instrument.bin" to stdout;
