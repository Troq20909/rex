# REX

A lightweight domain-specific language for high-performance scientific research modeling.
Transpiles `.rex` source files to native Zig code — producing small (<85KB), fast binaries
with zero runtime overhead.

---

## Features

- **7 scientific domains** — quantum chaos, spectral graph theory, random matrix theory, chemistry, genomics, fluid dynamics, thermal physics
- **Native Zig output** — compiles to small, fast binaries with no garbage collector
- **Lab-ready** — CSV export, DAQ piping, embedded no_std support
- **Zero dependencies** — transpiler is pure Python, runtime is a single Zig file

---

## Quick Start
```bash
# Transpile a Rex model
python rex_transpiler.py examples/example_chemistry.rex

# Compile and run (requires Zig 0.15.2)
zig run examples/example_chemistry.zig
```

---

## Domain Reference

| Keyword | Domain | What it computes |
|---|---|---|
| `magnetic_quantum_walk` | Quantum chaos | CTQW return probability with magnetic flux |
| `loschmidt_echo` | Quantum chaos | Gaussian fidelity decay |
| `participation_ratio` | Quantum chaos | Inverse participation ratio |
| `spacing_analysis` | Quantum chaos | GOE level spacing ratio |
| `chaos_stats` | Quantum chaos | Brody parameter |
| `laplacian_spectrum` | Spectral graphs | Trace of Laplacian |
| `fiedler_value` | Spectral graphs | Algebraic connectivity |
| `eigenvalue_stats` | Spectral graphs | Mean level spacing |
| `quantum_walk` | Spectral graphs | Bessel approximation on cycle |
| `rmt` | Random matrices | GOE Wigner surmise spacing |
| `spectrum` | Random matrices | Spectral radius estimate |
| `residue` | Complex analysis | Residue at simple pole |
| `contour_integral` | Complex analysis | Contour integral value |
| `molecule_energy` | Chemistry | Hückel MO pi energy |
| `reaction_rate` | Chemistry | Arrhenius rate at 298K |
| `quantum_chemistry_hamiltonian` | Chemistry | H2 STO-3G ground state |
| `dna_sequence` | Genomics | GC content |
| `crispr_offtarget` | Genomics | Off-target probability |
| `gene_expression` | Genomics | Fold change |
| `population_genetics` | Genomics | Neutral drift fixation |
| `vector_titer` | Gene therapy | AAV9 Poisson titer |
| `delivery_efficiency` | Gene therapy | Multi-compartment efficiency |
| `heat_equation_1d` | Fluid dynamics | FTCS finite difference |
| `burgers_1d` | Fluid dynamics | Upwind scheme L2 norm |
| `navier_stokes_2d_simple` | Fluid dynamics | Kinetic energy decay |
| `turbulence_stats` | Fluid dynamics | Kolmogorov -5/3 law |
| `thermo_efficiency` | Thermodynamics | Carnot efficiency |
| `phase_transition` | Thermodynamics | Mean-field order parameter |
| `validate` | Lab | Validation check |
| `export_csv` | Lab | CSV file export |
| `lab_daq_output` | Lab | DAQ instrument piping |
| `embedded_no_std` | Lab | no_std hardened binary |

---

## Language Reference
```rex
// Graphs
graph g = (undirected, [(0,1),(1,2),(2,0)]);
laplacian_spectrum g;
fiedler_value g;

// Quantum walks
magnetic_quantum_walk ctqw on g time 10 flux 0.314;
loschmidt_echo g;

// Random matrices
rmt ensemble size 100;

// Chemistry
molecule_energy benzene using huckel;
reaction_rate arrhenius;

// Genomics
dna_sequence brca1 = "ATGCGT";
crispr_offtarget brca1 mismatches 2 pam "NGG";

// Fluids
heat_equation_1d rod length 1.0 dt 0.01 steps 100;
burgers_1d shock velocity 1.0 viscosity 0.01;

// Lab output
export_csv all_results to "results.csv";
lab_daq_output "data.bin" to stdout;
```

---

## Installation

**Requirements:** Python 3.8+, Zig 0.15.2
```bash
git clone https://github.com/Troq20909/rex.git
cd rex
python rex_transpiler.py examples/example_full.rex
zig run examples/example_full.zig
```

---

## Project Structure
```
rex/
├── rex_transpiler.py   # Rex → Zig transpiler (pure Python)
├── rex_runtime.zig     # Runtime library (35 functions)
├── build.zig           # Zig build system integration
├── CONTRIBUTING.md     # How to add new keywords
├── LICENSE             # MIT
└── examples/
    ├── example_chemistry.rex
    ├── example_fluids.rex
    ├── example_full.rex
    ├── example_genomics.rex
    ├── example_lab.rex
    ├── example_quantum.rex
    └── example_rmt.rex
```

---

## Known Limitations

- Function bodies support only `let` and `print` statements
- Domain keywords must appear at top level
- Comments are single-line only (`//`)
- Zig 0.15.2 required

---

## Roadmap

- [ ] Type inference
- [ ] Multi-file models
- [ ] WebAssembly target
- [ ] Jupyter notebook integration
- [ ] More domain keywords

---

## License

MIT — use, modify, distribute freely for research.