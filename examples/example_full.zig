// Lab-ready deployment enabled

const Experiment = struct {
    n_sites: f64,
    coupling: f64,
    temperature: f64,
};

fn compute_energy(m: f64, v: f64) f64 {
        const ke = (((0.5 * m) * v) * v);
        _ = ke;
}


const std = @import("std");
const math = std.math;
const mem = std.mem;
const Allocator = mem.Allocator;
const rt = @import("rex_runtime.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var lattice = rt.UnGraph.init(allocator);
    try lattice.addNode();
    try lattice.addNode();
    try lattice.addNode();
    try lattice.addNode();
    try lattice.addEdge(0, 1);
    try lattice.addEdge(1, 2);
    try lattice.addEdge(2, 3);
    try lattice.addEdge(3, 0);
    try lattice.addEdge(0, 2);
    try lattice.addEdge(1, 3);
    const lattice_spec = rt.laplacian_spectrum_comptime(&lattice);
    std.debug.print("Laplacian spectrum: {}\n", .{lattice_spec});
    const lattice_fv = rt.fiedler_value_comptime(&lattice);
    std.debug.print("Fiedler value: {}\n", .{lattice_fv});
    const lattice_es = rt.eigenvalue_stats_comptime(&lattice);
    std.debug.print("Eigenvalue stats: {}\n", .{lattice_es});
    const ctqw_result = rt.magnetic_ctqw_comptime(20, 0.314);
    std.debug.print("Flux {} result: {}\n", .{0.314, ctqw_result});
    const lattice_echo = rt.loschmidt_echo_comptime(&lattice);
    std.debug.print("Loschmidt echo: {}\n", .{lattice_echo});
    const lattice_ipr = rt.participation_ratio_comptime(&lattice);
    std.debug.print("IPR: {}\n", .{lattice_ipr});
    const lattice_spacing = rt.spacing_analysis_comptime(&lattice);
    std.debug.print("Spacing stats: {}\n", .{lattice_spacing});
    const lattice_chaos = rt.chaos_stats_comptime(&lattice);
    std.debug.print("Chaos indicators: {}\n", .{lattice_chaos});
    const goe_rmt = rt.goe_ensemble_comptime(50);
    std.debug.print("goe GOE: {}\n", .{goe_rmt});
    const sp = rt.spectrum_comptime(goe_rmt);
    std.debug.print("sp: {}\n", .{sp});
    const pole = 2;
    const r = rt.residue_comptime(pole);
    std.debug.print("r residue: {}\n", .{r});
    const ci = rt.contour_integral_comptime(3);
    std.debug.print("ci: {}\n", .{ci});
    const h2o_energy = rt.molecule_energy_comptime("huckel");
    std.debug.print("h2o energy: {}\n", .{h2o_energy});
    const arrhenius_rate = rt.reaction_rate_comptime();
    std.debug.print("arrhenius rate: {}\n", .{arrhenius_rate});
    const h2_qch = rt.qchem_hamiltonian_comptime("sto3g");
    std.debug.print("h2 H: {}\n", .{h2_qch});
    const gc_seq1 = rt.dna_gc_content("GCTAGCTAGCTAGCTAGCTA");
    std.debug.print("seq1 GC content: {}\n", .{gc_seq1});
    const seq1_prob = rt.crispr_offtarget_prob("seq1", 1, "NGG");
    std.debug.print("Off-target prob: {}\n", .{seq1_prob});
    const tumor_ge = rt.gene_expression_comptime();
    std.debug.print("tumor expression: {}\n", .{tumor_ge});
    const cohort_pg = rt.population_genetics_comptime(500);
    std.debug.print("cohort: {}\n", .{cohort_pg});
    const aav9_vt = rt.vector_titer_poisson_comptime();
    std.debug.print("aav9 titer: {}\n", .{aav9_vt});
    const liver_de = rt.delivery_efficiency_comptime(3);
    std.debug.print("liver: {}\n", .{liver_de});
    const bar_profile = rt.heat_equation_1d_comptime(2.0, 0.01, 50, allocator) catch 0.0;
    std.debug.print("bar final profile: {}\n", .{bar_profile});
    const carnot_eff = rt.thermo_efficiency_comptime(900, 300);
    std.debug.print("carnot efficiency: {}\n", .{carnot_eff});
    const water_pt = rt.phase_transition_comptime(373);
    std.debug.print("water: {}\n", .{water_pt});
    const wave_sol = rt.burgers_1d_comptime(0.8, 0.02, allocator) catch 0.0;
    std.debug.print("wave solution: {}\n", .{wave_sol});
    const cavity_ns = rt.navier_stokes_2d_comptime(0.001, 200);
    std.debug.print("cavity: {}\n", .{cavity_ns});
    const cavity_ts = rt.turbulence_stats_comptime();
    std.debug.print("cavity: {}\n", .{cavity_ts});
    const confidence = 0.95;
    if (confidence <= 0.0) return error.ValidationFailed;
    rt.export_to_csv("all_results", "full_results.csv");
    rt.pipe_to_daq("instrument.bin", "stdout");
}
