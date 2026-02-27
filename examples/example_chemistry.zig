// Lab-ready deployment enabled


const std = @import("std");
const math = std.math;
const mem = std.mem;
const Allocator = mem.Allocator;
const rt = @import("rex_runtime.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    _ = allocator;
    const benzene_energy = rt.molecule_energy_comptime("huckel");
    std.debug.print("benzene energy: {}\n", .{benzene_energy});
    const h2o_energy = rt.molecule_energy_comptime("sto3g");
    std.debug.print("h2o energy: {}\n", .{h2o_energy});
    const arrhenius_rate = rt.reaction_rate_comptime();
    std.debug.print("arrhenius rate: {}\n", .{arrhenius_rate});
    const h2_qch = rt.qchem_hamiltonian_comptime("sto3g");
    std.debug.print("h2 H: {}\n", .{h2_qch});
}
