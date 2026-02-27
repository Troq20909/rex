// build.zig
// Rex build file – Zig 0.15.2
// Usage:
//   zig build run                          – debug build
//   zig build run -Doptimize=ReleaseFast   – optimized build (recommended)
//   zig build run -Dmodel=examples/example_full  – run a specific model
//   zig build strip                        – smallest binary for lab deployment

const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    // Allow the user to specify which model to run (default: model)
    const model_name = b.option(
        []const u8,
        "model",
        "Name of the transpiled model file (without .zig extension, default: model)",
    ) orelse "model";

    const src_file = b.fmt("{s}.zig", .{model_name});

    // Rex runtime module – shared by all models
    const runtime_mod = b.createModule(.{
        .root_source_file = b.path("rex_runtime.zig"),
    });

    // Main executable
    const exe = b.addExecutable(.{
        .name = "rex_model",
        .root_source_file = b.path(src_file),
        .target = target,
        .optimize = optimize,
    });
    exe.root_module.addImport("rex_runtime.zig", runtime_mod);
    b.installArtifact(exe);

    // Run step
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    const run_step = b.step("run", "Build and run the Rex model");
    run_step.dependOn(&run_cmd.step);

    // Strip step for lab/embedded deployment
    const strip_exe = b.addExecutable(.{
        .name = "rex_model_stripped",
        .root_source_file = b.path(src_file),
        .target = target,
        .optimize = .ReleaseFast,
    });
    strip_exe.root_module.addImport("rex_runtime.zig", runtime_mod);
    strip_exe.root_module.strip = true;

    const strip_step = b.step("strip", "Build a stripped ReleaseFast binary for lab deployment");
    const strip_install = b.addInstallArtifact(strip_exe, .{});
    strip_step.dependOn(&strip_install.step);
}
