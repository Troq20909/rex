import re

content = open('rex_runtime.zig').read()
original = content

# Fix 1: ArrayListUnmanaged
content = content.replace(
    'nodes: std.ArrayList(usize),',
    'nodes: std.ArrayListUnmanaged(usize),'
)
content = content.replace(
    'edges: std.ArrayList([2]usize),',
    'edges: std.ArrayListUnmanaged([2]usize),'
)
content = content.replace(
    '.nodes = std.ArrayList(usize).init(allocator),',
    '.nodes = std.ArrayListUnmanaged(usize){},'
)
content = content.replace(
    '.edges = std.ArrayList([2]usize).init(allocator),',
    '.edges = std.ArrayListUnmanaged([2]usize){},'
)
content = content.replace(
    'self.nodes.deinit();',
    'self.nodes.deinit(self.allocator);'
)
content = content.replace(
    'self.edges.deinit();',
    'self.edges.deinit(self.allocator);'
)
content = content.replace(
    'try self.nodes.append(id);',
    'try self.nodes.append(self.allocator, id);'
)
content = content.replace(
    'try self.edges.append([2]usize{ u, v });',
    'try self.edges.append(self.allocator, [2]usize{ u, v });'
)

# Fix 2: math.min -> @min
content = content.replace('math.min(', '@min(')

# Fix 3: math.fabs -> @abs
content = content.replace('math.fabs(', '@abs(')

# Fix 4: spectrum_comptime - fix signature to accept f64 not *const UnGraph
old_spectrum = '''pub fn spectrum_comptime(graph: *const UnGraph) f64 {
    const n = @as(f64, @floatFromInt(graph.nodeCount()));
    const m = @as(f64, @floatFromInt(graph.edgeCount()));
    if (n <= 0.0) return 0.0;
    return math.sqrt(2.0 * m / n);
}'''
new_spectrum = '''pub fn spectrum_comptime(val: f64) f64 {
    if (val <= 0.0) return 0.0;
    return math.sqrt(val);
}'''
content = content.replace(old_spectrum, new_spectrum)

# Fix 5: export_to_csv - new file API
old_csv = '''pub fn export_to_csv(label: []const u8, filename: []const u8) void {
    const file = std.fs.cwd().createFile(filename, .{}) catch {
        std.debug.print("Warning: could not create {s}\\n", .{filename});
        return;
    };
    defer file.close();
    const writer = file.writer();
    writer.print("label,value\\n{s},0.0\\n", .{label}) catch return;
    std.debug.print("Exported results to {s}\\n", .{filename});
}'''
new_csv = '''pub fn export_to_csv(label: []const u8, filename: []const u8) void {
    var buf: [256]u8 = undefined;
    const content = std.fmt.bufPrint(&buf, "label,value\\n{s},0.0\\n", .{label}) catch return;
    std.fs.cwd().writeFile(.{ .sub_path = filename, .data = content }) catch {
        std.debug.print("Warning: could not write {s}\\n", .{filename});
        return;
    };
    std.debug.print("Exported results to {s}\\n", .{filename});
}'''
content = content.replace(old_csv, new_csv)

open('rex_runtime.zig', 'w').write(content)

changes = sum([
    'ArrayListUnmanaged' in content,
    '@min(' in content,
    '@abs(' in content,
    'spectrum_comptime(val: f64)' in content,
    'bufPrint' in content,
])
print(f'Done. {changes}/5 fixes applied.')
if content == original:
    print('WARNING: No changes made - check if already patched or patterns differ.')
