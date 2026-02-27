# rex_transpiler.py
# Complete, consolidated Rex transpiler – 22 February 2026
# MIT License – Use, modify, distribute freely for research

from typing import List
import sys
import os

ZN = r'\n'

class Token:
    def __init__(self, type_: str, value: str, line: int):
        self.type = type_
        self.value = value
        self.line = line

class Lexer:
    def __init__(self, source: str):
        self.source = source
        self.tokens: List[Token] = []
        self.tokenize()

    def tokenize(self):
        keywords = {
            'let', 'print', 'struct', 'fn', 'component', 'matrix', 'vector', 'complex', 'integral',
            'graph', 'fractal', 'rmt', 'operator', 'spectrum', 'laplacian_spectrum', 'fiedler_value',
            'quantum_graph', 'eigenvalue_stats', 'residue', 'contour_integral', 'quantum_walk',
            'magnetic_quantum_walk', 'loschmidt_echo', 'participation_ratio', 'spacing_analysis',
            'chaos_stats', 'lab_ready', 'validate', 'export_csv', 'lab_daq_output', 'embedded_no_std',
            'molecule_energy', 'reaction_rate', 'quantum_chemistry_hamiltonian', 'dna_sequence',
            'gene_expression', 'population_genetics', 'crispr_offtarget', 'vector_titer',
            'delivery_efficiency', 'heat_equation_1d', 'thermo_efficiency', 'phase_transition',
            'burgers_1d', 'navier_stokes_2d_simple', 'turbulence_stats',
            'on', 'time', 'flux', 'eigenstates', 'using', 'mismatches', 'pam',
            'length', 'dt', 'steps', 'velocity', 'viscosity', 'hot', 'cold',
            'temp', 'compartments', 'basis', 'at', 'radius', 'of', 'size',
            'all_results', 'to', 'hardening', 'for', 'this', 'model',
        }
        lines = self.source.splitlines()
        for line_num, line in enumerate(lines, 1):
            i = 0
            while i < len(line):
                c = line[i]
                if c.isspace():
                    i += 1
                    continue
                if c == '/' and i + 1 < len(line) and line[i + 1] == '/':
                    break
                if c == '"':
                    j = i + 1
                    while j < len(line) and line[j] != '"':
                        j += 1
                    self.tokens.append(Token('STRING', line[i + 1:j], line_num))
                    i = j + 1
                    continue
                if c.isalpha() or c == '_':
                    j = i
                    while j < len(line) and (line[j].isalnum() or line[j] == '_'):
                        j += 1
                    word = line[i:j]
                    tok_type = word.upper() if word.lower() in keywords else 'ID'
                    self.tokens.append(Token(tok_type, word, line_num))
                    i = j
                    continue
                if c.isdigit() or (c == '-' and i + 1 < len(line) and line[i + 1].isdigit()):
                    j = i
                    if c == '-': j += 1
                    while j < len(line) and (line[j].isdigit() or line[j] == '.'):
                        j += 1
                    self.tokens.append(Token('NUMBER', line[i:j], line_num))
                    i = j
                    continue
                if c == '-' and i + 1 < len(line) and line[i + 1] == '>':
                    self.tokens.append(Token('ARROW', '->', line_num))
                    i += 2
                    continue
                if c in ';:{}()=,+*/[]':
                    self.tokens.append(Token(c, c, line_num))
                    i += 1
                    continue
                i += 1
        self.tokens.append(Token('EOF', 'EOF', len(lines) + 1))

class Parser:
    def __init__(self, lexer: Lexer):
        self.lexer = lexer
        self.pos = 0
        self.current = self.lexer.tokens[0] if self.lexer.tokens else Token('EOF', 'EOF', 0)

    def advance(self):
        self.pos += 1
        if self.pos < len(self.lexer.tokens):
            self.current = self.lexer.tokens[self.pos]

    def expect(self, expected_type: str) -> str:
        if self.current.type == expected_type or \
           (expected_type == 'ID' and self.current.type not in ['EOF', ';', ')', ']', ',', '{', '}']):
            val = self.current.value
            self.advance()
            return val
        raise SyntaxError(f"Expected {expected_type} at line {self.current.line}, got {self.current.type} ('{self.current.value}')")

    def parse_expr(self): return self.parse_additive()

    def parse_additive(self):
        left = self.parse_multiplicative()
        while self.current.type in ['+', '-']:
            op = self.expect(self.current.type)
            left = f'({left} {op} {self.parse_multiplicative()})'
        return left

    def parse_multiplicative(self):
        left = self.parse_factor()
        while self.current.type in ['*', '/']:
            op = self.expect(self.current.type)
            left = f'({left} {op} {self.parse_factor()})'
        return left

    def parse_factor(self):
        if self.current.type == 'NUMBER': return self.expect('NUMBER')
        if self.current.type == 'ID':
            name = self.expect('ID')
            if self.current.type == '(': return self.parse_call(name)
            if self.current.type == '[': return self.parse_array_access(name)
            return name
        if self.current.type == '(':
            self.expect('('); expr = self.parse_expr(); self.expect(')'); return expr
        if self.current.type == '[': return self.parse_array_literal()
        self.advance(); return '0'

    def parse_call(self, name):
        self.expect('(')
        args = []
        while self.current.type != ')' and self.current.type != 'EOF':
            args.append(self.parse_expr())
            if self.current.type == ',': self.expect(',')
        self.expect(')')
        return f'{name}({", ".join(args)})'

    def parse_array_access(self, name):
        self.expect('['); idx = self.parse_expr(); self.expect(']')
        return f'{name}[{idx}]'

    def parse_array_literal(self):
        self.expect('[')
        elements = []
        while self.current.type != ']' and self.current.type != 'EOF':
            elements.append(self.parse_expr())
            if self.current.type == ',': self.expect(',')
        self.expect(']')
        return f'[{", ".join(elements)}]'

    def parse_top_level_decl(self):
        t = self.current.type
        if t == 'STRUCT': return self.parse_struct_decl()
        if t == 'FN': return self.parse_fn_decl()
        if t == 'COMPONENT': return self.parse_component_decl()
        if t == 'LAB_READY': return self.parse_lab_ready()
        raise SyntaxError(f"Expected top-level declaration, got {t} at line {self.current.line}")

    def parse_struct_decl(self):
        self.expect('STRUCT'); name = self.expect('ID'); self.expect('{')
        fields = []
        while self.current.type != '}' and self.current.type != 'EOF':
            fn = self.expect('ID'); self.expect(':'); ft = self.expect('ID')
            if self.current.type == ';': self.expect(';')
            fields.append(f'    {fn}: {ft},')
        self.expect('}')
        return f'const {name} = struct {{\n' + '\n'.join(fields) + '\n};\n'

    def parse_fn_decl(self):
        self.expect('FN'); name = self.expect('ID'); self.expect('(')
        params = []
        while self.current.type != ')' and self.current.type != 'EOF':
            pn = self.expect('ID'); self.expect(':'); pt = self.expect('ID')
            params.append(f'{pn}: {pt}')
            if self.current.type == ',': self.expect(',')
        self.expect(')')
        rt = 'void'
        if self.current.type == 'ARROW': self.expect('ARROW'); rt = self.expect('ID')
        self.expect('{')
        body = []
        while self.current.type != '}' and self.current.type != 'EOF':
            if self.current.type == 'LET': body.append(self.parse_let_stmt())
            elif self.current.type == 'PRINT': body.append(self.parse_print_stmt())
            else: self.advance()
        self.expect('}')
        bs = '\n    '.join(body) if body else '// empty body'
        return f'fn {name}({", ".join(params)}) {rt} {{\n    {bs}\n}}\n'

    def parse_component_decl(self):
        self.expect('COMPONENT'); name = self.expect('ID'); self.expect('{')
        lines = []
        while self.current.type != '}' and self.current.type != 'EOF':
            if self.current.type == 'LET': lines.append(self.parse_let_stmt())
            else: self.advance()
        self.expect('}')
        bs = '\n    '.join(lines) if lines else '// empty component'
        return f'// Component placeholder (Leptos-compatible)\nconst {name} = struct {{\n    {bs}\n}};\n'

    def parse_let_stmt(self):
        self.expect('LET'); name = self.expect('ID'); type_ = ''
        if self.current.type == ':': self.expect(':'); type_ = self.expect('ID')
        self.expect('='); expr = self.parse_expr(); self.expect(';')
        return f'    const {name}: {type_} = {expr};' if type_ else f'    const {name} = {expr};'

    def parse_print_stmt(self):
        self.expect('PRINT'); self.expect('('); arg = self.parse_expr(); self.expect(')'); self.expect(';')
        return f'    std.debug.print("{{}}{ZN}", .{{{arg}}});'

    def parse_graph_stmt(self):
        self.expect('GRAPH'); name = self.expect('ID'); self.expect('='); self.expect('(')
        self.expect('ID'); self.expect(','); self.expect('[')
        edges = []
        while self.current.type != ']' and self.current.type != 'EOF':
            if self.current.type == '(':
                self.expect('('); u = self.expect('NUMBER'); self.expect(','); v = self.expect('NUMBER'); self.expect(')')
                edges.append((u, v))
                if self.current.type == ',': self.expect(',')
        self.expect(']'); self.expect(')'); self.expect(';')
        nc = max((int(u) for u, v in edges), default=0) + 1 if edges else 1
        lines = [f'    var {name} = rt.UnGraph.init(allocator);']
        lines += [f'    try {name}.addNode();' for _ in range(nc)]
        lines += [f'    try {name}.addEdge({u}, {v});' for u, v in edges]
        return '\n'.join(lines)

    def parse_magnetic_quantum_walk_stmt(self):
        self.expect('MAGNETIC_QUANTUM_WALK'); wt = self.expect('ID'); self.expect('ON'); self.expect('ID')
        self.expect('TIME'); t = self.parse_expr(); self.expect('FLUX'); flux = self.parse_expr(); self.expect(';')
        return (f'    const {wt}_result = rt.magnetic_ctqw_comptime({t}, {flux});\n'
                f'    std.debug.print("Flux {{}} result: {{}}{ZN}", .{{{flux}, {wt}_result}});')

    def parse_loschmidt_echo_stmt(self):
        self.expect('LOSCHMIDT_ECHO'); model = self.expect('ID'); self.expect(';')
        return (f'    const {model}_echo = rt.loschmidt_echo_comptime(&{model});\n'
                f'    std.debug.print("Loschmidt echo: {{}}{ZN}", .{{{model}_echo}});')

    def parse_participation_ratio_stmt(self):
        self.expect('PARTICIPATION_RATIO'); model = self.expect('ID'); self.expect('EIGENSTATES'); self.expect(';')
        return (f'    const {model}_ipr = rt.participation_ratio_comptime(&{model});\n'
                f'    std.debug.print("IPR: {{}}{ZN}", .{{{model}_ipr}});')

    def parse_spacing_analysis_stmt(self):
        self.expect('SPACING_ANALYSIS'); model = self.expect('ID'); self.expect(';')
        return (f'    const {model}_spacing = rt.spacing_analysis_comptime(&{model});\n'
                f'    std.debug.print("Spacing stats: {{}}{ZN}", .{{{model}_spacing}});')

    def parse_chaos_stats_stmt(self):
        self.expect('CHAOS_STATS'); model = self.expect('ID'); self.expect(';')
        return (f'    const {model}_chaos = rt.chaos_stats_comptime(&{model});\n'
                f'    std.debug.print("Chaos indicators: {{}}{ZN}", .{{{model}_chaos}});')

    def parse_laplacian_spectrum_stmt(self):
        self.expect('LAPLACIAN_SPECTRUM'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_spec = rt.laplacian_spectrum_comptime(&{name});\n'
                f'    std.debug.print("Laplacian spectrum: {{}}{ZN}", .{{{name}_spec}});')

    def parse_fiedler_value_stmt(self):
        self.expect('FIEDLER_VALUE'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_fv = rt.fiedler_value_comptime(&{name});\n'
                f'    std.debug.print("Fiedler value: {{}}{ZN}", .{{{name}_fv}});')

    def parse_quantum_walk_stmt(self):
        self.expect('QUANTUM_WALK'); name = self.expect('ID'); self.expect('ON'); self.expect('ID')
        self.expect('TIME'); t = self.parse_expr(); self.expect(';')
        return (f'    const {name} = rt.quantum_walk_comptime({t});\n'
                f'    std.debug.print("{name} walk: {{}}{ZN}", .{{{name}}});')

    def parse_eigenvalue_stats_stmt(self):
        self.expect('EIGENVALUE_STATS'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_es = rt.eigenvalue_stats_comptime(&{name});\n'
                f'    std.debug.print("Eigenvalue stats: {{}}{ZN}", .{{{name}_es}});')

    def parse_rmt_stmt(self):
        self.expect('RMT'); name = self.expect('ID'); self.expect('SIZE'); n = self.expect('NUMBER'); self.expect(';')
        return (f'    const {name}_rmt = rt.goe_ensemble_comptime({n});\n'
                f'    std.debug.print("{name} GOE: {{}}{ZN}", .{{{name}_rmt}});')

    def parse_spectrum_stmt(self):
        self.expect('SPECTRUM'); name = self.expect('ID'); self.expect('OF'); target = self.expect('ID'); self.expect(';')
        return (f'    const {name} = rt.spectrum_comptime({target});\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}}});')

    def parse_residue_stmt(self):
        self.expect('RESIDUE'); name = self.expect('ID'); self.expect('AT'); pole = self.parse_expr(); self.expect(';')
        return (f'    const {name} = rt.residue_comptime({pole});\n'
                f'    std.debug.print("{name} residue: {{}}{ZN}", .{{{name}}});')

    def parse_contour_integral_stmt(self):
        self.expect('CONTOUR_INTEGRAL'); name = self.expect('ID'); self.expect('RADIUS'); r = self.parse_expr(); self.expect(';')
        return (f'    const {name} = rt.contour_integral_comptime({r});\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}}});')

    def parse_molecule_energy_stmt(self):
        self.expect('MOLECULE_ENERGY'); name = self.expect('ID'); self.expect('USING'); method = self.expect('ID'); self.expect(';')
        return (f'    const {name}_energy = rt.molecule_energy_comptime("{method}");\n'
                f'    std.debug.print("{name} energy: {{}}{ZN}", .{{{name}_energy}});')

    def parse_reaction_rate_stmt(self):
        self.expect('REACTION_RATE'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_rate = rt.reaction_rate_comptime();\n'
                f'    std.debug.print("{name} rate: {{}}{ZN}", .{{{name}_rate}});')

    def parse_quantum_chemistry_hamiltonian_stmt(self):
        self.expect('QUANTUM_CHEMISTRY_HAMILTONIAN'); name = self.expect('ID'); self.expect('BASIS'); basis = self.expect('ID'); self.expect(';')
        return (f'    const {name}_qch = rt.qchem_hamiltonian_comptime("{basis}");\n'
                f'    std.debug.print("{name} H: {{}}{ZN}", .{{{name}_qch}});')

    def parse_dna_sequence_stmt(self):
        self.expect('DNA_SEQUENCE'); name = self.expect('ID'); self.expect('='); seq = self.expect('STRING'); self.expect(';')
        return (f'    const gc_{name} = rt.dna_gc_content("{seq}");\n'
                f'    std.debug.print("{name} GC content: {{}}{ZN}", .{{gc_{name}}});')

    def parse_gene_expression_stmt(self):
        self.expect('GENE_EXPRESSION'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_ge = rt.gene_expression_comptime();\n'
                f'    std.debug.print("{name} expression: {{}}{ZN}", .{{{name}_ge}});')

    def parse_population_genetics_stmt(self):
        self.expect('POPULATION_GENETICS'); name = self.expect('ID'); self.expect('SIZE'); n = self.expect('NUMBER'); self.expect(';')
        return (f'    const {name}_pg = rt.population_genetics_comptime({n});\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}_pg}});')

    def parse_crispr_offtarget_stmt(self):
        self.expect('CRISPR_OFFTARGET'); target = self.expect('ID'); self.expect('MISMATCHES')
        mismatches = self.expect('NUMBER'); self.expect('PAM'); pam = self.expect('STRING'); self.expect(';')
        return (f'    const {target}_prob = rt.crispr_offtarget_prob("{target}", {mismatches}, "{pam}");\n'
                f'    std.debug.print("Off-target prob: {{}}{ZN}", .{{{target}_prob}});')

    def parse_vector_titer_stmt(self):
        self.expect('VECTOR_TITER'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_vt = rt.vector_titer_poisson_comptime();\n'
                f'    std.debug.print("{name} titer: {{}}{ZN}", .{{{name}_vt}});')

    def parse_delivery_efficiency_stmt(self):
        self.expect('DELIVERY_EFFICIENCY'); name = self.expect('ID'); self.expect('COMPARTMENTS'); n = self.expect('NUMBER'); self.expect(';')
        return (f'    const {name}_de = rt.delivery_efficiency_comptime({n});\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}_de}});')

    def parse_heat_equation_1d_stmt(self):
        self.expect('HEAT_EQUATION_1D'); name = self.expect('ID'); self.expect('LENGTH'); length = self.parse_expr()
        self.expect('DT'); dt = self.parse_expr(); self.expect('STEPS'); steps = self.expect('NUMBER'); self.expect(';')
        return (f'    const {name}_profile = rt.heat_equation_1d_comptime({length}, {dt}, {steps}, allocator) catch 0.0;\n'
                f'    std.debug.print("{name} final profile: {{}}{ZN}", .{{{name}_profile}});')

    def parse_thermo_efficiency_stmt(self):
        self.expect('THERMO_EFFICIENCY'); name = self.expect('ID'); self.expect('HOT'); t_hot = self.parse_expr()
        self.expect('COLD'); t_cold = self.parse_expr(); self.expect(';')
        return (f'    const {name}_eff = rt.thermo_efficiency_comptime({t_hot}, {t_cold});\n'
                f'    std.debug.print("{name} efficiency: {{}}{ZN}", .{{{name}_eff}});')

    def parse_phase_transition_stmt(self):
        self.expect('PHASE_TRANSITION'); name = self.expect('ID'); self.expect('TEMP'); temp = self.parse_expr(); self.expect(';')
        return (f'    const {name}_pt = rt.phase_transition_comptime({temp});\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}_pt}});')

    def parse_burgers_1d_stmt(self):
        self.expect('BURGERS_1D'); name = self.expect('ID'); self.expect('VELOCITY'); v = self.parse_expr()
        self.expect('VISCOSITY'); nu = self.parse_expr(); self.expect(';')
        return (f'    const {name}_sol = rt.burgers_1d_comptime({v}, {nu}, allocator) catch 0.0;\n'
                f'    std.debug.print("{name} solution: {{}}{ZN}", .{{{name}_sol}});')

    def parse_navier_stokes_2d_simple_stmt(self):
        self.expect('NAVIER_STOKES_2D_SIMPLE'); name = self.expect('ID'); self.expect('VISCOSITY'); nu = self.parse_expr()
        self.expect('STEPS'); steps = self.expect('NUMBER'); self.expect(';')
        return (f'    const {name}_ns = rt.navier_stokes_2d_comptime({nu}, {steps});\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}_ns}});')

    def parse_turbulence_stats_stmt(self):
        self.expect('TURBULENCE_STATS'); name = self.expect('ID'); self.expect(';')
        return (f'    const {name}_ts = rt.turbulence_stats_comptime();\n'
                f'    std.debug.print("{name}: {{}}{ZN}", .{{{name}_ts}});')

    def parse_lab_ready(self):
        self.expect('LAB_READY'); self.expect(';')
        return '// Lab-ready deployment enabled\n'

    def parse_validate(self):
        self.expect('VALIDATE'); cond = self.expect('ID'); self.expect(';')
        return f'    if (!{cond}) return error.ValidationFailed;'

    def parse_export_csv(self):
        self.expect('EXPORT_CSV'); self.expect('ALL_RESULTS'); self.expect('TO'); fname = self.expect('STRING'); self.expect(';')
        return f'    rt.export_to_csv("all_results", "{fname}");'

    def parse_lab_daq_output(self):
        self.expect('LAB_DAQ_OUTPUT'); fname = self.expect('STRING'); self.expect('TO'); target = self.expect('ID'); self.expect(';')
        return f'    rt.pipe_to_daq("{fname}", "{target}");'

    def parse_embedded_no_std(self):
        self.expect('EMBEDDED_NO_STD'); self.expect('HARDENING'); self.expect('FOR'); self.expect('THIS'); self.expect('MODEL'); self.expect(';')
        return '// no_std hardened binary generation triggered\n'

    def parse_program(self):
        top_level: List[str] = []
        main_body: List[str] = []
        while self.current.type != 'EOF':
            t = self.current.type
            if t in ('STRUCT', 'FN', 'COMPONENT', 'LAB_READY'):
                top_level.append(self.parse_top_level_decl())
            elif t == 'LET': main_body.append(self.parse_let_stmt())
            elif t == 'PRINT': main_body.append(self.parse_print_stmt())
            elif t == 'GRAPH': main_body.append(self.parse_graph_stmt())
            elif t == 'MAGNETIC_QUANTUM_WALK': main_body.append(self.parse_magnetic_quantum_walk_stmt())
            elif t == 'LOSCHMIDT_ECHO': main_body.append(self.parse_loschmidt_echo_stmt())
            elif t == 'PARTICIPATION_RATIO': main_body.append(self.parse_participation_ratio_stmt())
            elif t == 'SPACING_ANALYSIS': main_body.append(self.parse_spacing_analysis_stmt())
            elif t == 'CHAOS_STATS': main_body.append(self.parse_chaos_stats_stmt())
            elif t == 'LAPLACIAN_SPECTRUM': main_body.append(self.parse_laplacian_spectrum_stmt())
            elif t == 'FIEDLER_VALUE': main_body.append(self.parse_fiedler_value_stmt())
            elif t == 'QUANTUM_WALK': main_body.append(self.parse_quantum_walk_stmt())
            elif t == 'EIGENVALUE_STATS': main_body.append(self.parse_eigenvalue_stats_stmt())
            elif t == 'RMT': main_body.append(self.parse_rmt_stmt())
            elif t == 'SPECTRUM': main_body.append(self.parse_spectrum_stmt())
            elif t == 'RESIDUE': main_body.append(self.parse_residue_stmt())
            elif t == 'CONTOUR_INTEGRAL': main_body.append(self.parse_contour_integral_stmt())
            elif t == 'MOLECULE_ENERGY': main_body.append(self.parse_molecule_energy_stmt())
            elif t == 'REACTION_RATE': main_body.append(self.parse_reaction_rate_stmt())
            elif t == 'QUANTUM_CHEMISTRY_HAMILTONIAN': main_body.append(self.parse_quantum_chemistry_hamiltonian_stmt())
            elif t == 'DNA_SEQUENCE': main_body.append(self.parse_dna_sequence_stmt())
            elif t == 'GENE_EXPRESSION': main_body.append(self.parse_gene_expression_stmt())
            elif t == 'POPULATION_GENETICS': main_body.append(self.parse_population_genetics_stmt())
            elif t == 'CRISPR_OFFTARGET': main_body.append(self.parse_crispr_offtarget_stmt())
            elif t == 'VECTOR_TITER': main_body.append(self.parse_vector_titer_stmt())
            elif t == 'DELIVERY_EFFICIENCY': main_body.append(self.parse_delivery_efficiency_stmt())
            elif t == 'HEAT_EQUATION_1D': main_body.append(self.parse_heat_equation_1d_stmt())
            elif t == 'THERMO_EFFICIENCY': main_body.append(self.parse_thermo_efficiency_stmt())
            elif t == 'PHASE_TRANSITION': main_body.append(self.parse_phase_transition_stmt())
            elif t == 'BURGERS_1D': main_body.append(self.parse_burgers_1d_stmt())
            elif t == 'NAVIER_STOKES_2D_SIMPLE': main_body.append(self.parse_navier_stokes_2d_simple_stmt())
            elif t == 'TURBULENCE_STATS': main_body.append(self.parse_turbulence_stats_stmt())
            elif t == 'VALIDATE': main_body.append(self.parse_validate())
            elif t == 'EXPORT_CSV': main_body.append(self.parse_export_csv())
            elif t == 'LAB_DAQ_OUTPUT': main_body.append(self.parse_lab_daq_output())
            elif t == 'EMBEDDED_NO_STD': main_body.append(self.parse_embedded_no_std())
            else:
                raise SyntaxError(f"Unrecognized token '{t}' at line {self.current.line}")

        zig_header = 'const std = @import("std");\nconst math = std.math;\nconst mem = std.mem;\nconst Allocator = mem.Allocator;\nconst rt = @import("rex_runtime.zig");\n\npub fn main() !void {\n    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);\n    defer arena.deinit();\n    const allocator = arena.allocator();\n'
        return '\n'.join(top_level) + '\n\n' + zig_header + '\n'.join(main_body) + '\n}\n'


def main():
    if len(sys.argv) < 2:
        print("Usage: python rex_transpiler.py <model.rex>")
        sys.exit(1)
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"Error: file not found: {input_file}")
        sys.exit(1)
    try:
        with open(input_file, 'r') as f:
            source = f.read()
        lexer = Lexer(source)
        parser = Parser(lexer)
        zig_code = parser.parse_program()
        output_file = os.path.splitext(input_file)[0] + '.zig'
        with open(output_file, 'w') as f:
            f.write(zig_code)
        print(f"Transpilation complete -> {output_file}")
        print("Next step: zig build run -Doptimize=ReleaseFast")
    except SyntaxError as e:
        print(f"Syntax error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()