# Contributing to REX

Thank you for your interest in contributing. Rex is a free, open research tool and all contributions are welcome — from bug fixes and new domain keywords to example models and documentation improvements.

---

## Ways to Contribute

- **Bug reports** – Found a parse error or incorrect Zig output? Open an issue with the `.rex` file that triggered it.
- **New domain keywords** – Have a scientific domain Rex doesn't cover? Propose it as an issue or submit a PR.
- **Example models** – Real `.rex` files from your research are enormously valuable to other researchers.
- **Documentation** – Clearer explanations, better examples, translated READMEs.
- **Zig backend improvements** – Better comptime patterns, more efficient output, additional Zig stdlib usage.

---

## Getting Started

1. Fork the repository and clone your fork:

```bash
git clone https://github.com/yourusername/rex.git
cd rex
```

2. Make sure you have Python 3.8+ and Zig 0.15.2 installed.

3. Run the full example to verify your setup:

```bash
python rex_transpiler.py examples/example_full.rex
zig build run -Doptimize=ReleaseFast
```

---

## Adding a New Domain Keyword

Adding a new keyword to Rex requires four steps, all in `rex_transpiler.py`:

**Step 1 – Add to the keyword set in `Lexer.tokenize()`:**
```python
keywords = {
    ...
    'your_new_keyword',
    ...
}
```

**Step 2 – Write a parse method in `Parser`:**
```python
def parse_your_new_keyword_stmt(self) -> str:
    self.expect('YOUR_NEW_KEYWORD')
    name = self.expect('ID')
    # parse any parameters your keyword needs
    self.expect(';')
    return (f'    const result = your_new_keyword_comptime({name});\n'
            f'    std.debug.print("Result: {{}}\n", .{{result}});')
```

**Step 3 – Wire it into `parse_program()`:**
```python
elif t == 'YOUR_NEW_KEYWORD':
    main_body.append(self.parse_your_new_keyword_stmt())
```

**Step 4 – Add an example `.rex` file in `examples/` showing realistic usage.**

---

## Code Style

- Python code follows standard PEP 8 conventions.
- Parse methods are named `parse_<keyword>_stmt()` for statement-level keywords.
- All Zig output uses `const` (not `let`), `std.debug.print`, and the Arena allocator pattern.
- Each parse method should handle exactly one keyword and end with `self.expect(';')`.
- Keep parse methods short and single-purpose — one keyword, one method.

---

## Submitting a Pull Request

1. Create a branch with a descriptive name:
```bash
git checkout -b add-fourier-transform-keyword
```

2. Make your changes and test them:
```bash
python rex_transpiler.py examples/your_new_example.rex
zig build run -Doptimize=ReleaseFast
```

3. Commit with a clear message:
```bash
git commit -m "Add fourier_transform keyword with comptime Zig emission"
```

4. Push and open a PR against `main`. Describe what the keyword does and include a sample `.rex` file in the PR description.

---

## Reporting Bugs

Please include in your issue:

- The `.rex` source file or a minimal reproduction
- The error message or incorrect output you received
- Your Python version (`python --version`) and Zig version (`zig version`)

---

## License

By contributing to Rex you agree that your contributions will be licensed under the MIT License.
