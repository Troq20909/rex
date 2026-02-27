// example_rmt.rex
// Random matrix theory and complex analysis
// Run: python rex_transpiler.py examples/example_rmt.rex

lab_ready;

rmt ensemble size 100;
spectrum sp of ensemble;

let pole = 1;
residue r at pole;
contour_integral ci radius 2;
