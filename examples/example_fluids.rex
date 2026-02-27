// example_fluids.rex
// Heat equation, thermodynamic efficiency, Burgers equation, Navier-Stokes
// Run: python rex_transpiler.py examples/example_fluids.rex

lab_ready;

heat_equation_1d rod length 1.0 dt 0.01 steps 100;
heat_equation_1d sensor length 0.5 dt 0.005 steps 200;

thermo_efficiency carnot hot 800 cold 300;
phase_transition water temp 373;

burgers_1d shock velocity 1.0 viscosity 0.01;
navier_stokes_2d_simple cavity viscosity 0.001 steps 500;
turbulence_stats cavity;
