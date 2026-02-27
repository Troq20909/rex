// example_lab.rex
// Lab deployment, DAQ output, embedded no_std hardening
// Run: python rex_transpiler.py examples/example_lab.rex

lab_ready;
embedded_no_std hardening for this model;

let threshold = 0.05;
validate threshold;

heat_equation_1d sensor length 0.5 dt 0.005 steps 200;

export_csv all_results to "results.csv";
lab_daq_output "sensor_data.bin" to stdout;
