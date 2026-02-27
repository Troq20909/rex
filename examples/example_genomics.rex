// example_genomics.rex
// DNA analysis, CRISPR off-target scoring, gene therapy delivery
// Run: python rex_transpiler.py examples/example_genomics.rex

lab_ready;

dna_sequence brca1 = "ATGCGTACGATCGATCGGCTAGC";
crispr_offtarget brca1 mismatches 2 pam "NGG";

dna_sequence tp53 = "GCTAGCTAGCTAGCTAGCTAGC";
crispr_offtarget tp53 mismatches 1 pam "NGG";

gene_expression tumor;
population_genetics cohort size 1000;

vector_titer aav9;
delivery_efficiency liver compartments 3;
