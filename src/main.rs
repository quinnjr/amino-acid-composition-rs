//! AminoAcidComposition CLI - Protein sequence analysis tool
//!
//! Usage: amino-acid-composition <input.fasta> <output.txt>

use amino_acid_composition_rs::AminoAcidCompositionPlugin;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input.fasta> <output.txt>", args[0]);
        eprintln!();
        eprintln!("Analyzes amino acid composition of protein sequences.");
        eprintln!("Input: FASTA format protein sequences");
        eprintln!("Output: Amino acid counts and percentages");
        process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    let mut plugin = AminoAcidCompositionPlugin::new();

    // Input phase
    if let Err(e) = plugin.input(input_file) {
        eprintln!("Error reading input file '{}': {}", input_file, e);
        process::exit(1);
    }

    eprintln!("Loaded {} sequences", plugin.sequences().len());

    // Run phase
    plugin.run();

    let overall = plugin.overall();
    eprintln!(
        "Analyzed {} amino acids ({} unknown characters)",
        overall.total(),
        overall.unknown()
    );

    // Output phase
    if let Err(e) = plugin.output(output_file) {
        eprintln!("Error writing output file '{}': {}", output_file, e);
        process::exit(1);
    }

    eprintln!("Results written to '{}'", output_file);
}
