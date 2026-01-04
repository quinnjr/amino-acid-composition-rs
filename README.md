# amino-acid-composition-rs

Rust implementation of the AminoAcidComposition plugin for PluMA - Protein sequence analysis and amino acid counting.

## Overview

This plugin analyzes protein sequences in FASTA format and computes the amino acid composition (frequency of each of the 20 standard amino acids).

## Installation

```bash
cargo install --path .
```

## Usage

### Command Line

```bash
amino-acid-composition input.fasta output.txt
```

### Input Format

Standard FASTA format:

```
>sp|P00533|EGFR_HUMAN Epidermal growth factor receptor
MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEV
VLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALA
>sp|P04637|P53_HUMAN Cellular tumor antigen p53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
```

### Output Format

```
Amino Acid Composition Analysis
================================

Total sequences: 2
Total amino acids: 600
Unknown characters: 0

Overall Composition:
--------------------
 AA            Name      Count        %
  A         Alanine         45    7.50%
  C        Cysteine         12    2.00%
...
```

### As a Library

```rust
use amino_acid_composition_rs::AminoAcidCompositionPlugin;

let mut plugin = AminoAcidCompositionPlugin::new();
plugin.input("proteins.fasta")?;
plugin.run();

let comp = plugin.overall();
println!("Total: {} amino acids", comp.total());
println!("Alanine: {:.2}%", comp.percentage('A'));

plugin.output("results.txt")?;
```

## Features

- Parses standard FASTA format
- Handles both upper and lowercase sequences
- Reports unknown/non-standard characters
- Parallel processing for multiple sequences
- Per-sequence and overall composition statistics

## Amino Acids

The 20 standard amino acids analyzed:

| Code | Name | Code | Name |
|------|------|------|------|
| A | Alanine | L | Leucine |
| C | Cysteine | M | Methionine |
| D | Aspartic acid | N | Asparagine |
| E | Glutamic acid | P | Proline |
| F | Phenylalanine | Q | Glutamine |
| G | Glycine | R | Arginine |
| H | Histidine | S | Serine |
| I | Isoleucine | T | Threonine |
| K | Lysine | V | Valine |
| W | Tryptophan | Y | Tyrosine |

## Building

```bash
cargo build --release
```

## Testing

```bash
cargo test
```

## Benchmarking

```bash
cargo bench
```

## References

- Original PluMA plugin: https://github.com/movingpictures83/AminoAcidComposition

## License

MIT
