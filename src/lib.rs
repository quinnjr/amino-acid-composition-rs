//! AminoAcidComposition - Protein sequence analysis and amino acid counting
//!
//! This plugin analyzes protein sequences in FASTA format and computes
//! the amino acid composition (frequency of each amino acid).

use rayon::prelude::*;
use std::path::Path;

/// The 20 standard amino acids
pub const AMINO_ACIDS: [char; 20] = [
    'A', // Alanine
    'C', // Cysteine
    'D', // Aspartic acid
    'E', // Glutamic acid
    'F', // Phenylalanine
    'G', // Glycine
    'H', // Histidine
    'I', // Isoleucine
    'K', // Lysine
    'L', // Leucine
    'M', // Methionine
    'N', // Asparagine
    'P', // Proline
    'Q', // Glutamine
    'R', // Arginine
    'S', // Serine
    'T', // Threonine
    'V', // Valine
    'W', // Tryptophan
    'Y', // Tyrosine
];

/// Full names for amino acids
pub const AMINO_ACID_NAMES: [&str; 20] = [
    "Alanine",
    "Cysteine",
    "Aspartic acid",
    "Glutamic acid",
    "Phenylalanine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Lysine",
    "Leucine",
    "Methionine",
    "Asparagine",
    "Proline",
    "Glutamine",
    "Arginine",
    "Serine",
    "Threonine",
    "Valine",
    "Tryptophan",
    "Tyrosine",
];

/// Amino acid composition results
#[derive(Debug, Clone, Default)]
pub struct Composition {
    /// Count of each amino acid (indexed by position in AMINO_ACIDS)
    counts: [u64; 20],
    /// Total amino acid count
    total: u64,
    /// Count of unknown/non-standard characters
    unknown: u64,
}

impl Composition {
    /// Create a new empty composition
    pub fn new() -> Self {
        Self::default()
    }

    /// Get count for a specific amino acid
    #[inline]
    pub fn count(&self, aa: char) -> u64 {
        if let Some(idx) = aa_to_index(aa) {
            self.counts[idx]
        } else {
            0
        }
    }

    /// Get percentage for a specific amino acid
    #[inline]
    pub fn percentage(&self, aa: char) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        if let Some(idx) = aa_to_index(aa) {
            (self.counts[idx] as f64 / self.total as f64) * 100.0
        } else {
            0.0
        }
    }

    /// Get total amino acid count
    #[inline]
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Get count of unknown characters
    #[inline]
    pub fn unknown(&self) -> u64 {
        self.unknown
    }

    /// Get raw counts array
    pub fn counts(&self) -> &[u64; 20] {
        &self.counts
    }

    /// Merge another composition into this one
    pub fn merge(&mut self, other: &Composition) {
        for i in 0..20 {
            self.counts[i] += other.counts[i];
        }
        self.total += other.total;
        self.unknown += other.unknown;
    }
}

/// Convert amino acid character to index (0-19)
#[inline(always)]
fn aa_to_index(aa: char) -> Option<usize> {
    match aa.to_ascii_uppercase() {
        'A' => Some(0),
        'C' => Some(1),
        'D' => Some(2),
        'E' => Some(3),
        'F' => Some(4),
        'G' => Some(5),
        'H' => Some(6),
        'I' => Some(7),
        'K' => Some(8),
        'L' => Some(9),
        'M' => Some(10),
        'N' => Some(11),
        'P' => Some(12),
        'Q' => Some(13),
        'R' => Some(14),
        'S' => Some(15),
        'T' => Some(16),
        'V' => Some(17),
        'W' => Some(18),
        'Y' => Some(19),
        _ => None,
    }
}

/// Convert amino acid byte to index (0-19) - faster version for byte processing
#[inline(always)]
fn aa_byte_to_index(b: u8) -> Option<usize> {
    // Convert to uppercase
    let b = if b >= b'a' && b <= b'z' { b - 32 } else { b };

    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'D' => Some(2),
        b'E' => Some(3),
        b'F' => Some(4),
        b'G' => Some(5),
        b'H' => Some(6),
        b'I' => Some(7),
        b'K' => Some(8),
        b'L' => Some(9),
        b'M' => Some(10),
        b'N' => Some(11),
        b'P' => Some(12),
        b'Q' => Some(13),
        b'R' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'V' => Some(17),
        b'W' => Some(18),
        b'Y' => Some(19),
        _ => None,
    }
}

/// A protein sequence with header and sequence data
#[derive(Debug, Clone)]
pub struct Sequence {
    /// FASTA header (without >)
    pub header: String,
    /// Protein sequence
    pub sequence: String,
}

/// AminoAcidComposition plugin for PluMA
pub struct AminoAcidCompositionPlugin {
    /// Loaded sequences
    sequences: Vec<Sequence>,
    /// Per-sequence compositions
    per_sequence: Vec<Composition>,
    /// Overall composition
    overall: Composition,
}

impl AminoAcidCompositionPlugin {
    /// Create a new empty plugin
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            per_sequence: Vec::new(),
            overall: Composition::new(),
        }
    }

    /// Load sequences from FASTA file
    pub fn input<P: AsRef<Path>>(&mut self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        self.sequences = parse_fasta(&content);
        Ok(())
    }

    /// Load sequences from string (for testing)
    pub fn input_string(&mut self, content: &str) {
        self.sequences = parse_fasta(content);
    }

    /// Compute amino acid composition
    pub fn run(&mut self) {
        // Parallel computation for each sequence
        if self.sequences.len() >= 10 {
            self.per_sequence = self
                .sequences
                .par_iter()
                .map(|seq| compute_composition(&seq.sequence))
                .collect();
        } else {
            self.per_sequence = self
                .sequences
                .iter()
                .map(|seq| compute_composition(&seq.sequence))
                .collect();
        }

        // Merge into overall composition
        self.overall = Composition::new();
        for comp in &self.per_sequence {
            self.overall.merge(comp);
        }
    }

    /// Write results to output file
    pub fn output<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(path)?;

        writeln!(file, "Amino Acid Composition Analysis")?;
        writeln!(file, "================================")?;
        writeln!(file)?;
        writeln!(file, "Total sequences: {}", self.sequences.len())?;
        writeln!(file, "Total amino acids: {}", self.overall.total())?;
        writeln!(file, "Unknown characters: {}", self.overall.unknown())?;
        writeln!(file)?;

        writeln!(file, "Overall Composition:")?;
        writeln!(file, "--------------------")?;
        writeln!(
            file,
            "{:>3} {:>15} {:>10} {:>8}",
            "AA", "Name", "Count", "%"
        )?;

        for (i, &aa) in AMINO_ACIDS.iter().enumerate() {
            let count = self.overall.counts[i];
            let pct = self.overall.percentage(aa);
            writeln!(
                file,
                "{:>3} {:>15} {:>10} {:>7.2}%",
                aa, AMINO_ACID_NAMES[i], count, pct
            )?;
        }

        // Per-sequence breakdown if multiple sequences
        if self.sequences.len() > 1 && self.sequences.len() <= 20 {
            writeln!(file)?;
            writeln!(file, "Per-Sequence Composition:")?;
            writeln!(file, "-------------------------")?;

            for (i, (seq, comp)) in self.sequences.iter().zip(&self.per_sequence).enumerate() {
                let header = if seq.header.len() > 40 {
                    format!("{}...", &seq.header[..37])
                } else {
                    seq.header.clone()
                };
                writeln!(file, "\n{}. {} (length: {})", i + 1, header, comp.total())?;

                for (j, &aa) in AMINO_ACIDS.iter().enumerate() {
                    let count = comp.counts[j];
                    if count > 0 {
                        let pct = comp.percentage(aa);
                        write!(file, "{}:{:.1}% ", aa, pct)?;
                    }
                }
                writeln!(file)?;
            }
        }

        Ok(())
    }

    /// Get overall composition
    pub fn overall(&self) -> &Composition {
        &self.overall
    }

    /// Get per-sequence compositions
    pub fn per_sequence(&self) -> &[Composition] {
        &self.per_sequence
    }

    /// Get loaded sequences
    pub fn sequences(&self) -> &[Sequence] {
        &self.sequences
    }
}

impl Default for AminoAcidCompositionPlugin {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse FASTA format content into sequences
pub fn parse_fasta(content: &str) -> Vec<Sequence> {
    let mut sequences = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('>') {
            // Save previous sequence if any
            if !current_header.is_empty() {
                sequences.push(Sequence {
                    header: current_header,
                    sequence: current_seq,
                });
            }
            current_header = line[1..].to_string();
            current_seq = String::new();
        } else {
            current_seq.push_str(line);
        }
    }

    // Don't forget the last sequence
    if !current_header.is_empty() {
        sequences.push(Sequence {
            header: current_header,
            sequence: current_seq,
        });
    }

    sequences
}

/// Compute amino acid composition for a single sequence
pub fn compute_composition(sequence: &str) -> Composition {
    let mut comp = Composition::new();

    for b in sequence.bytes() {
        // Skip whitespace
        if b.is_ascii_whitespace() {
            continue;
        }

        if let Some(idx) = aa_byte_to_index(b) {
            comp.counts[idx] += 1;
            comp.total += 1;
        } else {
            comp.unknown += 1;
        }
    }

    comp
}

/// Compute composition using SIMD-friendly chunk processing
pub fn compute_composition_fast(sequence: &[u8]) -> Composition {
    let mut comp = Composition::new();

    // Process in chunks for better cache utilization
    const CHUNK_SIZE: usize = 1024;

    for chunk in sequence.chunks(CHUNK_SIZE) {
        for &b in chunk {
            if b.is_ascii_whitespace() {
                continue;
            }

            if let Some(idx) = aa_byte_to_index(b) {
                comp.counts[idx] += 1;
                comp.total += 1;
            } else {
                comp.unknown += 1;
            }
        }
    }

    comp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aa_to_index() {
        assert_eq!(aa_to_index('A'), Some(0));
        assert_eq!(aa_to_index('a'), Some(0));
        assert_eq!(aa_to_index('Y'), Some(19));
        assert_eq!(aa_to_index('X'), None);
        assert_eq!(aa_to_index('*'), None);
    }

    #[test]
    fn test_compute_composition() {
        let seq = "ACDEFGHIKLMNPQRSTVWY";
        let comp = compute_composition(seq);

        // Each amino acid should appear once
        assert_eq!(comp.total(), 20);
        assert_eq!(comp.unknown(), 0);

        for &aa in &AMINO_ACIDS {
            assert_eq!(comp.count(aa), 1);
            assert!((comp.percentage(aa) - 5.0).abs() < 0.01);
        }
    }

    #[test]
    fn test_composition_with_unknown() {
        let seq = "ACDEFX*BZ";
        let comp = compute_composition(seq);

        assert_eq!(comp.total(), 5); // ACDEF
        assert_eq!(comp.unknown(), 4); // X, *, B, Z
    }

    #[test]
    fn test_parse_fasta() {
        let fasta = ">seq1 description\nACDEF\nGHIKL\n>seq2\nMNPQR";
        let seqs = parse_fasta(fasta);

        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].header, "seq1 description");
        assert_eq!(seqs[0].sequence, "ACDEFGHIKL");
        assert_eq!(seqs[1].header, "seq2");
        assert_eq!(seqs[1].sequence, "MNPQR");
    }

    #[test]
    fn test_plugin_workflow() {
        let fasta = ">test\nAAAAACCCCDDDDD";

        let mut plugin = AminoAcidCompositionPlugin::new();
        plugin.input_string(fasta);
        plugin.run();

        let comp = plugin.overall();
        assert_eq!(comp.count('A'), 5);
        assert_eq!(comp.count('C'), 4);
        assert_eq!(comp.count('D'), 5);
        assert_eq!(comp.total(), 14);
    }

    #[test]
    fn test_percentage_calculation() {
        let seq = "AAAA"; // 100% Alanine
        let comp = compute_composition(seq);

        assert_eq!(comp.percentage('A'), 100.0);
        assert_eq!(comp.percentage('C'), 0.0);
    }

    #[test]
    fn test_case_insensitive() {
        let seq1 = "ACDEF";
        let seq2 = "acdef";

        let comp1 = compute_composition(seq1);
        let comp2 = compute_composition(seq2);

        assert_eq!(comp1.total(), comp2.total());
        for &aa in &AMINO_ACIDS {
            assert_eq!(comp1.count(aa), comp2.count(aa));
        }
    }
}
