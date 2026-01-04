//! Benchmarks for AminoAcidComposition using Criterion v0.8.0

use amino_acid_composition_rs::{
    compute_composition, compute_composition_fast, AminoAcidCompositionPlugin,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

/// Generate a random-ish protein sequence
fn generate_sequence(len: usize) -> String {
    const AAS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut seq = String::with_capacity(len);

    for i in 0..len {
        let idx = (i * 7 + 13) % AAS.len();
        seq.push(AAS[idx] as char);
    }

    seq
}

/// Generate FASTA with multiple sequences
fn generate_fasta(num_seqs: usize, seq_len: usize) -> String {
    let mut fasta = String::new();

    for i in 0..num_seqs {
        fasta.push_str(&format!(">sequence_{}\n", i + 1));

        let seq = generate_sequence(seq_len);
        // Split into 80-char lines
        for chunk in seq.as_bytes().chunks(80) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
    }

    fasta
}

fn bench_composition_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("composition_single");

    for len in [100, 1000, 10000, 100000] {
        let seq = generate_sequence(len);

        group.bench_with_input(BenchmarkId::new("standard", len), &seq, |b, s| {
            b.iter(|| {
                let comp = compute_composition(s);
                black_box(comp.total())
            })
        });

        let seq_bytes = seq.as_bytes();
        group.bench_with_input(BenchmarkId::new("fast", len), seq_bytes, |b, s| {
            b.iter(|| {
                let comp = compute_composition_fast(s);
                black_box(comp.total())
            })
        });
    }

    group.finish();
}

fn bench_full_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_pipeline");

    // Different numbers of sequences
    for (num_seqs, seq_len) in [(1, 10000), (10, 1000), (100, 1000), (1000, 500)] {
        let fasta = generate_fasta(num_seqs, seq_len);
        let label = format!("{}x{}", num_seqs, seq_len);

        group.bench_with_input(BenchmarkId::new("parse_and_run", &label), &fasta, |b, f| {
            b.iter(|| {
                let mut plugin = AminoAcidCompositionPlugin::new();
                plugin.input_string(f);
                plugin.run();
                black_box(plugin.overall().total())
            })
        });
    }

    group.finish();
}

fn bench_parallel_vs_sequential(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_vs_sequential");

    // Many small sequences vs few large sequences
    let many_small = generate_fasta(1000, 500);
    let few_large = generate_fasta(10, 50000);

    group.bench_function("many_small_1000x500", |b| {
        b.iter(|| {
            let mut plugin = AminoAcidCompositionPlugin::new();
            plugin.input_string(&many_small);
            plugin.run();
            black_box(plugin.overall().total())
        })
    });

    group.bench_function("few_large_10x50000", |b| {
        b.iter(|| {
            let mut plugin = AminoAcidCompositionPlugin::new();
            plugin.input_string(&few_large);
            plugin.run();
            black_box(plugin.overall().total())
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_composition_single,
    bench_full_pipeline,
    bench_parallel_vs_sequential
);
criterion_main!(benches);
