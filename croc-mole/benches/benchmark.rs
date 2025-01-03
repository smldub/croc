use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use croc_mole::mol::{MolBuilder};

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut molTi = MolBuilder::new();
    let _ = molTi.set_atoms_str("Ti 0 0 0");
    let _ = molTi.set_basis("cc-pvdz");
    let molTi = molTi.build();
    let mut molTi2 = MolBuilder::new();
    let _ = molTi2.set_atoms_str("Ti 0 0 0; Ti 1 0 0");
    let _ = molTi2.set_basis("cc-pvdz");
    let molTi2 = molTi2.build();
    let mut molH2 = MolBuilder::new();
    let _ = molH2.set_atoms_str("H 0 0 0; H 1 0 0");
    let _ = molH2.set_basis("sto-3g");
    let molH2 = molH2.build();
    let mut group = c.benchmark_group("int1e_ovlp");
    for mol  in [molH2,molTi, molTi2].iter() {
        group.sample_size(2000);
        group.bench_function(BenchmarkId::from_parameter(mol.get_nao()), |b| b.iter(|| mol.get_int1e_ovlp()));
        //c.bench_function("int2e", |b| b.iter(|| mol.get_int2e()));
    }
}
criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
