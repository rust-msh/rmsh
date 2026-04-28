use criterion::{criterion_group, criterion_main, Criterion};
use glam::DVec3;
use rcad_modeling::make_box_brep;
use rcad_step::writer::ExportSelection;
use rcad_step::StepWriter;

fn all_selection() -> ExportSelection<'static> {
    ExportSelection { selected_faces: &[], selected_edges: &[] }
}

fn step_roundtrip_box(c: &mut Criterion) {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let step_string = StepWriter::write_string(&brep, all_selection());
    c.bench_function("step_write_box", |bench| {
        bench.iter(|| StepWriter::write_string(&brep, all_selection()))
    });
    c.bench_function("step_parse_box", |bench| {
        bench.iter(|| rcad_step::StepReader::parse_string(&step_string).unwrap())
    });
}

criterion_group!(benches, step_roundtrip_box);
criterion_main!(benches);
