use cfd_math::high_order::dg::{numerical_flux, FluxParams, FluxType};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use leto::Array1;

fn bench_flux_allocation(c: &mut Criterion) {
    let u_l = Array1::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let u_r = Array1::from_shape_vec([5], vec![1.1, 2.1, 3.1, 4.1, 5.1]).unwrap();
    let n = Array1::from_shape_vec([3], vec![1.0, 0.0, 0.0]).unwrap();
    let params = FluxParams::new(FluxType::LaxFriedrichs);

    c.bench_function("numerical_flux_allocation", |b| {
        b.iter(|| black_box(numerical_flux(&u_l, &u_r, &n, &params)))
    });
}

criterion_group!(benches, bench_flux_allocation);
criterion_main!(benches);
