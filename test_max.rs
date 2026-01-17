fn main() {
    let v = vec![1.0, 2.0, 3.0];
    let m = v.iter().cloned().reduce(f64::max);
    println!("{:?}", m);
}
