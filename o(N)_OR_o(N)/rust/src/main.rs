use rand::{thread_rng, Rng};

fn f(_t: f32, n: i32, _a: &Vec<f32>) -> (f32, u128) {

    let start:std::time::Instant = std::time::Instant::now();

    let mut _tmp:f32 = 0.0;

    let mut cnt: i32 = 0;
    let mut iter = _a.iter();
    while cnt <= n {
        _tmp = match iter.next() {
            Some(i) =>  _tmp + *i*_t.powi(cnt),
            None => std::f32::NAN
        };
        cnt += 1;
    }

    let end:std::time::Duration = start.elapsed();

    return (_tmp, end.as_millis());
}

fn g(_t: f32, n: i32, _a: &Vec<f32>) -> (f32, u128) {

    let start:std::time::Instant = std::time::Instant::now();

    let mut _tmp:f32 = 0.0;

    let mut cnt: i32 = n;
    let mut iter = _a.iter().rev();
    while cnt >= 0 {
        _tmp = match iter.next() {
            Some(i) =>  _tmp*_t + *i,
            None => std::f32::NAN
        };
        cnt -= 1;
    }

    let end:std::time::Duration = start.elapsed();

    return (_tmp, end.as_millis());
}

#[allow(non_upper_case_globals)]
fn main() {
    println!("Hello, world!");
    const N_MAX:usize = 100000000;
    let mut rng = thread_rng();
    let a: Vec<f32> = vec![rng.gen(); N_MAX + 1];

    let mut result_f: [f32; 4] = [0.0; 4];
    let mut runtime_f: [u128; 4] = [0; 4];
    let mut result_g: [f32; 4] = [0.0; 4];
    let mut runtime_g: [u128; 4] = [0; 4];
    
    const t: f32 = 0.9;

    const times: [i32; 4] = [100000, 1000000, 10000000, 100000000];

    for i in 0..4 {
        (result_f[i], runtime_f[i]) = f(t, times[i], &a);
        (result_g[i], runtime_g[i]) = g(t, times[i], &a);
    }

    for i in 0..4 { println!("{} {} {}", times[i], result_f[i], runtime_f[i]) }
    println!();
    for i in 0..4 { println!("{} {} {}", times[i], result_g[i], runtime_g[i]) }
    
}
