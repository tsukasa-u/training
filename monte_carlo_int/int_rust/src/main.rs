use rand::{thread_rng, Rng};

fn f(x: [f32; 10]) -> f32 {
    return x.iter().sum::<f32>().powf(2.0);
}

fn f_int( a: [(f32, f32); 10], n: i32) -> f32 {

    let mut rng = thread_rng();

    let mut sum: f32 = 0.0;
    for _t in 0..n {
            
        let mut x: [f32; 10] = [0.0; 10];
        for i in 0..10 {
            x[i] = rng.gen_range((a[i].0)..(a[i].1));
        }

        sum += f(x);
    }
    return a.iter().fold(1.0, |y, x| y*(x.1-x.0))*sum/(n as f32);
}

fn main() {
    
    let result: f32 = f_int(
        [
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0)
        ],
        10000000
    );

    print!("{}", result);
}
