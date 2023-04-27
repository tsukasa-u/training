// https://helve-blog.com/posts/math/quasi-newton-method-python/

fn main() {
    let mut x: [f64; 2] = [3.0, 3.0];
    let mut B: [[f64; 2]; 2] = [
        [1.0, 0.0],
        [0.0, 1.0]
    ];

    let tau: f64 = 0.9;
    let xi: f64 = 0.3;

    let mut grad: [f64; 2] = rhs_grad(&x);

    for i in 0..10 {
        qnm(&mut x, &mut grad, &mut B, tau, xi);
        // println!("{:?}", x);
        println!("{:?}", B);
    }
}

fn rhs(x: &[f64; 2]) -> f64 {
    return 2.0*x[0].powf(2.0) + x[1].powf(2.0) + x[0]*x[1];
}

fn rhs_grad(x: &[f64; 2]) -> [f64; 2] {
    return [4.0*x[0] + x[1], x[0] + 2.0*x[1]];
}

fn line_search(x: &mut [f64; 2], grad: &[f64; 2], d: &[f64; 2], f: f64, xi: f64, tau: f64) -> f64{
    let mut alpha: f64 = 1.0;
    let tmp: f64 = xi*(grad[0]*d[0] + grad[1]*d[1]);

    while rhs(&[x[0] + alpha*d[0], x[1] + alpha*d[1]]) > f + alpha*tmp {
        alpha *= tau;
    }

    // *x = [x[0] + alpha*d[0], x[1] + alpha*d[1]];

    return alpha;
}

fn qnm(x: &mut [f64; 2], grad: &mut [f64; 2], B: &mut[[f64; 2]; 2], tau: f64, xi: f64) {
    let f: f64 = rhs(x);
    let d: [f64; 2] = {
        let tmp: f64 = B[0][0]*B[1][1] - B[1][0]*B[0][1];
        [(B[1][1]*grad[0] - B[0][1]*grad[1])/-tmp, (-B[1][0]*grad[0] + B[0][0]*grad[1])/-tmp]
    };
    // println!("{:?}", d);

    let alpha: f64 = line_search(x, grad, &d, f, xi, tau);

    // println!("{:?}", x);
    *x = [(*x)[0] + alpha*d[0], x[1] + alpha*d[1]];
    // println!("{:?}", x);

    let grad_old: [f64; 2] = grad.clone();
    *grad = rhs_grad(x);

    let s: [f64; 2] = [alpha*d[0], alpha*d[1]];
    let y: [f64; 2] = [grad[0] - grad_old[0], grad[1] - grad_old[1]];
    let Bs: [f64; 2] = [B[0][0]*s[0] + B[0][1]*s[1], B[1][0]*s[0] + B[1][1]*s[1]];
    // let Bs: [f64; 2] = [B[0][0]*s[0] + B[1][0]*s[1], B[0][1]*s[0] + B[1][1]*s[1]];

    let sBs: f64 = s[0]*Bs[0] + s[1]*Bs[1];
    let sy: f64 = s[0]*y[0] + s[1]*y[1]; 
    *B = [
        [B[0][0] - Bs[0].powf(2.0)/sBs + y[0].powf(2.0)/sy, B[0][1] - Bs[1]*Bs[0]/sBs + y[1]*y[0]/sy],
        [B[1][0] - Bs[0]*Bs[1]/sBs + y[0]*y[1]/sy, B[1][1] - Bs[1].powf(2.0)/sBs + y[1].powf(2.0)/sy]
    ];
    // *B = [
    //     [B[0][0] - Bs[0].powf(2.0)/sBs + y[0].powf(2.0)/sy, B[1][0] - Bs[1]*Bs[0]/sBs + y[1]*y[0]/sy],
    //     [B[0][1] - Bs[0]*Bs[1]/sBs + y[0]*y[1]/sy, B[1][1] - Bs[1].powf(2.0)/sBs + y[1].powf(2.0)/sy]
    // ];
    // println!("{:?}", s);
}