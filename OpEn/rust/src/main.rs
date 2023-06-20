use optimization_engine::{panoc::*, *};

fn f_cost(n: &[i32], u: &[f64]) -> f64 {
    u[0]
}

fn f_grad(n: &[i32], u: &[f64], grad: &mut [f64]) {
    grad[0] = 1.0
}

pub fn procedure() {
    let tolerance = 1e-12;
    let n_dim_u = 2;
    let lbfgs_memory = 10;
    let max_iters = 100;


    let n: [i32; 1] = [32];
    let mut u = [-1.5, 0.9];

    // the cache is created only ONCE
    let mut panoc_cache = PANOCCache::new(n_dim_u, tolerance, lbfgs_memory);

    let mut idx = 0;
    while idx < 100 {

        // update the function definitions (`f` and `df`)
        let df = |u: &[f64], grad: &mut [f64]| -> Result<(), SolverError> {
            f_grad(&n, &u, grad);
            Ok(())
        };
        let f = |u: &[f64], c: &mut f64| -> Result<(), SolverError> {
            *c = f_cost(&n, &u);
            Ok(())
        };

        // define the bounds at every iteration
        let bounds_min: [f64; 2] = [0.0; 2];
        let bounds_max: [f64; 2] = [10.0; 2];
        let bounds = constraints::Rectangle::new(Some(&bounds_min), Some(&bounds_max));

        // the problem definition is updated at every iteration
        let problem = Problem::new(&bounds, df, f);

        // updated instance of the solver
        let mut panoc = PANOCOptimizer::new(problem, &mut panoc_cache).with_max_iter(max_iters);

        let status = panoc.solve(&mut u).unwrap();

        idx += 1;

        // print useful information
        println!(
            "parameters: (iters = {}",
            status.iterations()
        );
        println!("u = {:#.6?}", u);
    }
}

fn main() {
    sample1::sample();
}