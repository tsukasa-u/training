// Perform a forward FFT of size 1234
use rustfft::{FftPlanner, num_complex::Complex};

let mut planner = FftPlanner::new();
let fft = planner.plan_fft_forward(1234);

let mut buffer = vec![Complex{ re: 0.0f32, im: 0.0f32 }; 1234];
fft.process(&mut buffer);
let fft_clone = Arc::clone(&fft);