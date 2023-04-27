const N: usize = 256;
// const NT: i32 = 31847*16;
// const NTINT: i32 = 31847;
const NT: i32 = 6000*128;
const NTINT: i32 = 6000;

fn runfft( ar: &mut [f64; N], ai: &mut [f64; N], flag: bool) {
    use rustfft::{Fft, FftDirection, num_complex::Complex, algorithm::Radix4};

    let fft: Radix4<f64> = Radix4::new(N, if !flag { FftDirection::Forward } else {FftDirection::Inverse} );

    // let mut buffer: Vec<Complex<f64>> = vec![Complex{ re: 10.0f64, im: 0.0f64 }; N];
    // let mut buffer: Vec<Complex<f64>> = ai.iter().map(|&s| Complex{re: s, im: 0f64}).collect();
    let mut buffer: Vec<Complex<f64>> = ar.iter().zip(ai.iter()).map(|(&s, &t)| Complex{re: s, im: t}).collect();

    // for i in 0..N {
    //     buffer[i] = Complex{ re: ((i as f64)/32.0*std::f64::consts::PI).sin(), im: 0f64};
    // }
    fft.process(&mut buffer);

    // for i in buffer {
    //     println!("{}", i.re);
    // }
    // println!("{:?}", buffer);
    if flag {
        for i in 0..N {
            ar[i] = buffer[i].re/(N as f64);
            ai[i] = buffer[i].im/(N as f64);
        }
    } else {
        for i in 0..N {
            ar[i] = buffer[i].re;
            ai[i] = buffer[i].im;
        }
    }
}

fn rhs( ar: & mut [f64; N], ai: & mut [f64; N], akr: &mut [f64; N], aki: &mut [f64; N], delta: f64, dt: f64) {
/* INPUT  : ar , ai  in wave number space  */
/* OUTPUT : akr aki  = rhs *dt in wave number space  */
     
    let mut ur: [f64; N] = [0f64; N];
    let mut ui: [f64; N] = [0f64; N];
    let mut dur: [f64; N] = [0f64; N];
    let mut dui: [f64; N] = [0f64; N];
    let mut di: f64;
    let mut d3i: f64;

    let dx: f64 = 2.0/(N as f64);

    for i in 0..N {
        ur[i]=ar[i];
        ui[i]=ai[i];
    }
    //  fft1(ur,ui,N,ITEL,1);
    runfft(&mut ur, &mut ui, true);

    for i in 0..(N/2+1)  {
        di=i as f64;
	    dur[i]=-1.0*std::f64::consts::PI*di*ai[i];
        dui[i]= 1.0*std::f64::consts::PI*di*ar[i];
    }
    for i in 1..(N/2) {
        dur[N-i]=dur[i];
        dui[N-i]=-dui[i];
    }
    // fft1(dur,dui,N,ITEL,1);
    runfft(&mut dur, &mut dui, true);

    for i in 0..N {
        dur[i]=ur[i]*dur[i];
        dui[i]=ui[i]*dui[i];
    }

    // fft1(dur,dui,N,ITEL,0);
    runfft(&mut dur, &mut dui, false);
    for i in 0..(N/2+1) {
        di=i as f64;
        d3i=di*di*di*std::f64::consts::PI*std::f64::consts::PI*std::f64::consts::PI;
        akr[i]=-dt*(dur[i]+1.0*d3i*delta*ai[i]);
        aki[i]=-dt*(dui[i]-1.0*d3i*delta*ar[i]);
    }
    for i in 1..(N/2) {
        di=i as f64;
        akr[N-i]=akr[i];
        aki[N-i]=-aki[i];
    }

     return;
}

fn rk( ar: & mut [f64; N], ai: & mut [f64; N], delta: f64, dx: f64, dt: f64) {
    let mut aro: [f64; N] = [0f64; N];
    let mut aio: [f64; N] = [0f64; N];
    let mut ak1r: [f64; N] = [0f64; N];
    let mut ak1i: [f64; N] = [0f64; N];
    let mut ak2r: [f64; N] = [0f64; N];
    let mut ak2i: [f64; N] = [0f64; N];
    let mut ak3r: [f64; N] = [0f64; N];
    let mut ak3i: [f64; N] = [0f64; N];
    let mut ak4r: [f64; N] = [0f64; N];
    let mut ak4i: [f64; N] = [0f64; N];
    for i in 0..N {
        aro[i]=ar[i];
        aio[i]=ai[i];
    }
    rhs(ar, ai, &mut ak1r, &mut ak1i, delta, dt);
     for i in 0..N {
        ar[i]=aro[i]+ak1r[i]*0.5;
        ai[i]=aio[i]+ak1i[i]*0.5;
    }
    rhs(ar,ai,&mut ak2r,&mut ak2i,delta,dt);
    for i in 0..N {
        ar[i]=aro[i]+ak2r[i]*0.5;
        ai[i]=aio[i]+ak2i[i]*0.5;
    }
    rhs(ar,ai,&mut ak3r,&mut ak3i,delta,dt);
    for i in 0..N {
        ar[i]=aro[i]+ak3r[i];
        ai[i]=aio[i]+ak3i[i];
    }
    rhs(ar,ai,&mut ak4r,&mut ak4i,delta,dt);
    for i in 0..N {
        ar[i]=aro[i]+(ak1r[i]+2.0*ak2r[i]+2.0*ak3r[i]+ak4r[i])/6.0;
        ai[i]=aio[i]+(ak1i[i]+2.0*ak2i[i]+2.0*ak3i[i]+ak4i[i])/6.0;
        // println!("{} {}",ar[i], aro[i])
    }
    return;
}

fn initc(ar: &mut [f64; N], ai: &mut [f64; N], dx: f64) {
    let mut u: [f64; N] = [0f64; N];
    for i in 0..N {
        u[i]=(std::f64::consts::PI*dx*(i as f64)).cos();
    }
    for i in 0..N {
        ar[i]=u[i];
        ai[i]=0.0;
    }
    //   fft1(ar,ai,N,ITEL,0);
    runfft(ar, ai, false);
    return;
}

fn main() {
    let u: [f64; N];
    let mut ar: [f64; N] = [0f64; N];
    let mut ai: [f64; N] = [0f64; N];
    let mut vr: [f64; N] = [0f64; N];
    let mut vi: [f64; N] = [0f64; N];
    let dx: f64;
    let dt: f64;
    let delta: f64;

    dx=2.0/(N as f64);
    dt=0.00001;
    delta=0.022*0.022;

    println!("t,x,u"); 

    initc(&mut ar, &mut ai, dx);

    for j in 0..N {
        vr[j]=ar[j];
        vi[j]=ai[j];
    }
    // fft1(vr,vi,N,ITEL,1);
    runfft(&mut vr,&mut vi, true);

    for j in 0..N {
        println!("{}, {}, {}", 0, j, vr[j]); 
    }  
    println!("");

    for i in 1..NT {
        rk(&mut ar, &mut ai,delta,dx,dt);
        if i%NTINT==0 {

            for j in 0..N {
                vr[j]=ar[j];vi[j]=ai[j];
            }
            // fft1(vr,vi,N,ITEL,1);
            runfft(&mut vr,&mut vi, true);
            
            for j in 0..N {
                println!("{}, {}, {}", i, j, vr[j]); 
            }  
            println!("");

        }
    }

    return;
}
