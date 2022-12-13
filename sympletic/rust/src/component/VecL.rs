pub mod mod_vec_n {
    #[derive(Copy, Clone)]
    pub struct VecL<T, const N: usize>([T; N]);
    
    impl<T, const N: usize> std::ops::Add<VecL<T, N>> for VecL<T, N> where T: std::ops::Add<Output = T> + Copy + Default {
        type Output = VecL<T, N>;
        fn add(self, rhs: Self) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] + rhs.0[x];
            }
            return Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::AddAssign<VecL<T, N>> for VecL<T, N> where T: std::ops::Add<Output = T> + Copy + Default {
        fn add_assign(&mut self, rhs: Self) {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] + rhs.0[x];
            }
            *self = Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::Sub<VecL<T, N>> for VecL<T, N> where T: std::ops::Sub<Output = T> + Copy + Default {
        type Output = VecL<T, N>;
        fn sub(self, rhs: Self) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] - rhs.0[x];
            }
            return Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::SubAssign<VecL<T, N>> for VecL<T, N> where T: std::ops::Sub<Output = T> + Copy + Default {
        fn sub_assign(&mut self, rhs: Self) {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] - rhs.0[x];
            }
            *self = Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::Neg for VecL<T, N> where T: std::ops::Neg<Output = T> + Copy + Default {
        type Output = Self;

        fn neg(self) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = -self.0[x];
            }
            return Self(tmp);
        }
    }

    // impl<T, const N: usize> std::ops::Mul<VecL<T, N>> for T where T: std::ops::Mul<Output = T> {
    //     type Output = VecL<T, N>;
    //     fn mul(self, rhs: T) -> Self::Output {
    //         let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

    //         let mut iter = self.0.iter();
    //         for x in iter {
    //             tmp.push(*x * rhs);
    //         }
    //         return  tmp.into();
    //     }
    // }

    impl<T, A, const N: usize> std::ops::Mul<A> for VecL<T, N> where T: std::ops::Mul<A, Output=T> + Copy + Default, A: Copy + Default {
        type Output = VecL<T, N>;
        fn mul(self, rhs: A) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] * rhs;
            }
            return Self(tmp);
        }
    }

    impl<T, A, const N: usize> std::ops::MulAssign<A> for VecL<T, N> where T: std::ops::Mul<A, Output = T> + Copy + Default, A: Copy + Default {
        fn mul_assign(&mut self, rhs: A) {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] * rhs;
            }
            *self = Self(tmp);
        }
    }

    impl<T, A, const N: usize> std::ops::Div<A> for VecL<T, N> where T: std::ops::Div<A, Output = T> + Copy + Default, A: Copy + Default {
        type Output = VecL<T, N>;
        fn div(self, rhs: A) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self.0[x] / rhs;
            }
            return Self(tmp);
        }
    }
}
pub use mod_vec_n::VecL;