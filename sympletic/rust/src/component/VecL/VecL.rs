pub mod mod_vec_n {
    #[derive(Copy, Clone)]
    pub struct VecL<T, const N: usize>([T; N]);

    #[allow(dead_code)]
    impl<T, const N: usize> VecL<T, N> {
        pub fn new() -> Self where T: Default + Copy {
            return Self([Default::default(); N]);
        }

        pub fn iter(&self) -> std::slice::Iter<'_, T> {
            return self.0.iter();
        }
    }

    impl<T, const N: usize> std::ops::Deref for VecL<T, N> {
        type Target = [T; N];
        fn deref(&self) -> &Self::Target {
            return & self.0
        }
    }

    impl<T, const N: usize> std::ops::DerefMut for VecL<T, N> {
        fn deref_mut(&mut self) -> &mut Self::Target {
            &mut self.0
        }
    }
    
    impl<T, const N: usize> std::ops::Add<VecL<T, N>> for VecL<T, N> where T: std::ops::Add<Output = T> + Copy + Default {
        type Output = VecL<T, N>;
        fn add(self, rhs: Self) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self[x] + rhs[x];
            }
            return Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::AddAssign<VecL<T, N>> for VecL<T, N> where T: std::ops::Add<Output = T> + Copy + Default {
        fn add_assign(&mut self, rhs: Self) {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self[x] + rhs[x];
            }
            *self = Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::Sub<VecL<T, N>> for VecL<T, N> where T: std::ops::Sub<Output = T> + Copy + Default {
        type Output = VecL<T, N>;
        fn sub(self, rhs: Self) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self[x] - rhs[x];
            }
            return Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::SubAssign<VecL<T, N>> for VecL<T, N> where T: std::ops::Sub<Output = T> + Copy + Default {
        fn sub_assign(&mut self, rhs: Self) {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self[x] - rhs[x];
            }
            *self = Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::Neg for VecL<T, N> where T: std::ops::Neg<Output = T> + Copy + Default {
        type Output = Self;

        fn neg(self) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = -self[x];
            }
            return Self(tmp);
        }
    }

    // impl<T, const N: usize> std::ops::Mul<VecL<T, N>> for T where T: std::ops::Mul<Output = T> {
    //     type Output = VecL<T, N>;
    //     fn mul(self, rhs: T) -> Self::Output {
    //         let mut tmp: Vec<T> = Vec::with_capacity(self.len());

    //         let mut iter = self.iter();
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
                *y = self[x] * rhs;
            }
            return Self(tmp);
        }
    }

    impl<T, A, const N: usize> std::ops::MulAssign<A> for VecL<T, N> where T: std::ops::Mul<A, Output = T> + Copy + Default, A: Copy + Default {
        fn mul_assign(&mut self, rhs: A) {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self[x] * rhs;
            }
            *self = Self(tmp);
        }
    }

    impl<T, A, const N: usize> std::ops::Div<A> for VecL<T, N> where T: std::ops::Div<A, Output = T> + Copy + Default, A: Copy + Default {
        type Output = VecL<T, N>;
        fn div(self, rhs: A) -> Self::Output {
            let mut tmp: [T; N] = [Default::default(); N];
            for (x, y) in tmp.iter_mut().enumerate() {
                *y = self[x] / rhs;
            }
            return Self(tmp);
        }
    }

    impl<T, const N: usize> std::ops::Index<usize> for VecL<T, N> {
        type Output = T;
    
        fn index(&self, _index: usize) -> &Self::Output {
            return &((*self).0[_index]);
        }
    }

    impl<T, const N: usize> std::ops::IndexMut<usize> for VecL<T, N> {
    
        fn index_mut(&mut self, _index: usize) -> &mut Self::Output {
            return &mut ((*self).0[_index]);
        }
    }
}
pub use mod_vec_n::VecL;