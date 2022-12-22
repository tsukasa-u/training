pub mod mod_vec3 {
    #[derive(Copy, Clone)]
    pub struct Vec3(pub f32, pub f32, pub f32) ;

    impl std::ops::Add<Vec3> for Vec3 {
        type Output = Vec3;

        fn add(self, rhs: Vec3) -> Vec3 {
            return Vec3(
                self.0 + rhs.0,
                self.1 + rhs.1,
                self.2 + rhs.2
            );
        }
    }

    impl std::ops::AddAssign<Vec3> for Vec3 {
        fn add_assign(&mut self, rhs: Vec3) {
            *self =  Self(
                self.0 + rhs.0,
                self.1 + rhs.1,
                self.2 + rhs.2
            );
        }
    }

    impl std::ops::Sub<Vec3> for Vec3 {
        type Output = Vec3;

        fn sub(self, rhs: Vec3) -> Vec3 {
            return Vec3(
                self.0 - rhs.0,
                self.1 - rhs.1,
                self.2 - rhs.2
            );
        }
    }

    impl std::ops::SubAssign<Vec3> for Vec3 {
        fn sub_assign(&mut self, rhs: Vec3) {
            *self =  Self(
                self.0 - rhs.0,
                self.1 - rhs.1,
                self.2 - rhs.2
            );
        }
    }

    impl std::ops::Mul<Vec3> for f32 {
        type Output = Vec3;

        fn mul(self, rhs: Vec3) -> Vec3 {
            return Vec3(
                self * rhs.0,
                self * rhs.1,
                self * rhs.2
            );
        }
    }

    impl std::ops::Mul<f32> for Vec3 {
        type Output = Vec3;

        fn mul(self, rhs: f32) -> Vec3 {
            return Vec3(
                self.0 * rhs,
                self.1 * rhs,
                self.2 * rhs
            );
        }
    }
}
pub use mod_vec3::Vec3;