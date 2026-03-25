use std::ops;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            x: x,
            y: y,
            z: z
        }
    }

    pub fn cofm1() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0
        }
    }

    pub fn cofm2() -> Self {
        Self {
            x: 1.0,
            y: 0.0,
            z: 0.0
        }
    }

    pub fn set(&mut self, x: f64, y: f64, z: f64) -> () {
        self.x = x;
        self.y = y;
        self.z = z;
    }

    // Normalises the vector
    pub fn unit(&mut self) {
        let norm_squared = self.x.powi(2) + self.y.powi(2) + self.z.powi(2);
        let norm = norm_squared.sqrt();
        self.x /= norm;
        self.y /= norm;
        self.z /= norm;
    }

    // Returns the length of the vector
    pub fn length(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    // Returns the squared length of the vector
    pub fn sqr(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    // Returns the dot product of two vectors
    pub fn dot(&self, other: &Vec3) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }

    // Returns the cross product of two vectors
    pub fn cross(&self, other: &Vec3) -> Vec3 {
        let temp_x = self.y*other.z - self.z*other.y;
        let temp_y = self.z*other.x - self.x*other.z;
        let temp_z = self.x*other.y - self.y*other.x;
        Vec3::new(temp_x, temp_y, temp_z)
    }
    
}

// in-place multiplication by f64
impl ops::MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

// in-place division by f64
impl ops::DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

// in-place addition of vector
impl ops::AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, rhs: Vec3) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}


// in-place subtraction of vector
impl ops::SubAssign<Vec3> for Vec3 {
    fn sub_assign(&mut self, rhs: Vec3) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

// Sum of two vectors
impl ops::Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

// Difference between two vectors
impl ops::Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// Multiplication of Vec3 by f64
impl ops::Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self{
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

// Multiplication of f64 by Vec3
impl ops::Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        Vec3{
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

// Division of Vec3 by f64
impl ops::Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self{
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

// Division of f64 by Vec3
impl ops::Div<Vec3> for f64 {
    type Output = Vec3;

    fn div(self, rhs: Vec3) -> Vec3 {
        Vec3{
            x: self / rhs.x,
            y: self / rhs.y,
            z: self / rhs.z,
        }
    }
}

// Negative of a vector
impl ops::Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self {
        Vec3{
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
