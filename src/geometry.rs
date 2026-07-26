use std::ops::Sub;

#[derive(Debug, Clone, Copy)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn distance(&self, other: &Self) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)).sqrt()
    }
}

#[derive(Debug)]
pub struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn magnitude(&self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn angle(&self, other: &Self) -> f64 {
        self.cross(other).magnitude().atan2(self.dot(other))
    }
}

impl Sub for Point3 {
    type Output = Vector3;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: rhs.x - self.x,
            y: rhs.y - self.y,
            z: rhs.z - self.z,
        }
    }
}
