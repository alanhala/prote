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

    pub fn plane(points: &[Point3]) -> (Plane3, f64) {
        let n = points.len() as f64;
        let centroid = Point3::new(
            points.iter().map(|p| p.x).sum::<f64>() / n,
            points.iter().map(|p| p.y).sum::<f64>() / n,
            points.iter().map(|p| p.z).sum::<f64>() / n,
        );
        let mut covariance = nalgebra::Matrix3::<f64>::zeros();
        for point in points {
            let distance = nalgebra::Vector3::new(point.x - centroid.x, point.y - centroid.y, point.z - centroid.z);
            covariance += distance * distance.transpose();
        }
        covariance /= n;
        let eigen = nalgebra::SymmetricEigen::new(covariance);
        let (min_index, _) = eigen
            .eigenvalues
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .expect("eigen can't be empty");
        let normal_column = eigen.eigenvectors.column(min_index);
        let normal = Vector3::new(normal_column[0], normal_column[1], normal_column[2]);
        let deviation = points
            .iter()
            .map(|p| {
                let offset = Vector3::new(p.x - centroid.x, p.y - centroid.y, p.z - centroid.z);
                offset.dot(&normal).abs()
            })
            .fold(0.0, f64::max);
        (Plane3 { centroid, normal }, deviation)
    }
}

#[derive(Debug, Clone)]
pub struct Plane3 {
    centroid: Point3,
    normal: Vector3,
}

#[derive(Debug, Clone)]
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
