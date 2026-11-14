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
        let (min_index, _) = eigen.eigenvalues.iter().enumerate().min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).expect("eigen can't be empty");
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
    pub centroid: Point3,
    pub normal: Vector3,
}

impl Plane3 {
    pub fn angle(&self, other: &Self) -> f64 {
        let angle = self.normal.angle(&other.normal);
        angle.min(std::f64::consts::PI - angle)
    }

    /// The angle between the plane's normal and a direction, whichever way the normal points: 0 when
    /// the direction is perpendicular to the plane, 90° when it lies in the plane.
    pub fn normal_angle(&self, direction: &Vector3) -> f64 {
        let angle = self.normal.angle(direction);
        angle.min(std::f64::consts::PI - angle)
    }

    /// The line where the two planes meet, or `None` when they are parallel. The line's `point` is
    /// its closest point to this plane's centroid.
    pub fn intersection(&self, other: &Self) -> Option<Line3> {
        let cross = self.normal.cross(&other.normal);
        let length = cross.magnitude();
        if length < 1e-9 {
            return None;
        }
        let direction = Vector3::new(cross.x / length, cross.y / length, cross.z / length);
        // Lies in this plane and is perpendicular to the line, so it is the shortest way from the centroid to the line.
        let toward_line = direction.cross(&self.normal);
        let step = other.normal.dot(&(other.centroid - self.centroid)) / other.normal.dot(&toward_line);
        let point =
            Point3::new(self.centroid.x + toward_line.x * step, self.centroid.y + toward_line.y * step, self.centroid.z + toward_line.z * step);
        Some(Line3 { point, direction })
    }
}

#[derive(Debug, Clone)]
pub struct Line3 {
    pub point: Point3,
    pub direction: Vector3,
}

impl Line3 {
    pub fn distance(&self, point: &Point3) -> f64 {
        (*point - self.point).cross(&self.direction).magnitude()
    }
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
        Self { x: self.y * other.z - self.z * other.y, y: self.z * other.x - self.x * other.z, z: self.x * other.y - self.y * other.x }
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
        Self::Output { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planes_meet_in_a_line_closest_to_the_first_centroid() {
        // The floor (z = 0) and a wall (x = 2) meet in the line x = 2, z = 0.
        let floor = Plane3 { centroid: Point3::new(0.0, 0.0, 0.0), normal: Vector3::new(0.0, 0.0, 1.0) };
        let wall = Plane3 { centroid: Point3::new(2.0, 5.0, 7.0), normal: Vector3::new(1.0, 0.0, 0.0) };
        let line = floor.intersection(&wall).unwrap();
        assert!(line.point.distance(&Point3::new(2.0, 0.0, 0.0)) < 1e-9);
        assert!((line.distance(&floor.centroid) - 2.0).abs() < 1e-9);
        assert!((line.distance(&wall.centroid) - 7.0).abs() < 1e-9);
    }

    #[test]
    fn parallel_planes_do_not_meet() {
        let floor = Plane3 { centroid: Point3::new(0.0, 0.0, 0.0), normal: Vector3::new(0.0, 0.0, 1.0) };
        let ceiling = Plane3 { centroid: Point3::new(0.0, 0.0, 3.0), normal: Vector3::new(0.0, 0.0, -1.0) };
        assert!(floor.intersection(&ceiling).is_none());
    }
}
