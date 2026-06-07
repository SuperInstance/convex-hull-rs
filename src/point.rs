//! Point representation and geometric primitives.
//!
//! Provides the [`Point`] struct and orientation tests used by all convex hull algorithms.

use std::fmt;

/// A 2D point with `f64` coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    /// x-coordinate.
    pub x: f64,
    /// y-coordinate.
    pub y: f64,
}

impl Point {
    /// Create a new point.
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Squared Euclidean distance to another point.
    #[inline]
    pub fn dist_sq(&self, other: &Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Euclidean distance to another point.
    #[inline]
    pub fn distance(&self, other: &Point) -> f64 {
        self.dist_sq(other).sqrt()
    }

    /// Cross product of vectors (self → b) × (self → c).
    ///
    /// Positive if counter-clockwise, negative if clockwise, zero if collinear.
    #[inline]
    pub fn cross(&self, b: &Point, c: &Point) -> f64 {
        (b.x - self.x) * (c.y - self.y) - (b.y - self.y) * (c.x - self.x)
    }

    /// Orientation test for the ordered triple (self, b, c).
    #[inline]
    pub fn orientation(&self, b: &Point, c: &Point) -> Orientation {
        let val = self.cross(b, c);
        if val > 1e-10 {
            Orientation::CounterClockwise
        } else if val < -1e-10 {
            Orientation::Clockwise
        } else {
            Orientation::Collinear
        }
    }

    /// Polar angle from self to other, in [0, 2π).
    pub fn polar_angle(&self, other: &Point) -> f64 {
        let angle = (other.y - self.y).atan2(other.x - self.x);
        if angle < 0.0 { angle + 2.0 * std::f64::consts::PI } else { angle }
    }

    /// Subtract another point from this one (vector subtraction).
    #[inline]
    pub fn sub(&self, other: &Point) -> Point {
        Point::new(self.x - other.x, self.y - other.y)
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.4}, {:.4})", self.x, self.y)
    }
}

/// Orientation of three points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Counter-clockwise turn (positive cross product).
    CounterClockwise,
    /// Clockwise turn (negative cross product).
    Clockwise,
    /// Points are collinear.
    Collinear,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_new() {
        let p = Point::new(3.0, 4.0);
        assert_eq!(p.x, 3.0);
        assert_eq!(p.y, 4.0);
    }

    #[test]
    fn test_distance() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(3.0, 4.0);
        assert!((a.distance(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dist_sq() {
        let a = Point::new(1.0, 2.0);
        let b = Point::new(4.0, 6.0);
        assert!((a.dist_sq(&b) - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_orientation_ccw() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(1.0, 0.0);
        let c = Point::new(0.0, 1.0);
        assert_eq!(a.orientation(&b, &c), Orientation::CounterClockwise);
    }

    #[test]
    fn test_orientation_cw() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(0.0, 1.0);
        let c = Point::new(1.0, 0.0);
        assert_eq!(a.orientation(&b, &c), Orientation::Clockwise);
    }

    #[test]
    fn test_orientation_collinear() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(1.0, 1.0);
        let c = Point::new(2.0, 2.0);
        assert_eq!(a.orientation(&b, &c), Orientation::Collinear);
    }

    #[test]
    fn test_polar_angle() {
        let origin = Point::new(0.0, 0.0);
        let right = Point::new(1.0, 0.0);
        let up = Point::new(0.0, 1.0);
        assert!((origin.polar_angle(&right) - 0.0).abs() < 1e-10);
        assert!((origin.polar_angle(&up) - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_cross_product() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(1.0, 0.0);
        let c = Point::new(0.0, 1.0);
        assert!(a.cross(&b, &c) > 0.0);
    }
}
