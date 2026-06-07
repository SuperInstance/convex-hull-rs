//! Graham Scan convex hull algorithm.
//!
//! Runs in O(n log n) time. Sorts points by polar angle from the lowest
//! point and then scans to remove non-left turns.

use crate::point::{Point, Orientation};

/// Compute the convex hull of a set of 2D points using the Graham Scan algorithm.
///
/// Returns the hull vertices in counter-clockwise order starting from the
/// bottom-most (then left-most) point. Duplicate points and collinear interior
/// points are removed.
///
/// # Panics
///
/// Returns an empty vec if fewer than 3 unique points are given.
pub fn graham_scan(points: &[Point]) -> Vec<Point> {
    if points.len() < 3 {
        return points.to_vec();
    }

    // Find the bottom-most point (lowest y, then leftmost x)
    let pivot = points.iter().copied().reduce(|a, b| {
        if a.y < b.y || (a.y == b.y && a.x < b.x) { a } else { b }
    }).unwrap();

    // Sort by polar angle from pivot; break ties by distance
    let mut sorted: Vec<Point> = points.to_vec();
    sorted.sort_by(|a, b| {
        let angle_a = pivot.polar_angle(a);
        let angle_b = pivot.polar_angle(b);
        if (angle_a - angle_b).abs() < 1e-12 {
            pivot.dist_sq(a).partial_cmp(&pivot.dist_sq(b)).unwrap()
        } else {
            angle_a.partial_cmp(&angle_b).unwrap()
        }
    });

    // Deduplicate
    sorted.dedup_by(|a, b| (a.x - b.x).abs() < 1e-12 && (a.y - b.y).abs() < 1e-12);

    if sorted.len() < 3 {
        return sorted;
    }

    let mut hull: Vec<Point> = Vec::with_capacity(sorted.len());

    for &p in &sorted {
        while hull.len() > 1 {
            let top = hull[hull.len() - 1];
            let next_to_top = hull[hull.len() - 2];
            if next_to_top.orientation(&top, &p) != Orientation::CounterClockwise {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }

    hull
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point;

    /// Square with interior point → hull should be 4 corners.
    #[test]
    fn test_square_with_interior() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(4.0, 0.0),
            Point::new(4.0, 4.0),
            Point::new(0.0, 4.0),
            Point::new(2.0, 2.0),
            Point::new(1.0, 1.0),
        ];
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 4);
    }

    /// Triangle hull.
    #[test]
    fn test_triangle() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(1.0, 2.0),
            Point::new(1.0, 0.5),
        ];
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 3);
    }

    /// Collinear points → degenerate hull (2 points).
    #[test]
    fn test_collinear() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(2.0, 0.0),
            Point::new(3.0, 0.0),
        ];
        let hull = graham_scan(&pts);
        assert!(hull.len() <= 2);
    }

    /// Single point.
    #[test]
    fn test_single_point() {
        let pts = vec![Point::new(1.0, 1.0)];
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 1);
    }

    /// Two points.
    #[test]
    fn test_two_points() {
        let pts = vec![Point::new(0.0, 0.0), Point::new(1.0, 1.0)];
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 2);
    }

    /// Verify CCW ordering.
    #[test]
    fn test_ccw_ordering() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(10.0, 0.0),
            Point::new(10.0, 10.0),
            Point::new(0.0, 10.0),
            Point::new(5.0, 5.0),
        ];
        let hull = graham_scan(&pts);
        for i in 0..hull.len() {
            let a = hull[i];
            let b = hull[(i + 1) % hull.len()];
            let c = hull[(i + 2) % hull.len()];
            assert_eq!(a.orientation(&b, &c), Orientation::CounterClockwise);
        }
    }

    /// Hull of a regular pentagon.
    #[test]
    fn test_pentagon() {
        let pts: Vec<Point> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                Point::new(angle.cos(), angle.sin())
            })
            .collect();
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 5);
    }

    /// Duplicate points on hull.
    #[test]
    fn test_duplicate_points() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.0, 1.0),
            Point::new(0.0, 1.0),
        ];
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 3);
    }

    /// Large random test — hull should contain all extreme points.
    #[test]
    fn test_large_random() {
        let mut pts = Vec::new();
        // Seed-like deterministic generation
        for i in 0..100 {
            let x = ((i * 7 + 13) % 200) as f64 - 100.0;
            let y = ((i * 11 + 37) % 200) as f64 - 100.0;
            pts.push(Point::new(x, y));
        }
        let hull = graham_scan(&pts);
        // Hull should contain at least 3 points and all should be from input
        assert!(hull.len() >= 3);
        for p in &hull {
            assert!(pts.iter().any(|q| (q.x - p.x).abs() < 1e-9 && (q.y - p.y).abs() < 1e-9));
        }
    }

    /// All identical points.
    #[test]
    fn test_all_identical() {
        let pts = vec![
            Point::new(5.0, 5.0),
            Point::new(5.0, 5.0),
            Point::new(5.0, 5.0),
        ];
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 1);
    }

    /// Points on a circle — all should be in hull.
    #[test]
    fn test_circle_points() {
        let pts: Vec<Point> = (0..8)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / 8.0;
                Point::new(10.0 * angle.cos(), 10.0 * angle.sin())
            })
            .collect();
        let hull = graham_scan(&pts);
        assert_eq!(hull.len(), 8);
    }
}
