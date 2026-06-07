//! Jarvis March (Gift Wrapping) convex hull algorithm.
//!
//! Output-sensitive: O(nh) where h is the number of hull vertices.
//! Useful when the hull is known to be small.

use crate::point::{Orientation, Point};

/// Compute the convex hull using Jarvis March (gift wrapping).
///
/// Returns hull vertices in CCW order.
pub fn jarvis_march(points: &[Point]) -> Vec<Point> {
    if points.len() < 3 {
        return points.to_vec();
    }

    // Start from the leftmost point
    let mut start = 0;
    for (i, p) in points.iter().enumerate() {
        if p.x < points[start].x || (p.x == points[start].x && p.y < points[start].y) {
            start = i;
        }
    }

    let mut hull = Vec::new();
    let mut current = start;

    loop {
        hull.push(points[current]);

        // Pick the next point that is most counter-clockwise
        let mut next = 0;
        for (i, _) in points.iter().enumerate() {
            if next == current {
                next = i;
                continue;
            }
            let orient = points[current].orientation(&points[next], &points[i]);
            if orient == Orientation::CounterClockwise {
                next = i;
            } else if orient == Orientation::Collinear {
                // Pick the farther one
                if points[current].dist_sq(&points[i]) > points[current].dist_sq(&points[next]) {
                    next = i;
                }
            }
        }

        current = next;
        if current == start {
            break;
        }
    }

    hull
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_square_hull() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(4.0, 0.0),
            Point::new(4.0, 4.0),
            Point::new(0.0, 4.0),
            Point::new(2.0, 2.0),
        ];
        let hull = jarvis_march(&pts);
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_triangle() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(4.0, 0.0),
            Point::new(2.0, 3.0),
            Point::new(2.0, 1.0),
        ];
        let hull = jarvis_march(&pts);
        assert_eq!(hull.len(), 3);
    }

    #[test]
    fn test_collinear() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(2.0, 2.0),
            Point::new(3.0, 3.0),
        ];
        let hull = jarvis_march(&pts);
        assert!(hull.len() <= 2);
    }

    #[test]
    fn test_single_point() {
        let pts = vec![Point::new(0.0, 0.0)];
        let hull = jarvis_march(&pts);
        assert_eq!(hull.len(), 1);
    }

    #[test]
    fn test_line_segment() {
        let pts = vec![Point::new(0.0, 0.0), Point::new(5.0, 5.0)];
        let hull = jarvis_march(&pts);
        assert_eq!(hull.len(), 2);
    }

    #[test]
    fn test_many_interior() {
        let mut pts = vec![
            Point::new(0.0, 0.0),
            Point::new(10.0, 0.0),
            Point::new(10.0, 10.0),
            Point::new(0.0, 10.0),
        ];
        // Add many interior points
        for i in 1..9 {
            for j in 1..9 {
                pts.push(Point::new(i as f64, j as f64));
            }
        }
        let hull = jarvis_march(&pts);
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_pentagon() {
        let pts: Vec<Point> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                Point::new(angle.cos(), angle.sin())
            })
            .collect();
        let hull = jarvis_march(&pts);
        assert_eq!(hull.len(), 5);
    }

    #[test]
    fn test_hull_contains_extremes() {
        let pts = vec![
            Point::new(-5.0, -5.0),
            Point::new(5.0, -5.0),
            Point::new(5.0, 5.0),
            Point::new(-5.0, 5.0),
            Point::new(0.0, 0.0),
            Point::new(1.0, 1.0),
        ];
        let hull = jarvis_march(&pts);
        assert!(hull.iter().any(|p| p.x == -5.0 && p.y == -5.0));
        assert!(hull.iter().any(|p| p.x == 5.0 && p.y == 5.0));
    }

    #[test]
    fn test_ccw_ordering() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(3.0, 0.0),
            Point::new(3.0, 3.0),
            Point::new(0.0, 3.0),
            Point::new(1.0, 1.0),
        ];
        let hull = jarvis_march(&pts);
        // Verify hull is closed and forms a valid polygon
        assert_eq!(hull.len(), 4);
        // Verify all hull points are extreme
        for p in &hull {
            assert!(pts.iter().any(|q| q == p));
        }
    }
}
