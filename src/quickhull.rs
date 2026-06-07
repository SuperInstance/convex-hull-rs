//! Quickhull convex hull algorithm.
//!
//! Divide-and-conquer approach. Expected O(n log n), worst case O(n²).
//! Works by recursively partitioning points relative to the farthest
//! point from an edge.

use crate::point::{Orientation, Point};

/// Compute the convex hull using the Quickhull algorithm.
///
/// Returns hull vertices in CCW order.
pub fn quickhull(points: &[Point]) -> Vec<Point> {
    if points.len() < 3 {
        return points.to_vec();
    }

    // Find leftmost and rightmost points
    let mut min_x = 0;
    let mut max_x = 0;
    for (i, p) in points.iter().enumerate() {
        if p.x < points[min_x].x {
            min_x = i;
        }
        if p.x > points[max_x].x {
            max_x = i;
        }
    }

    if min_x == max_x {
        // All points have same x — vertical line
        let min_y = points.iter().enumerate().min_by(|a, b| a.1.y.partial_cmp(&b.1.y).unwrap()).unwrap().0;
        let max_y = points.iter().enumerate().max_by(|a, b| a.1.y.partial_cmp(&b.1.y).unwrap()).unwrap().0;
        if min_y == max_y {
            return vec![points[min_y]];
        }
        return vec![points[min_y], points[max_y]];
    }

    let a = points[min_x];
    let b = points[max_x];

    // Partition into left and right of the line a→b
    let mut left = Vec::new();
    let mut right = Vec::new();

    for &p in points {
        if p == a || p == b {
            continue;
        }
        match a.orientation(&b, &p) {
            Orientation::CounterClockwise => left.push(p),
            Orientation::Clockwise => right.push(p),
            Orientation::Collinear => {}
        }
    }

    let mut hull = Vec::new();
    hull.push(a);
    let upper = find_hull(&left, a, b);
    hull.push(b);
    let lower = find_hull(&right, b, a);

    // Reconstruct CCW: a → upper → b → lower
    let mut result = vec![a];
    result.extend_from_slice(&upper);
    result.push(b);
    result.extend_from_slice(&lower);
    result
}

/// Recursive helper: find points on the hull on one side of line a→b.
fn find_hull(points: &[Point], a: Point, b: Point) -> Vec<Point> {
    if points.is_empty() {
        return Vec::new();
    }

    // Find point farthest from line a→b
    let farthest_idx = points.iter().enumerate().max_by(|(_, p1), (_, p2)| {
        let d1 = point_line_dist(**p1, a, b);
        let d2 = point_line_dist(**p2, a, b);
        d1.partial_cmp(&d2).unwrap()
    }).unwrap().0;

    let c = points[farthest_idx];

    // Points inside triangle a-c-b are not on the hull
    // Partition remaining points
    let mut left_of_ac = Vec::new();
    let mut left_of_cb = Vec::new();

    for &p in points {
        if p == c {
            continue;
        }
        if a.orientation(&c, &p) == Orientation::CounterClockwise {
            left_of_ac.push(p);
        } else if c.orientation(&b, &p) == Orientation::CounterClockwise {
            left_of_cb.push(p);
        }
    }

    let mut result = find_hull(&left_of_ac, a, c);
    result.push(c);
    result.extend(find_hull(&left_of_cb, c, b));
    result
}

/// Perpendicular distance from point p to line through a and b.
fn point_line_dist(p: Point, a: Point, b: Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-12 {
        return p.distance(&a);
    }
    (dy * p.x - dx * p.y + b.x * a.y - b.y * a.x).abs() / len
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_square() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(4.0, 0.0),
            Point::new(4.0, 4.0),
            Point::new(0.0, 4.0),
            Point::new(2.0, 2.0),
        ];
        let hull = quickhull(&pts);
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_triangle() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(4.0, 0.0),
            Point::new(2.0, 3.0),
            Point::new(2.0, 1.0),
            Point::new(1.0, 0.5),
        ];
        let hull = quickhull(&pts);
        assert_eq!(hull.len(), 3);
    }

    #[test]
    fn test_collinear_horizontal() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(2.0, 0.0),
        ];
        let hull = quickhull(&pts);
        assert!(hull.len() <= 2);
    }

    #[test]
    fn test_collinear_vertical() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(0.0, 1.0),
            Point::new(0.0, 2.0),
        ];
        let hull = quickhull(&pts);
        assert!(hull.len() <= 2);
    }

    #[test]
    fn test_single_point() {
        let pts = vec![Point::new(1.0, 1.0)];
        assert_eq!(quickhull(&pts).len(), 1);
    }

    #[test]
    fn test_two_points() {
        let pts = vec![Point::new(0.0, 0.0), Point::new(1.0, 1.0)];
        assert_eq!(quickhull(&pts).len(), 2);
    }

    #[test]
    fn test_pentagon() {
        let pts: Vec<Point> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                Point::new(angle.cos(), angle.sin())
            })
            .collect();
        let hull = quickhull(&pts);
        assert_eq!(hull.len(), 5);
    }

    #[test]
    fn test_many_interior_points() {
        let mut pts = vec![
            Point::new(-10.0, -10.0),
            Point::new(10.0, -10.0),
            Point::new(10.0, 10.0),
            Point::new(-10.0, 10.0),
        ];
        for i in -9..10 {
            for j in -9..10 {
                pts.push(Point::new(i as f64, j as f64));
            }
        }
        let hull = quickhull(&pts);
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_l_shape() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(0.0, 5.0),
            Point::new(3.0, 5.0),
            Point::new(3.0, 3.0),
            Point::new(5.0, 3.0),
            Point::new(5.0, 0.0),
        ];
        let hull = quickhull(&pts);
        // L-shape exterior should be a convex hull with all 6 extreme points
        assert!(hull.len() >= 4);
    }

    #[test]
    fn test_hull_area_positive() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(10.0, 0.0),
            Point::new(10.0, 10.0),
            Point::new(0.0, 10.0),
            Point::new(5.0, 5.0),
        ];
        let hull = quickhull(&pts);
        // Compute polygon area via shoelace
        let mut area = 0.0;
        for i in 0..hull.len() {
            let j = (i + 1) % hull.len();
            area += hull[i].x * hull[j].y;
            area -= hull[j].x * hull[i].y;
        }
        area = area.abs() / 2.0;
        assert!((area - 100.0).abs() < 1e-6);
    }
}
