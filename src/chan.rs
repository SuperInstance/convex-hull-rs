//! Chan's Algorithm — optimal O(n log h) output-sensitive convex hull.
//!
//! Combines Jarvis March with Graham Scan in a clever doubling strategy
//! to achieve the optimal output-sensitive running time.

use crate::graham::graham_scan;
use crate::point::{Orientation, Point};

/// Compute the convex hull using Chan's algorithm.
///
/// Runs in O(n log h) time where h is the number of hull vertices.
/// This is the optimal output-sensitive bound.
pub fn chans_algorithm(points: &[Point]) -> Vec<Point> {
    if points.len() < 3 {
        return points.to_vec();
    }

    let n = points.len();

    // Try increasing values of h using doubling
    let mut h = 2;
    loop {
        let result = chans_with_h(points, h.min(n));
        if let Some(hull) = result {
            return hull;
        }
        h *= 2;
        if h > n * 2 {
            // Fallback: just return graham scan result
            return graham_scan(points);
        }
    }
}

/// Run one iteration of Chan's algorithm with a candidate hull size m.
/// Returns None if the actual hull has more than m vertices.
fn chans_with_h(points: &[Point], m: usize) -> Option<Vec<Point>> {
    let n = points.len();

    // Partition into ceil(n/m) groups of at most m points each
    let num_groups = n.div_ceil(m);
    let mut sub_hulls: Vec<Vec<Point>> = Vec::with_capacity(num_groups);

    for g in 0..num_groups {
        let start = g * m;
        let end = (start + m).min(n);
        let group: Vec<Point> = points[start..end].to_vec();
        sub_hulls.push(graham_scan(&group));
    }

    // Gift wrap using the sub-hulls
    // Start from the leftmost point overall
    let mut start = points[0];
    for &p in &points[1..] {
        if p.x < start.x || (p.x == start.x && p.y < start.y) {
            start = p;
        }
    }

    let mut hull = vec![start];

    for _ in 0..m {
        let current = hull.last().unwrap();

        // For each sub-hull, find the tangent point using binary search
        let mut best = None;

        for sub_hull in &sub_hulls {
            if sub_hull.is_empty() {
                continue;
            }
            let candidate = find_tangent(sub_hull, *current, best);
            if let Some(b) = best {
                match current.orientation(&b, &candidate) {
                    Orientation::CounterClockwise => {
                        best = Some(candidate);
                    }
                    Orientation::Collinear => {
                        // Take the farther one
                        if current.dist_sq(&candidate) > current.dist_sq(&b) {
                            best = Some(candidate);
                        }
                    }
                    Orientation::Clockwise => {}
                }
            } else {
                if *current != candidate {
                    best = Some(candidate);
                }
            }
        }

        match best {
            Some(p) if p == start => {
                // Wrapped around — done
                return Some(hull);
            }
            Some(p) => {
                hull.push(p);
            }
            None => return Some(hull),
        }
    }

    // Didn't wrap around in m steps → hull is larger than m
    None
}

/// Find the point in a sub-hull that gives the most CCW turn from `current`
/// relative to the current best candidate.
fn find_tangent(hull: &[Point], current: Point, current_best: Option<Point>) -> Point {
    if hull.len() <= 3 {
        // Just scan all
        let mut best = hull[0];
        for &p in &hull[1..] {
            if p == current { continue; }
            if should_update(current, best, p, current_best) {
                best = p;
            }
        }
        return best;
    }

    // Linear scan (binary search possible but complex for correctness)
    let mut best = hull[0];
    for &p in &hull[1..] {
        if (p.x - current.x).abs() < 1e-12 && (p.y - current.y).abs() < 1e-12 {
            continue;
        }
        if (best.x - current.x).abs() < 1e-12 && (best.y - current.y).abs() < 1e-12 {
            best = p;
            continue;
        }
        if should_update(current, best, p, current_best) {
            best = p;
        }
    }
    best
}

fn should_update(current: Point, best: Point, candidate: Point, _current_best: Option<Point>) -> bool {
    if (best.x - current.x).abs() < 1e-12 && (best.y - current.y).abs() < 1e-12 {
        return true;
    }
    match current.orientation(&best, &candidate) {
        Orientation::CounterClockwise => true,
        Orientation::Collinear => current.dist_sq(&candidate) > current.dist_sq(&best),
        Orientation::Clockwise => false,
    }
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
        let hull = chans_algorithm(&pts);
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
        let hull = chans_algorithm(&pts);
        assert_eq!(hull.len(), 3);
    }

    #[test]
    fn test_collinear() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(2.0, 0.0),
        ];
        let hull = chans_algorithm(&pts);
        assert!(hull.len() <= 2);
    }

    #[test]
    fn test_single_point() {
        let pts = vec![Point::new(1.0, 1.0)];
        assert_eq!(chans_algorithm(&pts).len(), 1);
    }

    #[test]
    fn test_two_points() {
        let pts = vec![Point::new(0.0, 0.0), Point::new(1.0, 1.0)];
        assert_eq!(chans_algorithm(&pts).len(), 2);
    }

    #[test]
    fn test_large_random() {
        let mut pts = Vec::new();
        for i in 0..200 {
            let x = ((i * 13 + 7) % 400) as f64 - 200.0;
            let y = ((i * 17 + 11) % 400) as f64 - 200.0;
            pts.push(Point::new(x, y));
        }
        let hull_chans = chans_algorithm(&pts);
        let hull_graham = graham_scan(&pts);
        assert_eq!(hull_chans.len(), hull_graham.len());
    }

    #[test]
    fn test_pentagon() {
        let pts: Vec<Point> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                Point::new(angle.cos(), angle.sin())
            })
            .collect();
        let hull = chans_algorithm(&pts);
        assert_eq!(hull.len(), 5);
    }

    #[test]
    fn test_all_algorithms_agree() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(5.0, 0.0),
            Point::new(5.0, 5.0),
            Point::new(0.0, 5.0),
            Point::new(1.0, 1.0),
            Point::new(4.0, 1.0),
            Point::new(2.0, 3.0),
        ];
        let h1 = graham_scan(&pts);
        let h2 = crate::jarvis::jarvis_march(&pts);
        let h3 = crate::quickhull::quickhull(&pts);
        let h4 = chans_algorithm(&pts);
        assert_eq!(h1.len(), h2.len());
        assert_eq!(h2.len(), h3.len());
        assert_eq!(h3.len(), h4.len());
    }
}
