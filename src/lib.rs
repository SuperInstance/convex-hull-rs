//! # convex-hull-rs
//!
//! A pure-Rust library implementing multiple convex hull algorithms for 2D point sets.
//!
//! # Algorithms
//!
//! - **Graham Scan** — O(n log n) expected, classic sort-and-scan approach
//! - **Jarvis March** (Gift Wrapping) — O(nh) where h is hull size
//! - **Quickhull** — O(n log n) expected, divide-and-conquer
//! - **Chan's Algorithm** — O(n log h) optimal output-sensitive algorithm
//!
//! # Example
//!
//! ```
//! use convex_hull_rs::{graham::graham_scan, point::Point};
//!
//! let points = vec![
//!     Point::new(0.0, 0.0),
//!     Point::new(1.0, 0.0),
//!     Point::new(0.0, 1.0),
//!     Point::new(0.5, 0.5),
//! ];
//! let hull = graham_scan(&points);
//! assert_eq!(hull.len(), 3);
//! ```

pub mod point;
pub mod graham;
pub mod jarvis;
pub mod quickhull;
pub mod chan;

pub use point::Point;
