# convex-hull-rs

Convex hull algorithms in pure Rust: Graham scan, Jarvis march, Quickhull, Chan's algorithm.

## Features

- **Graham Scan** — O(n log n) classic sort-and-scan
- **Jarvis March** — O(nh) output-sensitive gift wrapping
- **Quickhull** — O(n log n) expected divide-and-conquer
- **Chan's Algorithm** — O(n log h) optimal output-sensitive

No external dependencies. Pure `std` Rust.

## Usage

```rust
use convex_hull_rs::{graham::graham_scan, point::Point};

let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(0.0, 1.0),
    Point::new(0.5, 0.5),
];
let hull = graham_scan(&points);
assert_eq!(hull.len(), 3);
```

## License

MIT OR Apache-2.0
