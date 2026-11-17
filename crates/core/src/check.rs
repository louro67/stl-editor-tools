//! Mesh diagnostics: what's wrong with an input, and proof that an output is a clean solid.

use serde::Serialize;

use crate::contour::FxMap;
use crate::mesh::{norm, sub, Mesh};

fn find(p: &mut [u32], mut a: u32) -> u32 {
    while p[a as usize] != a {
        p[a as usize] = p[p[a as usize] as usize];
        a = p[a as usize];
    }
    a
}

fn union(p: &mut [u32], a: u32, b: u32) {
    let (ra, rb) = (find(p, a), find(p, b));
    if ra != rb {
        p[ra.max(rb) as usize] = ra.min(rb);
    }
}

struct Edges {
    boundary: Vec<(u32, u32)>,
    nonmanifold: usize,
    inconsistent: usize,
}

fn edges(m: &Mesh) -> Edges {
    let mut half: FxMap<u64, u32> = FxMap::default();
    half.reserve(m.f.len() * 3);
    for t in &m.f {
        for k in 0..3 {
            let (a, b) = (t[k] as u64, t[(k + 1) % 3] as u64);
            *half.entry(a << 32 | b).or_insert(0) += 1;
        }
    }
    let mut out = Edges {
        boundary: Vec::new(),
        nonmanifold: 0,
        inconsistent: 0,
    };
    for (&key, &n) in &half {
        let (a, b) = ((key >> 32) as u32, key as u32);
        let rev = half
            .get(&((b as u64) << 32 | a as u64))
            .copied()
            .unwrap_or(0);
        if a < b || rev == 0 {
            let total = n + rev;
            if total == 1 {
                out.boundary.push((a, b));
            } else if total > 2 {
                out.nonmanifold += 1;
            } else if n != 1 || rev != 1 {
                out.inconsistent += 1;
            }
        }
    }
    out
}

fn components(m: &Mesh) -> usize {
    let mut p: Vec<u32> = (0..m.v.len() as u32).collect();
    for t in &m.f {
        union(&mut p, t[0], t[1]);
        union(&mut p, t[1], t[2]);
    }
    let mut used = vec![false; m.v.len()];
    for t in &m.f {
        for &i in t {
            used[i as usize] = true;
        }
    }
    (0..m.v.len())
        .filter(|&i| used[i] && find(&mut p, i as u32) == i as u32)
        .count()
}

/// What an input mesh looks like.
#[derive(Clone, Debug, Serialize)]
pub struct Inspection {
    pub triangles: usize,
    pub vertices: usize,
    pub size_mm: [f64; 3],
    pub watertight: bool,
    pub winding_consistent: bool,
    pub components: usize,
    pub boundary_edges: usize,
    pub nonmanifold_edges: usize,
    pub holes: usize,
    pub largest_hole_perimeters_mm: Vec<f64>,
}

pub fn inspect(m: &Mesh) -> Inspection {
    let e = edges(m);
    // Boundary loops: connected groups of boundary edges.
    let mut p: Vec<u32> = (0..m.v.len() as u32).collect();
    for &(a, b) in &e.boundary {
        union(&mut p, a, b);
    }
    let mut per: FxMap<u32, f64> = FxMap::default();
    for &(a, b) in &e.boundary {
        let r = find(&mut p, a);
        *per.entry(r).or_insert(0.0) += norm(sub(m.v[a as usize], m.v[b as usize]));
    }
    let mut perims: Vec<f64> = per.values().map(|x| (x * 100.0).round() / 100.0).collect();
    perims.sort_by(|a, b| b.total_cmp(a));
    let holes = perims.len();
    perims.truncate(10);
    let s = m.size();
    Inspection {
        triangles: m.f.len(),
        vertices: m.v.len(),
        size_mm: s.map(|x| (x * 1000.0).round() / 1000.0),
        watertight: e.boundary.is_empty() && e.nonmanifold == 0,
        winding_consistent: e.inconsistent == 0,
        components: components(m),
        boundary_edges: e.boundary.len(),
        nonmanifold_edges: e.nonmanifold,
        holes,
        largest_hole_perimeters_mm: perims,
    }
}

/// Checks on an output.
#[derive(Clone, Debug, Serialize)]
pub struct Verification {
    /// Every edge has exactly two faces, oriented consistently.
    pub watertight: bool,
    /// Watertight with positive volume: a printable solid.
    pub is_volume: bool,
    pub bodies: usize,
    pub volume_mm3: f64,
    pub euler: i64,
}

pub fn verify(m: &Mesh) -> Verification {
    let e = edges(m);
    let watertight = e.boundary.is_empty() && e.nonmanifold == 0 && e.inconsistent == 0;
    let volume = m.volume();
    let nedges = (m.f.len() * 3 / 2) as i64;
    Verification {
        watertight,
        is_volume: watertight && volume > 0.0,
        bodies: components(m),
        volume_mm3: (volume * 1000.0).round() / 1000.0,
        euler: m.v.len() as i64 - nedges + m.f.len() as i64,
    }
}
