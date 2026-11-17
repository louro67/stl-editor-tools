//! Indexed triangle mesh and small vector helpers.

use std::collections::HashMap;

pub type V3 = [f64; 3];

#[inline]
pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub fn mul(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
pub fn norm2(a: V3) -> f64 {
    dot(a, a)
}

/// An indexed triangle mesh. Faces are counter-clockwise seen from outside.
#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub v: Vec<V3>,
    pub f: Vec<[u32; 3]>,
}

impl Mesh {
    /// Build from a triangle soup, merging vertices with identical coordinates and dropping
    /// faces that use a vertex twice.
    pub fn from_soup(tris: impl IntoIterator<Item = [V3; 3]>) -> Mesh {
        let mut map: HashMap<[u64; 3], u32> = HashMap::new();
        let mut m = Mesh::default();
        for t in tris {
            let mut idx = [0u32; 3];
            for (k, p) in t.iter().enumerate() {
                let key = [canon(p[0]), canon(p[1]), canon(p[2])];
                let next = m.v.len() as u32;
                idx[k] = *map.entry(key).or_insert_with(|| {
                    m.v.push(*p);
                    next
                });
            }
            m.push_face(idx);
        }
        m
    }

    /// Add a face unless it is degenerate in its indices.
    pub fn push_face(&mut self, t: [u32; 3]) {
        if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
            self.f.push(t);
        }
    }

    pub fn bounds(&self) -> (V3, V3) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &self.v {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        (lo, hi)
    }

    pub fn size(&self) -> V3 {
        let (lo, hi) = self.bounds();
        sub(hi, lo)
    }

    pub fn scale(&mut self, s: f64) {
        for p in &mut self.v {
            *p = mul(*p, s);
        }
    }

    pub fn translate(&mut self, d: V3) {
        for p in &mut self.v {
            *p = add(*p, d);
        }
    }

    pub fn tri(&self, i: usize) -> [V3; 3] {
        let [a, b, c] = self.f[i];
        [self.v[a as usize], self.v[b as usize], self.v[c as usize]]
    }

    /// Signed volume (positive for an outward-oriented closed mesh).
    pub fn volume(&self) -> f64 {
        self.f
            .iter()
            .map(|t| {
                let (a, b, c) = (
                    self.v[t[0] as usize],
                    self.v[t[1] as usize],
                    self.v[t[2] as usize],
                );
                dot(a, cross(b, c))
            })
            .sum::<f64>()
            / 6.0
    }

    /// Drop vertices no face uses.
    pub fn compact(&mut self) {
        let mut remap = vec![u32::MAX; self.v.len()];
        let mut v = Vec::with_capacity(self.v.len());
        for t in &mut self.f {
            for i in t.iter_mut() {
                let r = &mut remap[*i as usize];
                if *r == u32::MAX {
                    *r = v.len() as u32;
                    v.push(self.v[*i as usize]);
                }
                *i = *r;
            }
        }
        self.v = v;
    }

    /// Per-vertex neighbor lists (sorted, unique).
    pub fn vertex_neighbors(&self) -> Vec<Vec<u32>> {
        let mut nb: Vec<Vec<u32>> = vec![Vec::new(); self.v.len()];
        for t in &self.f {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                nb[a as usize].push(b);
                nb[b as usize].push(a);
            }
        }
        for n in &mut nb {
            n.sort_unstable();
            n.dedup();
        }
        nb
    }
}

fn canon(x: f64) -> u64 {
    if x == 0.0 {
        0
    } else {
        x.to_bits()
    }
}
