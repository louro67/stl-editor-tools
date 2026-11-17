//! Bounding volume hierarchy over triangles: closest-point distance and a fast
//! generalized winding number (Barnes-Hut with a dipole far field, after Barill et al.,
//! "Fast Winding Numbers for Soups and Clouds", SIGGRAPH 2018).

use crate::mesh::{add, cross, dot, mul, norm, norm2, sub, Mesh, V3};

const LEAF: usize = 4;
const INV_4PI: f64 = 1.0 / (4.0 * std::f64::consts::PI);

#[derive(Clone, Copy, Default)]
struct Node {
    lo: V3,
    hi: V3,
    /// Leaf: first triangle. Interior: index of the left child (right is left + 1).
    first: u32,
    /// Triangles in a leaf; 0 for an interior node.
    count: u32,
    /// Sum of area vectors (half cross products) below this node.
    n: V3,
    /// Area-weighted centroid.
    c: V3,
    /// Radius of a ball around `c` holding every triangle below.
    r: f64,
    area: f64,
}

pub struct Bvh {
    nodes: Vec<Node>,
    tris: Vec<[V3; 3]>,
}

impl Bvh {
    pub fn new(m: &Mesh) -> Bvh {
        let tris: Vec<[V3; 3]> = (0..m.f.len()).map(|i| m.tri(i)).collect();
        let cent: Vec<V3> = tris
            .iter()
            .map(|t| mul(add(add(t[0], t[1]), t[2]), 1.0 / 3.0))
            .collect();
        let mut idx: Vec<u32> = (0..tris.len() as u32).collect();
        let mut nodes = vec![Node::default()];
        let mut stack = vec![(0usize, 0usize, tris.len())];
        while let Some((ni, s, e)) = stack.pop() {
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            let (mut clo, mut chi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for &t in &idx[s..e] {
                for p in tris[t as usize] {
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k]);
                        hi[k] = hi[k].max(p[k]);
                    }
                }
                let c = cent[t as usize];
                for k in 0..3 {
                    clo[k] = clo[k].min(c[k]);
                    chi[k] = chi[k].max(c[k]);
                }
            }
            nodes[ni].lo = lo;
            nodes[ni].hi = hi;
            if e - s <= LEAF {
                nodes[ni].first = s as u32;
                nodes[ni].count = (e - s) as u32;
                continue;
            }
            let ext = sub(chi, clo);
            let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                0
            } else if ext[1] >= ext[2] {
                1
            } else {
                2
            };
            let mid = (s + e) / 2;
            idx[s..e].select_nth_unstable_by(mid - s, |a, b| {
                cent[*a as usize][axis].total_cmp(&cent[*b as usize][axis])
            });
            let left = nodes.len();
            nodes.push(Node::default());
            nodes.push(Node::default());
            nodes[ni].first = left as u32;
            nodes[ni].count = 0;
            stack.push((left, s, mid));
            stack.push((left + 1, mid, e));
        }
        let tris: Vec<[V3; 3]> = idx.iter().map(|&i| tris[i as usize]).collect();
        // Dipole moments bottom-up: children always have larger indices than their parent.
        for ni in (0..nodes.len()).rev() {
            let nd = nodes[ni];
            if nd.count > 0 {
                let (mut n, mut c, mut area) = ([0.0; 3], [0.0; 3], 0.0);
                let ts = &tris[nd.first as usize..(nd.first + nd.count) as usize];
                for t in ts {
                    let a = mul(cross(sub(t[1], t[0]), sub(t[2], t[0])), 0.5);
                    let w = norm(a);
                    n = add(n, a);
                    c = add(c, mul(add(add(t[0], t[1]), t[2]), w / 3.0));
                    area += w;
                }
                c = if area > 0.0 {
                    mul(c, 1.0 / area)
                } else {
                    mul(add(nd.lo, nd.hi), 0.5)
                };
                let r = ts
                    .iter()
                    .flatten()
                    .map(|p| norm(sub(*p, c)))
                    .fold(0.0, f64::max);
                nodes[ni].n = n;
                nodes[ni].c = c;
                nodes[ni].r = r;
                nodes[ni].area = area;
            } else {
                let (l, r) = (nodes[nd.first as usize], nodes[nd.first as usize + 1]);
                let area = l.area + r.area;
                let c = if area > 0.0 {
                    mul(add(mul(l.c, l.area), mul(r.c, r.area)), 1.0 / area)
                } else {
                    mul(add(nd.lo, nd.hi), 0.5)
                };
                let rad = (norm(sub(l.c, c)) + l.r).max(norm(sub(r.c, c)) + r.r);
                nodes[ni].n = add(l.n, r.n);
                nodes[ni].c = c;
                nodes[ni].r = rad;
                nodes[ni].area = area;
            }
        }
        Bvh { nodes, tris }
    }

    /// Generalized winding number at `q`: about 1 inside, 0 outside, fractional near holes.
    /// `beta` trades accuracy for speed (2 is a good default).
    pub fn winding(&self, q: V3, beta: f64) -> f64 {
        let mut w = 0.0;
        let mut stack = [0u32; 96];
        let mut sp = 1;
        while sp > 0 {
            sp -= 1;
            let nd = &self.nodes[stack[sp] as usize];
            let d = sub(nd.c, q);
            let d2 = norm2(d);
            if d2 > beta * beta * nd.r * nd.r {
                w += dot(nd.n, d) / (d2 * d2.sqrt());
                continue;
            }
            if nd.count > 0 {
                for t in &self.tris[nd.first as usize..(nd.first + nd.count) as usize] {
                    w += solid_angle(t, q);
                }
            } else {
                stack[sp] = nd.first;
                stack[sp + 1] = nd.first + 1;
                sp += 2;
            }
        }
        w * INV_4PI
    }

    /// Squared distance from `q` to the nearest triangle, or `max2` if none is closer.
    pub fn dist2(&self, q: V3, max2: f64) -> f64 {
        let mut best = max2;
        let mut stack = [0u32; 96];
        let mut sp = 1;
        while sp > 0 {
            sp -= 1;
            let nd = &self.nodes[stack[sp] as usize];
            if box_dist2(nd, q) >= best {
                continue;
            }
            if nd.count > 0 {
                for t in &self.tris[nd.first as usize..(nd.first + nd.count) as usize] {
                    best = best.min(point_tri_dist2(q, t));
                }
            } else {
                let (a, b) = (nd.first, nd.first + 1);
                let (da, db) = (
                    box_dist2(&self.nodes[a as usize], q),
                    box_dist2(&self.nodes[b as usize], q),
                );
                // Visit the nearer child first (pushed last).
                let (near, far) = if da <= db { (a, b) } else { (b, a) };
                stack[sp] = far;
                stack[sp + 1] = near;
                sp += 2;
            }
        }
        best
    }
}

fn box_dist2(nd: &Node, q: V3) -> f64 {
    let mut d = 0.0;
    for k in 0..3 {
        let e = (nd.lo[k] - q[k]).max(0.0).max(q[k] - nd.hi[k]);
        d += e * e;
    }
    d
}

/// Signed solid angle of triangle `t` seen from `q` (Van Oosterom & Strackee).
#[inline]
fn solid_angle(t: &[V3; 3], q: V3) -> f64 {
    let (a, b, c) = (sub(t[0], q), sub(t[1], q), sub(t[2], q));
    let (la, lb, lc) = (norm(a), norm(b), norm(c));
    let num = dot(a, cross(b, c));
    let den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    2.0 * num.atan2(den)
}

/// Squared distance from `p` to triangle `t` (Ericson, Real-Time Collision Detection 5.1.5).
pub fn point_tri_dist2(p: V3, t: &[V3; 3]) -> f64 {
    let (a, b, c) = (t[0], t[1], t[2]);
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return norm2(ap);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return norm2(bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return norm2(sub(p, add(a, mul(ab, v))));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return norm2(cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return norm2(sub(p, add(a, mul(ac, w))));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return norm2(sub(p, add(b, mul(sub(c, b), w))));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    norm2(sub(p, add(a, add(mul(ab, v), mul(ac, w)))))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube() -> Mesh {
        let v = vec![
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 1.],
            [1., 1., 1.],
            [0., 1., 1.],
        ];
        let f = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ];
        Mesh { v, f }
    }

    #[test]
    fn winding_inside_outside() {
        let b = Bvh::new(&cube());
        assert!((b.winding([0.5, 0.5, 0.5], 2.0) - 1.0).abs() < 1e-9);
        let w = b.winding([3.0, 0.5, 0.5], 2.0);
        assert!(w.abs() < 1e-2, "{w}"); // dipole far field: about 1% at beta = 2
        assert!((b.winding([0.5, 0.5, 0.999], 2.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn winding_open_box_is_fractional_at_hole() {
        let mut m = cube();
        m.f.truncate(10); // drop the x=0 side
        let b = Bvh::new(&m);
        let w = b.winding([0.0, 0.5, 0.5], 2.0);
        assert!(w > 0.3 && w < 0.7, "{w}");
    }

    #[test]
    fn distance() {
        let b = Bvh::new(&cube());
        assert!((b.dist2([2.0, 0.5, 0.5], f64::INFINITY) - 1.0).abs() < 1e-12);
        assert!((b.dist2([0.5, 0.5, 0.25], f64::INFINITY) - 0.0625).abs() < 1e-12);
    }
}
