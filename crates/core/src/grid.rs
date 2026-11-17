//! The voxel grid: an adaptive inside/outside mask from the winding number, then
//! connected-component cleanup (keep the largest solid, fill enclosed voids).

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use crate::bvh::Bvh;
use crate::mesh::V3;
use crate::par;
use crate::pipeline::Monitor;

/// Block edge in nodes for the adaptive winding-number pass.
pub const B: usize = 16;

#[derive(Clone, Copy, Debug)]
pub struct Grid {
    pub n: [usize; 3],
    pub origin: V3,
    pub h: f64,
}

impl Grid {
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.n[0] * (j + self.n[1] * k)
    }
    #[inline]
    pub fn pos(&self, i: usize, j: usize, k: usize) -> V3 {
        [
            self.origin[0] + self.h * i as f64,
            self.origin[1] + self.h * j as f64,
            self.origin[2] + self.h * k as f64,
        ]
    }
    pub fn len(&self) -> usize {
        self.n[0] * self.n[1] * self.n[2]
    }
}

/// Inside (1) / outside (0) per grid node: winding number >= `iso`.
///
/// The winding number is smooth (harmonic) away from the input, so it is sampled on a coarse
/// lattice and refined octree-style: a box is split only when it is within `near` of the input
/// or its corners disagree or come close to `iso`; every other box takes its corners' side.
/// Only a thin shell around the surface (and around hole patches) is evaluated per node.
/// Returns the mask and the number of winding-number evaluations.
pub fn inside_mask(
    bvh: &Bvh,
    g: &Grid,
    iso: f64,
    near: f64,
    mon: &dyn Monitor,
) -> Result<(Vec<u8>, usize), ()> {
    let nb = [g.n[0].div_ceil(B), g.n[1].div_ceil(B), g.n[2].div_ceil(B)];
    let cn = [nb[0] + 1, nb[1] + 1, nb[2] + 1];
    let cnode = |c: usize, ax: usize| (c * B).min(g.n[ax] - 1);

    // Coarse lattice, one z layer per task.
    let coarse: Vec<Vec<f64>> = par::map(cn[2], |ck| {
        let mut out = Vec::with_capacity(cn[0] * cn[1]);
        for cj in 0..cn[1] {
            for ci in 0..cn[0] {
                out.push(bvh.winding(g.pos(cnode(ci, 0), cnode(cj, 1), cnode(ck, 2)), BETA));
            }
        }
        out
    });
    if mon.cancelled() {
        return Err(());
    }
    let cw = |ci: usize, cj: usize, ck: usize| coarse[ck][ci + cn[0] * cj];
    let mut evals = cn[0] * cn[1] * cn[2];

    let mut mask = vec![0u8; g.len()];
    let slab = g.n[0] * g.n[1] * B;
    let done = AtomicUsize::new(0);
    let fine = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    par::chunks_mut(&mut mask, slab, |bk, chunk| {
        if stop.load(Ordering::Relaxed) || mon.cancelled() {
            stop.store(true, Ordering::Relaxed);
            return;
        }
        let mut r = Refiner {
            bvh,
            g,
            iso,
            near,
            chunk,
            k0: bk * B,
            evals: 0,
        };
        for bj in 0..nb[1] {
            for bi in 0..nb[0] {
                let mut w = [0.0; 8];
                for (c, x) in w.iter_mut().enumerate() {
                    *x = cw(bi + (c & 1), bj + ((c >> 1) & 1), bk + ((c >> 2) & 1));
                }
                r.box_(bi * B, bj * B, bk * B, B, w);
            }
        }
        fine.fetch_add(r.evals, Ordering::Relaxed);
        let d = done.fetch_add(1, Ordering::Relaxed) + 1;
        mon.progress(d as f32 / nb[2] as f32, "winding number");
    });
    if stop.load(Ordering::Relaxed) {
        return Err(());
    }
    evals += fine.load(Ordering::Relaxed);
    // The outermost layer is always outside so the surface closes.
    clear_border(&mut mask, g);
    Ok((mask, evals))
}

const BETA: f64 = 2.0;

/// Octree refinement inside one z slab of the mask.
struct Refiner<'a> {
    bvh: &'a Bvh,
    g: &'a Grid,
    iso: f64,
    near: f64,
    chunk: &'a mut [u8],
    k0: usize,
    evals: usize,
}

impl Refiner<'_> {
    fn w(&mut self, i: usize, j: usize, k: usize) -> f64 {
        self.evals += 1;
        let g = self.g;
        self.bvh.winding(
            g.pos(i.min(g.n[0] - 1), j.min(g.n[1] - 1), k.min(g.n[2] - 1)),
            BETA,
        )
    }

    fn set(&mut self, i: usize, j: usize, k: usize, inside: bool) {
        let g = self.g;
        if i < g.n[0] && j < g.n[1] && k < g.n[2] {
            self.chunk[i + g.n[0] * (j + g.n[1] * (k - self.k0))] = inside as u8;
        }
    }

    /// The box of `s` nodes per axis at (i0, j0, k0), with winding numbers `w` at its corners
    /// (corner bit x = 1, y = 2, z = 4; the far corners lie at +s).
    fn box_(&mut self, i0: usize, j0: usize, k0: usize, s: usize, w: [f64; 8]) {
        let g = self.g;
        if i0 >= g.n[0] || j0 >= g.n[1] || k0 >= g.n[2] {
            return;
        }
        let ins = w.iter().filter(|&&x| x >= self.iso).count();
        let doubtful = (ins != 0 && ins != 8) || w.iter().any(|&x| (x - self.iso).abs() < 0.25);
        let near = || {
            let hs = 0.5 * s as f64 * g.h;
            let c = g.pos(i0, j0, k0);
            let reach = 3f64.sqrt() * hs + self.near;
            self.bvh
                .dist2([c[0] + hs, c[1] + hs, c[2] + hs], reach * reach)
                < reach * reach
        };
        if !doubtful && !near() {
            for k in k0..k0 + s {
                for j in j0..j0 + s {
                    for i in i0..i0 + s {
                        self.set(i, j, k, ins == 8);
                    }
                }
            }
            return;
        }
        if s == 1 {
            self.set(i0, j0, k0, w[0] >= self.iso);
            return;
        }
        if s == 2 {
            for c in 0..8usize {
                let (i, j, k) = (i0 + (c & 1), j0 + ((c >> 1) & 1), k0 + ((c >> 2) & 1));
                let x = if c == 0 { w[0] } else { self.w(i, j, k) };
                self.set(i, j, k, x >= self.iso);
            }
            return;
        }
        // Split in eight: winding numbers on the 3 x 3 x 3 lattice, corners reused.
        let h = s / 2;
        let mut l = [0.0; 27];
        for c in 0..27usize {
            let (a, b, d) = (c % 3, (c / 3) % 3, c / 9);
            l[c] = if a != 1 && b != 1 && d != 1 {
                w[(a / 2) | ((b / 2) << 1) | ((d / 2) << 2)]
            } else {
                self.w(i0 + a * h, j0 + b * h, k0 + d * h)
            };
        }
        for child in 0..8usize {
            let (cx, cy, cz) = (child & 1, (child >> 1) & 1, (child >> 2) & 1);
            let mut cw = [0.0; 8];
            for (c, x) in cw.iter_mut().enumerate() {
                let (a, b, d) = (cx + (c & 1), cy + ((c >> 1) & 1), cz + ((c >> 2) & 1));
                *x = l[a + 3 * b + 9 * d];
            }
            self.box_(i0 + cx * h, j0 + cy * h, k0 + cz * h, h, cw);
        }
    }
}

fn clear_border(mask: &mut [u8], g: &Grid) {
    let [n0, n1, n2] = g.n;
    for k in 0..n2 {
        for j in 0..n1 {
            if k == 0 || j == 0 || k == n2 - 1 || j == n1 - 1 {
                let r = g.idx(0, j, k);
                mask[r..r + n0].fill(0);
            } else {
                mask[g.idx(0, j, k)] = 0;
                mask[g.idx(n0 - 1, j, k)] = 0;
            }
        }
    }
}

/// Runs of equal mask value along x, labelled by 6-connected component.
struct Runs {
    /// Row r = j + n1 * k holds runs `start[r]..start[r + 1]`.
    start: Vec<usize>,
    x0: Vec<u32>,
    x1: Vec<u32>,
    parent: Vec<u32>,
}

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
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        p[hi as usize] = lo;
    }
}

fn label(mask: &[u8], g: &Grid, value: u8) -> Runs {
    let [n0, n1, n2] = g.n;
    let rows = n1 * n2;
    let mut start = Vec::with_capacity(rows + 1);
    let (mut x0, mut x1) = (Vec::new(), Vec::new());
    for r in 0..rows {
        start.push(x0.len());
        let row = &mask[r * n0..(r + 1) * n0];
        let mut i = 0;
        while i < n0 {
            if row[i] == value {
                let s = i;
                while i < n0 && row[i] == value {
                    i += 1;
                }
                x0.push(s as u32);
                x1.push(i as u32);
            } else {
                i += 1;
            }
        }
    }
    start.push(x0.len());
    let mut parent: Vec<u32> = (0..x0.len() as u32).collect();
    for k in 0..n2 {
        for j in 0..n1 {
            let r = j + n1 * k;
            let link = |q: usize, parent: &mut Vec<u32>| {
                let (mut a, mut b) = (start[r], start[q]);
                while a < start[r + 1] && b < start[q + 1] {
                    if x0[a] < x1[b] && x0[b] < x1[a] {
                        union(parent, a as u32, b as u32);
                    }
                    if x1[a] < x1[b] {
                        a += 1;
                    } else {
                        b += 1;
                    }
                }
            };
            if j > 0 {
                link(r - 1, &mut parent);
            }
            if k > 0 {
                link(r - n1, &mut parent);
            }
        }
    }
    Runs {
        start,
        x0,
        x1,
        parent,
    }
}

/// Keep the largest inside component and fill outside pockets that don't reach the border.
/// Returns (inside voxels dropped, void voxels filled, inside components found).
pub fn clean(mask: &mut [u8], g: &Grid) -> (usize, usize, usize) {
    let n0 = g.n[0];
    let n1 = g.n[1];
    let mut runs = label(mask, g, 1);
    let nr = runs.x0.len();
    let mut size = vec![0usize; nr];
    for a in 0..nr {
        let root = find(&mut runs.parent, a as u32) as usize;
        size[root] += (runs.x1[a] - runs.x0[a]) as usize;
    }
    let comps = (0..nr).filter(|&a| runs.parent[a] as usize == a).count();
    let keep = (0..nr).max_by_key(|&a| size[a]).unwrap_or(0) as u32;
    let mut dropped = 0;
    for r in 0..runs.start.len() - 1 {
        for a in runs.start[r]..runs.start[r + 1] {
            if find(&mut runs.parent, a as u32) != keep {
                mask[r * n0 + runs.x0[a] as usize..r * n0 + runs.x1[a] as usize].fill(0);
                dropped += (runs.x1[a] - runs.x0[a]) as usize;
            }
        }
    }

    let mut out = label(mask, g, 0);
    let no = out.x0.len();
    let mut border = vec![false; no];
    for r in 0..out.start.len() - 1 {
        let (j, k) = (r % n1, r / n1);
        for a in out.start[r]..out.start[r + 1] {
            if j == 0
                || k == 0
                || j == n1 - 1
                || k == g.n[2] - 1
                || out.x0[a] == 0
                || out.x1[a] as usize == n0
            {
                let root = find(&mut out.parent, a as u32) as usize;
                border[root] = true;
            }
        }
    }
    let mut filled = 0;
    for r in 0..out.start.len() - 1 {
        for a in out.start[r]..out.start[r + 1] {
            let root = find(&mut out.parent, a as u32) as usize;
            if !border[root] {
                mask[r * n0 + out.x0[a] as usize..r * n0 + out.x1[a] as usize].fill(1);
                filled += (out.x1[a] - out.x0[a]) as usize;
            }
        }
    }
    (dropped, filled, comps)
}
