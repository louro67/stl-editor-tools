//! Biharmonic fairing of hole patches: vertices that are far from every input triangle were
//! invented to close a hole, so they are moved to minimise the bi-Laplacian energy with the
//! rest of the surface held fixed. Patches then continue the curvature around them instead of
//! keeping the voxel steps of the inside/outside mask.

use crate::bvh::Bvh;
use crate::mesh::Mesh;
use crate::par;

/// Compressed adjacency.
pub struct Adjacency {
    pub start: Vec<u32>,
    pub nb: Vec<u32>,
}

impl Adjacency {
    pub fn new(m: &Mesh) -> Adjacency {
        let nv = m.v.len();
        let mut deg = vec![0u32; nv + 1];
        for t in &m.f {
            for k in 0..3 {
                deg[t[k] as usize] += 2;
            }
        }
        let mut start = vec![0u32; nv + 1];
        for i in 0..nv {
            start[i + 1] = start[i] + deg[i];
        }
        let mut fill = start.clone();
        let mut nb = vec![0u32; start[nv] as usize];
        for t in &m.f {
            for k in 0..3 {
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                nb[fill[a as usize] as usize] = b;
                nb[fill[a as usize] as usize + 1] = c;
                fill[a as usize] += 2;
            }
        }
        // Sort and dedup each list in place, then compact.
        let mut out_start = vec![0u32; nv + 1];
        let mut w = 0usize;
        for i in 0..nv {
            let (s, e) = (start[i] as usize, start[i + 1] as usize);
            let mut l: Vec<u32> = nb[s..e].to_vec();
            l.sort_unstable();
            l.dedup();
            nb[w..w + l.len()].copy_from_slice(&l);
            w += l.len();
            out_start[i + 1] = w as u32;
        }
        nb.truncate(w);
        Adjacency {
            start: out_start,
            nb,
        }
    }
    #[inline]
    pub fn of(&self, i: usize) -> &[u32] {
        &self.nb[self.start[i] as usize..self.start[i + 1] as usize]
    }
}

/// Fair the vertices farther than `thr` from the input (`bvh`), grown by `grow` rings.
/// Returns how many vertices moved.
pub fn fair(m: &mut Mesh, bvh: &Bvh, thr: f64, grow: usize, max_iter: usize) -> usize {
    let nv = m.v.len();
    let far: Vec<bool> = par::map(nv, |i| bvh.dist2(m.v[i], thr * thr) >= thr * thr);
    let adj = Adjacency::new(m);
    // 0 = fixed and untouched, 1 = free, 2 = first fixed ring, 3 = second fixed ring.
    let mut role = vec![0u8; nv];
    let mut front: Vec<u32> = (0..nv as u32).filter(|&i| far[i as usize]).collect();
    if front.is_empty() {
        return 0;
    }
    for &i in &front {
        role[i as usize] = 1;
    }
    for _ in 0..grow {
        let mut next = Vec::new();
        for &i in &front {
            for &j in adj.of(i as usize) {
                if role[j as usize] == 0 {
                    role[j as usize] = 1;
                    next.push(j);
                }
            }
        }
        front = next;
    }
    for ring in [2u8, 3] {
        let src: Vec<u32> = (0..nv as u32)
            .filter(|&i| role[i as usize] == ring - 1)
            .collect();
        for i in src {
            for &j in adj.of(i as usize) {
                if role[j as usize] == 0 {
                    role[j as usize] = ring;
                }
            }
        }
    }
    // Local numbering: free vertices first, then rings.
    let mut order: Vec<u32> = (0..nv as u32).filter(|&i| role[i as usize] == 1).collect();
    let nfree = order.len();
    order.extend((0..nv as u32).filter(|&i| role[i as usize] == 2));
    let nkr = order.len(); // vertices where K x is needed
    order.extend((0..nv as u32).filter(|&i| role[i as usize] == 3));
    let mut local = vec![u32::MAX; nv];
    for (l, &g) in order.iter().enumerate() {
        local[g as usize] = l as u32;
    }
    let nbr: Vec<Vec<u32>> = order[..nkr]
        .iter()
        .map(|&g| {
            adj.of(g as usize)
                .iter()
                .map(|&j| local[j as usize])
                .collect()
        })
        .collect();
    let n = order.len();
    let x0: Vec<[f64; 3]> = order.iter().map(|&g| m.v[g as usize]).collect();

    // K x over the first `nkr` local vertices, for x given on all `n`; neighbours outside the
    // work set never occur for these vertices.
    let lap = |x: &[[f64; 3]], upto: usize| -> Vec<[f64; 3]> {
        par::map(upto, |i| {
            let l = &nbr[i];
            let d = l.len() as f64;
            let mut s = [d * x[i][0], d * x[i][1], d * x[i][2]];
            for &j in l {
                let xj = x[j as usize];
                s[0] -= xj[0];
                s[1] -= xj[1];
                s[2] -= xj[2];
            }
            s
        })
    };
    // M p = (K K p) on free vertices, p zero outside them.
    let apply = |p: &[[f64; 3]]| -> Vec<[f64; 3]> {
        let mut full = vec![[0.0; 3]; n];
        full[..nfree].copy_from_slice(p);
        let mut y = lap(&full, nkr);
        y.resize(n, [0.0; 3]);
        lap(&y, nfree)
    };
    let mut y0 = lap(&x0, nkr);
    y0.resize(n, [0.0; 3]);
    let ky0 = lap(&y0, nfree);
    let b: Vec<[f64; 3]> = ky0.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
    let dinv: Vec<f64> = (0..nfree)
        .map(|i| {
            let d = nbr[i].len() as f64;
            1.0 / (d * d + d)
        })
        .collect();

    // Preconditioned conjugate gradient, three right-hand sides at once.
    let mut x = vec![[0.0; 3]; nfree];
    let mut r = b.clone();
    let mut z: Vec<[f64; 3]> = (0..nfree)
        .map(|i| [r[i][0] * dinv[i], r[i][1] * dinv[i], r[i][2] * dinv[i]])
        .collect();
    let mut p = z.clone();
    let dot3 = |a: &[[f64; 3]], b: &[[f64; 3]]| -> [f64; 3] {
        let mut s = [0.0; 3];
        for (u, v) in a.iter().zip(b) {
            for k in 0..3 {
                s[k] += u[k] * v[k];
            }
        }
        s
    };
    let bnorm = dot3(&b, &b);
    let mut rz = dot3(&r, &z);
    for _ in 0..max_iter {
        let ap = apply(&p);
        let pap = dot3(&p, &ap);
        let mut alpha = [0.0; 3];
        for k in 0..3 {
            alpha[k] = if pap[k] > 0.0 { rz[k] / pap[k] } else { 0.0 };
        }
        for i in 0..nfree {
            for k in 0..3 {
                x[i][k] += alpha[k] * p[i][k];
                r[i][k] -= alpha[k] * ap[i][k];
            }
        }
        let rr = dot3(&r, &r);
        if (0..3).all(|k| rr[k] <= 1e-10 * bnorm[k].max(1e-30)) {
            break;
        }
        for i in 0..nfree {
            for k in 0..3 {
                z[i][k] = r[i][k] * dinv[i];
            }
        }
        let rz_new = dot3(&r, &z);
        for k in 0..3 {
            let beta = if rz[k] > 0.0 { rz_new[k] / rz[k] } else { 0.0 };
            for i in 0..nfree {
                p[i][k] = z[i][k] + beta * p[i][k];
            }
        }
        rz = rz_new;
    }
    for i in 0..nfree {
        let g = order[i] as usize;
        for k in 0..3 {
            m.v[g][k] += x[i][k];
        }
    }
    nfree
}
