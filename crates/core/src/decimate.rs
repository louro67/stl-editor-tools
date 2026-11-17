//! Quadric error metric edge collapse (Garland & Heckbert 1997) that keeps a closed manifold
//! closed and manifold: collapses must pass the link condition and may not flip or flatten faces.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::mesh::{add, cross, dot, mul, norm, norm2, sub, Mesh, V3};

type Quadric = [f64; 10];

/// Triangles worse than this shape quality make a collapse proportionally more expensive.
const QUALITY: f64 = 0.3; // a2 ab ac ad b2 bc bd c2 cd d2

fn plane_quadric(n: V3, d: f64, w: f64) -> Quadric {
    let [a, b, c] = n;
    [
        a * a * w,
        a * b * w,
        a * c * w,
        a * d * w,
        b * b * w,
        b * c * w,
        b * d * w,
        c * c * w,
        c * d * w,
        d * d * w,
    ]
}

fn qadd(a: &Quadric, b: &Quadric) -> Quadric {
    let mut r = [0.0; 10];
    for i in 0..10 {
        r[i] = a[i] + b[i];
    }
    r
}

fn qeval(q: &Quadric, p: V3) -> f64 {
    let [x, y, z] = p;
    q[0] * x * x
        + 2.0 * q[1] * x * y
        + 2.0 * q[2] * x * z
        + 2.0 * q[3] * x
        + q[4] * y * y
        + 2.0 * q[5] * y * z
        + 2.0 * q[6] * y
        + q[7] * z * z
        + 2.0 * q[8] * z
        + q[9]
}

/// Point minimising `q`, if the system is well conditioned.
fn qsolve(q: &Quadric) -> Option<V3> {
    let a = [[q[0], q[1], q[2]], [q[1], q[4], q[5]], [q[2], q[5], q[7]]];
    let b = [-q[3], -q[6], -q[8]];
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let scale = (q[0] + q[4] + q[7]).powi(3);
    if det.abs() <= 1e-9 * scale.abs().max(1e-300) {
        return None;
    }
    let inv = 1.0 / det;
    let col = |c: usize| {
        let mut m = a;
        for r in 0..3 {
            m[r][c] = b[r];
        }
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    Some([col(0) * inv, col(1) * inv, col(2) * inv])
}

struct Entry {
    cost: f64,
    u: u32,
    v: u32,
    su: u32,
    sv: u32,
}
impl PartialEq for Entry {
    fn eq(&self, o: &Self) -> bool {
        self.cost == o.cost
    }
}
impl Eq for Entry {}
impl PartialOrd for Entry {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Entry {
    fn cmp(&self, o: &Self) -> Ordering {
        o.cost.total_cmp(&self.cost) // min-heap
    }
}

struct State {
    pos: Vec<V3>,
    faces: Vec<[u32; 3]>,
    alive: Vec<bool>,
    vf: Vec<Vec<u32>>,
    q: Vec<Quadric>,
    stamp: Vec<u32>,
    vdead: Vec<bool>,
    lambda: f64,
}

impl State {
    fn neighbors(&self, u: u32, out: &mut Vec<u32>) {
        out.clear();
        for &f in &self.vf[u as usize] {
            for &w in &self.faces[f as usize] {
                if w != u {
                    out.push(w);
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }

    fn candidate(&self, u: u32, v: u32) -> (f64, V3) {
        let q = qadd(&self.q[u as usize], &self.q[v as usize]);
        let (pu, pv) = (self.pos[u as usize], self.pos[v as usize]);
        let mid = mul(add(pu, pv), 0.5);
        let len2 = norm2(sub(pu, pv));
        let mut best = (qeval(&q, mid), mid);
        if let Some(x) = qsolve(&q) {
            // Stay near the edge: far optima come from nearly flat quadrics.
            if norm2(sub(x, mid)) <= len2 {
                let c = qeval(&q, x);
                if c < best.0 {
                    best = (c, x);
                }
            }
        }
        for p in [pu, pv] {
            let c = qeval(&q, p);
            if c < best.0 {
                best = (c, p);
            }
        }
        (best.0.max(0.0) + self.lambda * len2, best.1)
    }

    fn push(
        &self,
        heap: &mut BinaryHeap<Entry>,
        u: u32,
        v: u32,
        nu: &mut Vec<u32>,
        nv: &mut Vec<u32>,
    ) {
        if let Some((cost, _)) = self.evaluate(u, v, nu, nv) {
            heap.push(Entry {
                cost,
                u,
                v,
                su: self.stamp[u as usize],
                sv: self.stamp[v as usize],
            });
        }
    }

    /// Collapse v into u at x. `nu` holds u's neighbours from `valid`.
    fn collapse(&mut self, u: u32, v: u32, x: V3, nu: &[u32], faces_left: &mut usize) {
        let vfaces = std::mem::take(&mut self.vf[v as usize]);
        for f in vfaces {
            let t = &mut self.faces[f as usize];
            if t.contains(&u) {
                self.alive[f as usize] = false;
                *faces_left -= 1;
            } else {
                for k in t.iter_mut() {
                    if *k == v {
                        *k = u;
                    }
                }
                self.vf[u as usize].push(f);
            }
        }
        let alive = &self.alive;
        self.vf[u as usize].retain(|&f| alive[f as usize]);
        // The two opposite vertices lost a face each; they are among u's neighbours.
        for &w in nu {
            let l = &mut self.vf[w as usize];
            if l.iter().any(|&f| !alive[f as usize]) {
                l.retain(|&f| alive[f as usize]);
            }
        }
        self.pos[u as usize] = x;
        self.q[u as usize] = qadd(&self.q[u as usize], &self.q[v as usize]);
        self.vdead[v as usize] = true;
        self.stamp[u as usize] += 1;
    }

    /// Whether collapsing v into u at `x` keeps the mesh manifold and unflipped.
    fn valid(&self, u: u32, v: u32, x: V3, nu: &mut Vec<u32>, nv: &mut Vec<u32>) -> Option<f64> {
        self.neighbors(u, nu);
        self.neighbors(v, nv);
        let shared = self.vf[u as usize]
            .iter()
            .filter(|&&f| self.faces[f as usize].contains(&v))
            .count();
        if shared != 2 {
            return None;
        }
        // Link condition: common neighbours are exactly the two opposite vertices.
        let mut common = 0;
        let (mut i, mut j) = (0, 0);
        while i < nu.len() && j < nv.len() {
            match nu[i].cmp(&nv[j]) {
                Ordering::Less => i += 1,
                Ordering::Greater => j += 1,
                Ordering::Equal => {
                    common += 1;
                    i += 1;
                    j += 1;
                }
            }
        }
        if common != 2 || nu.len() + nv.len() < 8 {
            return None;
        }
        let mut minq: f64 = 1.0;
        for (a, other) in [(u, v), (v, u)] {
            for &f in &self.vf[a as usize] {
                let t = self.faces[f as usize];
                if t.contains(&other) {
                    continue;
                }
                let p = [
                    self.pos[t[0] as usize],
                    self.pos[t[1] as usize],
                    self.pos[t[2] as usize],
                ];
                let n0 = cross(sub(p[1], p[0]), sub(p[2], p[0]));
                let mut q = p;
                for k in 0..3 {
                    if t[k] == a {
                        q[k] = x;
                    }
                }
                let n1 = cross(sub(q[1], q[0]), sub(q[2], q[0]));
                let (l0, l1) = (norm(n0), norm(n1));
                if l1 <= 1e-12 * l0.max(1e-300) || dot(n0, n1) < 0.5 * l0 * l1 {
                    return None;
                }
                // Shape quality: 1 for equilateral, 0 for degenerate (2 sqrt(3) area / sum e^2).
                let e2 = norm2(sub(q[1], q[0])) + norm2(sub(q[2], q[1])) + norm2(sub(q[0], q[2]));
                minq = minq.min(2.0 * 3f64.sqrt() * l1 / e2.max(1e-300));
            }
        }
        Some(minq)
    }

    /// Full cost of collapsing (u, v): quadric error plus a length term, raised for collapses
    /// that would make poorly shaped triangles. None when the collapse is not allowed.
    fn evaluate(&self, u: u32, v: u32, nu: &mut Vec<u32>, nv: &mut Vec<u32>) -> Option<(f64, V3)> {
        let (cost, x) = self.candidate(u, v);
        let q = self.valid(u, v, x, nu, nv)?;
        if q < 0.02 {
            return None;
        }
        let penalty = if q < QUALITY { QUALITY / q } else { 1.0 };
        Some((cost * penalty, x))
    }
}

/// Reduce `m` (closed, manifold) to about `target` faces.
pub fn decimate(m: &mut Mesh, target: usize, cancelled: &dyn Fn() -> bool) -> Result<(), ()> {
    if m.f.len() <= target {
        return Ok(());
    }
    let nv = m.v.len();
    let mut vf: Vec<Vec<u32>> = vec![Vec::with_capacity(6); nv];
    let mut q = vec![[0.0; 10]; nv];
    let mut area_sum = 0.0;
    for (i, t) in m.f.iter().enumerate() {
        let (a, b, c) = (m.v[t[0] as usize], m.v[t[1] as usize], m.v[t[2] as usize]);
        let n = cross(sub(b, a), sub(c, a));
        let l = norm(n);
        area_sum += 0.5 * l;
        if l > 0.0 {
            let nn = mul(n, 1.0 / l);
            let pq = plane_quadric(nn, -dot(nn, a), 0.5 * l);
            for &k in t {
                q[k as usize] = qadd(&q[k as usize], &pq);
            }
        }
        for &k in t {
            vf[k as usize].push(i as u32);
        }
    }
    let lambda = 1e-3 * area_sum / m.f.len() as f64;
    let mut s = State {
        pos: std::mem::take(&mut m.v),
        faces: std::mem::take(&mut m.f),
        alive: Vec::new(),
        vf,
        q,
        stamp: vec![0; nv],
        vdead: vec![false; nv],
        lambda,
    };
    s.alive = vec![true; s.faces.len()];
    let mut faces_left = s.faces.len();
    let (mut nu, mut nvb) = (Vec::new(), Vec::new());
    let mut steps = 0u64;

    // Phase 1, threshold sweeps: walk the faces in memory order (spatially coherent) and
    // collapse each face's cheapest edge when its cost is below a threshold picked from a
    // sample of current costs. Close to greedy order, without a queue over millions of edges.
    let switch = target.saturating_mul(2).max(target + 1000);
    let mut rng = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    let nf = s.faces.len();
    let cheapest = |s: &State, t: [u32; 3]| -> (f64, u32, u32) {
        let mut best = (f64::INFINITY, 0, 0);
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let (c, _) = s.candidate(a, b);
            if c < best.0 {
                best = (c, a, b);
            }
        }
        best
    };
    while faces_left > switch {
        let mut sample: Vec<f64> = Vec::with_capacity(4096);
        let mut tries = 0;
        while sample.len() < 4096 && tries < 1 << 20 {
            tries += 1;
            let f = (next() % nf as u64) as usize;
            if s.alive[f] {
                sample.push(cheapest(&s, s.faces[f]).0);
            }
        }
        if sample.is_empty() {
            break;
        }
        sample.sort_by(f64::total_cmp);
        let frac = ((faces_left - switch) as f64 / faces_left as f64).min(0.3);
        let tau = sample[((sample.len() - 1) as f64 * frac) as usize];
        let before = faces_left;
        for f in 0..nf {
            if faces_left <= switch {
                break;
            }
            steps += 1;
            if steps % 65536 == 0 && cancelled() {
                return Err(());
            }
            if !s.alive[f] {
                continue;
            }
            let (c, u, v) = cheapest(&s, s.faces[f]);
            if c > tau {
                continue;
            }
            if let Some((c, x)) = s.evaluate(u, v, &mut nu, &mut nvb) {
                if c <= tau {
                    s.collapse(u, v, x, &nu, &mut faces_left);
                }
            }
        }
        if faces_left == before {
            break;
        }
    }

    // Phase 2: exact greedy order with a heap.
    let mut heap = BinaryHeap::with_capacity(faces_left * 3 / 2);
    for i in 0..s.faces.len() {
        if !s.alive[i] {
            continue;
        }
        let t = s.faces[i];
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if a < b {
                s.push(&mut heap, a, b, &mut nu, &mut nvb);
            }
        }
    }
    let mut around = Vec::new();
    let (mut tu, mut tv) = (Vec::new(), Vec::new());
    while faces_left > target {
        let Some(e) = heap.pop() else { break };
        steps += 1;
        if steps % 65536 == 0 && cancelled() {
            return Err(());
        }
        let (u, v) = (e.u, e.v);
        if s.vdead[u as usize]
            || s.vdead[v as usize]
            || s.stamp[u as usize] != e.su
            || s.stamp[v as usize] != e.sv
        {
            continue;
        }
        let Some((_, x)) = s.evaluate(u, v, &mut nu, &mut nvb) else {
            continue;
        };
        s.collapse(u, v, x, &nu, &mut faces_left);
        s.neighbors(u, &mut around);
        for &w in &around {
            s.push(&mut heap, u, w, &mut tu, &mut tv);
        }
    }
    let mut out = Mesh {
        v: s.pos,
        f: Vec::with_capacity(faces_left),
    };
    for (i, t) in s.faces.iter().enumerate() {
        if s.alive[i] {
            out.f.push(*t);
        }
    }
    out.compact();
    *m = out;
    Ok(())
}
