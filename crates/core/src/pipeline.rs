//! The whole pipeline, parameters and report.

use serde::Serialize;

use crate::bvh::Bvh;
use crate::check::{verify, Verification};
use crate::grid::{self, Grid};
use crate::mesh::{sub, Mesh, V3};
use crate::{contour, decimate, fair};

/// Tuning knobs. Lengths are in output millimetres (after `scale`).
#[derive(Clone, Debug, Serialize)]
pub struct Params {
    /// Uniform scale applied to the input first (0.1 = 10%).
    pub scale: f64,
    /// Grid spacing in mm. 0 or less picks the finest spacing that fits `max_voxels`.
    pub voxel: f64,
    /// Offset the surface outward by this many mm so thin parts survive slicing.
    pub min_feature: f64,
    /// Winding-number threshold for "inside" (0..1). Lower fills more.
    pub iso: f64,
    /// Triangle budget of the output.
    pub target_faces: usize,
    /// Vertices farther than `fair_tau` voxels from the input count as hole patches and are
    /// faired. 0 turns fairing off.
    pub fair_tau: f64,
    /// Extra rings of vertices faired around each patch.
    pub fair_grow: usize,
    /// Refuse grids with more nodes than this.
    pub max_voxels: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            scale: 1.0,
            voxel: 0.125,
            min_feature: 0.0,
            iso: 0.5,
            target_faces: 400_000,
            fair_tau: 0.5,
            fair_grow: 2,
            max_voxels: 400_000_000,
        }
    }
}

/// Progress, logging and cancellation. Calls may come from worker threads.
pub trait Monitor: Sync {
    fn progress(&self, _fraction: f32, _stage: &str) {}
    fn log(&self, _message: &str) {}
    fn cancelled(&self) -> bool {
        false
    }
}

pub struct NoMonitor;
impl Monitor for NoMonitor {}

#[derive(Debug, Clone)]
pub enum Error {
    InvalidParam(String),
    GridTooLarge {
        dims: [usize; 3],
        voxels: usize,
        max: usize,
        suggest: f64,
    },
    Cancelled,
    Failed(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Error::InvalidParam(s) => write!(f, "{s}"),
            Error::GridTooLarge { dims, voxels, max, suggest } => write!(
                f,
                "grid {dims:?} = {}M voxels exceeds the limit of {}M; use a voxel of {suggest:.3} mm or more, or raise the limit",
                voxels / 1_000_000,
                max / 1_000_000
            ),
            Error::Cancelled => write!(f, "cancelled"),
            Error::Failed(s) => write!(f, "{s}"),
        }
    }
}

impl std::error::Error for Error {}

#[derive(Clone, Debug, Serialize)]
pub struct Report {
    pub ok: bool,
    pub triangles: usize,
    pub watertight: bool,
    pub is_volume: bool,
    pub bodies: usize,
    pub size_mm: [f64; 3],
    pub volume_mm3: f64,
    pub grid: [usize; 3],
    pub voxel_mm: f64,
    pub winding_evaluations: usize,
    pub components_found: usize,
    pub components_dropped_voxels: usize,
    pub voids_filled_voxels: usize,
    pub contour_triangles: usize,
    pub faired_vertices: usize,
    pub params: Params,
}

fn band(h: f64, min_feature: f64) -> usize {
    ((min_feature / h).ceil() as usize).min(1 << 20) + 2
}

/// Grid size for a model of `size` mm at spacing `h`.
pub fn grid_dims(size: V3, h: f64, min_feature: f64) -> [usize; 3] {
    let pad = band(h, min_feature) + 2;
    size.map(|s| ((s / h).ceil() as usize).saturating_add(2 * pad + 1))
}

fn nodes(d: [usize; 3]) -> f64 {
    d[0] as f64 * d[1] as f64 * d[2] as f64
}

/// The finest spacing whose grid has at most `budget` nodes.
pub fn auto_voxel(size: V3, budget: usize, min_feature: f64) -> f64 {
    let vol = size[0].max(1e-9) * size[1].max(1e-9) * size[2].max(1e-9);
    let mut h = (vol / budget as f64).cbrt().max(1e-6);
    while nodes(grid_dims(size, h, min_feature)) > budget as f64 {
        h *= 1.02;
    }
    (h * 1000.0).ceil() / 1000.0
}

fn check(p: &Params) -> Result<(), Error> {
    let bad = |s: &str| Err(Error::InvalidParam(s.into()));
    if !(p.scale.is_finite() && p.scale > 0.0) {
        return bad("scale must be a number greater than 0");
    }
    if !p.voxel.is_finite() {
        return bad("voxel must be a finite number (0 for automatic)");
    }
    if !(p.min_feature.is_finite() && p.min_feature >= 0.0) {
        return bad("min_feature must be 0 or more");
    }
    if !(p.iso > 0.0 && p.iso < 1.0) {
        return bad("iso must be between 0 and 1");
    }
    if p.target_faces < 100 {
        return bad("target_faces must be at least 100");
    }
    if !(p.fair_tau.is_finite() && p.fair_tau >= 0.0) {
        return bad("fair_tau must be 0 or more");
    }
    if p.max_voxels < 1000 {
        return bad("max_voxels must be at least 1000");
    }
    Ok(())
}

/// Rebuild `input` as one watertight solid.
pub fn solidify(input: &Mesh, p: &Params, mon: &dyn Monitor) -> Result<(Mesh, Report), Error> {
    check(p)?;
    if input.f.is_empty() {
        return Err(Error::Failed("the mesh has no triangles".into()));
    }
    let mut m = input.clone();
    m.scale(p.scale);
    let (lo, hi) = m.bounds();
    let size = sub(hi, lo);
    let h = if p.voxel > 0.0 {
        p.voxel
    } else {
        auto_voxel(size, p.max_voxels, p.min_feature)
    };
    let k_band = band(h, p.min_feature);
    let pad = (k_band + 2) as f64;
    let dims = grid_dims(size, h, p.min_feature);
    // In f64: the product overflows usize on 32-bit WebAssembly for absurd voxel sizes.
    let count = nodes(dims);
    if count > p.max_voxels as f64 {
        let voxels = count.min(usize::MAX as f64) as usize;
        let suggest = auto_voxel(size, p.max_voxels, p.min_feature);
        return Err(Error::GridTooLarge {
            dims,
            voxels,
            max: p.max_voxels,
            suggest,
        });
    }
    let voxels = count as usize;
    let g = Grid {
        n: dims,
        origin: [lo[0] - pad * h, lo[1] - pad * h, lo[2] - pad * h],
        h,
    };
    mon.log(&format!(
        "input {} triangles; grid {:?} = {:.1}M nodes at {:.3} mm",
        m.f.len(),
        dims,
        voxels as f64 / 1e6,
        h
    ));

    mon.progress(0.0, "indexing");
    let bvh = Bvh::new(&m);
    let (mut mask, evals) =
        grid::inside_mask(&bvh, &g, p.iso, 0.5 * h, mon).map_err(|_| Error::Cancelled)?;
    mon.log(&format!(
        "winding number: {:.1}M evaluations ({:.0}% of the grid)",
        evals as f64 / 1e6,
        100.0 * evals as f64 / voxels as f64
    ));

    let (dropped, filled, comps) = grid::clean(&mut mask, &g);
    mon.log(&format!("{comps} solid components; kept the largest, dropped {dropped} voxels, filled {filled} void voxels"));
    if mask.iter().all(|&x| x == 0) {
        return Err(Error::Failed(
            "nothing is inside; the mesh may be too small for this voxel size".into(),
        ));
    }
    if mon.cancelled() {
        return Err(Error::Cancelled);
    }

    mon.progress(0.0, "signed distance");
    let far = (k_band + 1) as f64 * h;
    let field = contour::field(&bvh, &g, &mask, far, p.min_feature);
    mon.progress(0.0, "contouring");
    let mut out = contour::march(&field);
    drop(field);
    drop(mask);
    let contour_triangles = out.f.len();
    mon.log(&format!(
        "marching tetrahedra: {contour_triangles} triangles"
    ));
    if mon.cancelled() {
        return Err(Error::Cancelled);
    }

    let mut faired = 0;
    if p.fair_tau > 0.0 {
        mon.progress(0.0, "fairing hole patches");
        faired = fair::fair(
            &mut out,
            &bvh,
            p.min_feature + p.fair_tau * h,
            p.fair_grow,
            400,
        );
        mon.log(&format!(
            "faired {faired} hole-patch vertices ({:.1}%)",
            100.0 * faired as f64 / out.v.len().max(1) as f64
        ));
    }
    if mon.cancelled() {
        return Err(Error::Cancelled);
    }

    mon.progress(0.0, "decimating");
    decimate::decimate(&mut out, p.target_faces, &|| mon.cancelled())
        .map_err(|_| Error::Cancelled)?;
    let v: Verification = verify(&out);
    let s = out.size();
    mon.log(&format!(
        "output {} triangles, watertight={}, bodies={}, size={:.2} x {:.2} x {:.2} mm",
        out.f.len(),
        v.watertight,
        v.bodies,
        s[0],
        s[1],
        s[2]
    ));
    let report = Report {
        ok: v.is_volume && v.bodies == 1,
        triangles: out.f.len(),
        watertight: v.watertight,
        is_volume: v.is_volume,
        bodies: v.bodies,
        size_mm: s.map(|x| (x * 1000.0).round() / 1000.0),
        volume_mm3: v.volume_mm3,
        grid: dims,
        voxel_mm: h,
        winding_evaluations: evals,
        components_found: comps,
        components_dropped_voxels: dropped,
        voids_filled_voxels: filled,
        contour_triangles,
        faired_vertices: faired,
        params: Params {
            voxel: h,
            ..p.clone()
        },
    };
    Ok((out, report))
}
