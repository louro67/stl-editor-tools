//! mesh-solidify: turn a broken triangle mesh (holes, open shells, overlapping or
//! self-intersecting parts) into one watertight, manifold solid.
//!
//! It does not patch holes. It rebuilds the solid from an inside/outside test:
//!
//! 1. **Generalized winding number** of the input, evaluated adaptively on a grid
//!    (full resolution only near the surface and where the field is ambiguous).
//!    Holes and overlaps don't matter to it.
//! 2. **Largest solid** kept; enclosed voids filled.
//! 3. **Narrow-band signed distance** to the input, so the rebuilt surface follows the
//!    original faces instead of voxel steps.
//! 4. **Marching tetrahedra**, which is crack-free and manifold by construction.
//! 5. **Biharmonic fairing** of the vertices that patch holes, so patches follow the
//!    curvature of the surrounding surface.
//! 6. **Quadric decimation** to a triangle budget, keeping the mesh closed and manifold.
//!
//! The crate has no I/O or platform dependencies beyond mesh file parsing, so it builds for
//! WebAssembly (the SolidCraft plugin) as well as natively (the `mesh-solidify` CLI, which
//! enables the `parallel` feature).

#![allow(clippy::needless_range_loop)]

pub mod bvh;
pub mod check;
mod contour;
mod decimate;
mod fair;
mod grid;
pub mod io;
pub mod mesh;
mod par;
mod pipeline;

pub use check::{inspect, Inspection, Verification};
pub use mesh::Mesh;
pub use pipeline::{auto_voxel, grid_dims, solidify, Error, Monitor, NoMonitor, Params, Report};

/// Simplify a closed manifold mesh to about `target` triangles, keeping it closed and manifold.
pub fn simplify(m: &mut Mesh, target: usize) {
    let _ = decimate::decimate(m, target, &|| false);
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
