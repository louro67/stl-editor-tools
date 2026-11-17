#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy", "scipy", "scikit-image", "trimesh", "libigl>=2.6",
#   "pymeshlab", "manifold3d", "lxml", "rtree",
# ]
# ///
"""mesh-solidify: turn a broken mesh (holes, open shells, self-intersections)
into one watertight, printable solid.

Pipeline
  1. load and scale to print size (all tolerances are in print mm)
  2. generalized winding number on a voxel grid -> robust inside/outside
     (handles holes and overlapping shells without patching anything)
  3. keep the largest connected solid, drop specks
  4. narrow-band signed distance: sign from the winding number, magnitude from
     the exact distance to the input triangles -> smooth, no voxel terraces
  5. marching cubes + manifold repair
  6. biharmonic fairing of hole patches (vertices far from every input triangle)
     so patches continue the curvature of the surface around them
  7. quadric decimation, then verify watertight / manifold

Run it
  uv run solidify.py in.stl out.3mf --scale 0.1          # deps auto-installed
  mesh-solidify in.stl out.3mf --scale 0.1 --json        # after pip install
  mesh-solidify in.stl --inspect --json                  # diagnose only

Exit codes: 0 ok, 1 error, 2 bad arguments, 3 output failed verification,
            4 grid too large (raise --voxel or pass --max-voxels).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

__version__ = "0.1.1"

DEFAULTS = dict(scale=1.0, voxel=0.125, min_feature=0.0, iso=0.5,
                target_faces=400_000, fair_tau=0.5, fair_grow=2, max_voxels=400_000_000)
EXIT_OK, EXIT_ERR, EXIT_ARGS, EXIT_VERIFY, EXIT_GRID = 0, 1, 2, 3, 4
OUT_FORMATS = (".stl", ".3mf", ".obj", ".ply", ".off", ".glb")


# ---------------------------------------------------------------- arguments
def _num(kind, lo=None, hi=None, lo_open=False, hi_open=False):
    """argparse type: number of `kind` within bounds, with a readable error."""
    def parse(s):
        try:
            x = kind(s)
        except ValueError:
            # let ints be written as 4e8 or 400_000_000.0
            try:
                f = float(s)
            except ValueError:
                f = None
            if kind is not int or f is None or not f.is_integer():
                raise argparse.ArgumentTypeError(f"{s!r} is not a valid {kind.__name__}")
            x = int(f)
        if x != x or x in (float("inf"), float("-inf")):
            raise argparse.ArgumentTypeError(f"{s!r} must be finite")
        if lo is not None and (x <= lo if lo_open else x < lo):
            raise argparse.ArgumentTypeError(f"{s} must be {'>' if lo_open else '>='} {lo}")
        if hi is not None and (x >= hi if hi_open else x > hi):
            raise argparse.ArgumentTypeError(f"{s} must be {'<' if hi_open else '<='} {hi}")
        return x
    return parse


def _scale(s):
    """Scale factor; also accepts a percentage like '10%'."""
    s = s.strip()
    if s.endswith("%"):
        return _num(float, 0, lo_open=True)(s[:-1]) / 100
    return _num(float, 0, lo_open=True)(s)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mesh-solidify",
        description="Turn a broken mesh (holes, open shells, self-intersections) into one "
                    "watertight, printable solid via generalized winding numbers.",
        epilog="Examples:\n"
               "  mesh-solidify model.stl fixed.3mf --scale 10%%\n"
               "  mesh-solidify model.stl fixed.stl --scale 0.1 --voxel 0.1 --json\n"
               "  mesh-solidify model.stl --inspect --json\n\n"
               "Exit codes: 0 ok, 1 error, 2 bad arguments, 3 output failed verification, "
               "4 grid too large.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("src", help="input mesh (STL, OBJ, PLY, 3MF, OFF, GLB ...)")
    p.add_argument("dst", nargs="?", help=f"output mesh; format from extension {', '.join(OUT_FORMATS)}")
    p.add_argument("--scale", type=_scale, default=DEFAULTS["scale"],
                   help="uniform scale applied first, e.g. 0.1 or 10%% (default 1)")
    p.add_argument("--voxel", type=_num(float, 0, lo_open=True), default=DEFAULTS["voxel"],
                   help="grid spacing in output mm; smaller = more detail, more RAM/time (default 0.125)")
    p.add_argument("--min-feature", type=_num(float, 0), default=DEFAULTS["min_feature"],
                   help="uniform outward offset in output mm to fatten thin parts (default 0 = off)")
    p.add_argument("--iso", type=_num(float, 0, 1, lo_open=True, hi_open=True), default=DEFAULTS["iso"],
                   help="winding-number threshold, 0<iso<1; lower closes more (default 0.5)")
    p.add_argument("--target-faces", type=_num(int, 1), default=DEFAULTS["target_faces"],
                   help="decimate to at most this many triangles (default 400000)")
    p.add_argument("--fair-tau", type=_num(float, 0), default=DEFAULTS["fair_tau"],
                   help="vertices farther than this many voxels from the input are hole patch "
                        "and get faired; 0 disables fairing (default 0.5)")
    p.add_argument("--fair-grow", type=_num(int, 0), default=DEFAULTS["fair_grow"],
                   help="extra vertex rings around each patch to refair (default 2)")
    p.add_argument("--max-voxels", type=_num(int, 1), default=DEFAULTS["max_voxels"],
                   help="refuse grids larger than this (default 4e8, about 6 GB RAM peak)")
    p.add_argument("--inspect", action="store_true",
                   help="only report the input's defects (holes, components, size); no dst needed")
    p.add_argument("--json", action="store_true",
                   help="print one JSON object with the result on stdout (progress stays on stderr)")
    p.add_argument("-q", "--quiet", action="store_true", help="no progress output")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


# ---------------------------------------------------------------- helpers
class Log:
    def __init__(self, quiet):
        self.t0, self.quiet = time.time(), quiet

    def __call__(self, *s):
        if not self.quiet:
            print(f"[{time.time() - self.t0:6.1f}s]", *s, file=sys.stderr, flush=True)


def _import_deps():
    try:
        import numpy, scipy, skimage, trimesh, igl, pymeshlab, manifold3d  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            f"mesh-solidify: missing dependency '{e.name}'.\n"
            "Install with:  pip install numpy scipy scikit-image trimesh libigl pymeshlab manifold3d lxml rtree\n"
            "or run with:   uv run solidify.py ...   (installs them automatically)")


def _boundary_loops(mesh):
    """Hole boundary loops as (n_edges, perimeter) via connected components of boundary edges."""
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    e = mesh.edges_sorted
    u, c = np.unique(e, axis=0, return_counts=True)
    b = u[c == 1]
    if not len(b):
        return []
    nodes, inv = np.unique(b, return_inverse=True)
    inv = inv.reshape(-1, 2)
    g = coo_matrix((np.ones(len(inv)), (inv[:, 0], inv[:, 1])), shape=(len(nodes), len(nodes)))
    _, lab = connected_components(g, directed=False)
    ln = np.linalg.norm(mesh.vertices[b[:, 0]] - mesh.vertices[b[:, 1]], axis=1)
    el = lab[inv[:, 0]]
    loops = [(int((el == k).sum()), float(ln[el == k].sum())) for k in range(lab.max() + 1)]
    return sorted(loops, key=lambda x: -x[1])


def inspect(src, scale=1.0):
    import numpy as np
    import trimesh
    m = trimesh.load(src, force="mesh")
    m.apply_scale(scale)
    e = m.edges_sorted
    _, c = np.unique(e, axis=0, return_counts=True)
    loops = _boundary_loops(m)
    parts = m.split(only_watertight=False)
    return {
        "file": os.path.abspath(src), "scale": scale,
        "triangles": int(len(m.faces)), "vertices": int(len(m.vertices)),
        "size_mm": [round(float(x), 3) for x in m.extents],
        "watertight": bool(m.is_watertight), "winding_consistent": bool(m.is_winding_consistent),
        "components": len(parts),
        "boundary_edges": int((c == 1).sum()), "nonmanifold_edges": int((c > 2).sum()),
        "holes": len(loops),
        "largest_hole_perimeters_mm": [round(p, 2) for _, p in loops[:10]],
    }


# ---------------------------------------------------------------- pipeline
def solidify(src, dst, *, scale=1.0, voxel=0.125, min_feature=0.0, iso=0.5, target_faces=400_000,
             fair_tau=0.5, fair_grow=2, max_voxels=400_000_000, log=None):
    """Run the full pipeline; returns a result dict. Raises GridTooLarge / RuntimeError."""
    import numpy as np
    import trimesh, igl, pymeshlab as ml, manifold3d as mf
    import scipy.sparse as sp, scipy.sparse.linalg as spl
    from scipy import ndimage
    from skimage.measure import marching_cubes
    log = log or Log(True)

    m = trimesh.load(src, force="mesh")
    m.apply_scale(scale)
    V = np.ascontiguousarray(m.vertices, float)
    F = np.ascontiguousarray(m.faces, np.int64)
    if not len(F):
        raise RuntimeError("input has no triangles")
    h = voxel
    k_band = 2 + int(np.ceil(min_feature / h))   # band must cover the offset
    pad = k_band + 2
    lo = V.min(0) - pad * h
    dims = np.ceil((V.max(0) + pad * h - lo) / h).astype(int) + 1
    nvox = int(dims.prod())
    if nvox > max_voxels:
        # round up (and allow for padding) so the suggestion actually fits
        suggest = math.ceil(voxel * (nvox / max_voxels) ** (1 / 3) * 1.02 * 1000) / 1000
        raise GridTooLarge(f"grid {dims.tolist()} = {nvox/1e6:.0f}M voxels exceeds --max-voxels "
                           f"{max_voxels/1e6:.0f}M; use --voxel {suggest:.3f} or raise --max-voxels")
    log(f"input {len(F)} tris, watertight={m.is_watertight}; grid {dims.tolist()} = {nvox/1e6:.0f}M voxels")

    # 2. winding number, slab by slab (bounded memory per call)
    xs, ys, zs = [lo[i] + h * np.arange(dims[i]) for i in range(3)]
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    XY = np.c_[X.ravel(), Y.ravel()]
    W = np.empty(dims, np.float32)
    step = max(1, int(4e6 // len(XY)))
    for k in range(0, dims[2], step):
        zz = zs[k:k + step]
        Q = np.vstack([np.c_[XY, np.full(len(XY), z)] for z in zz])
        W[:, :, k:k + len(zz)] = igl.fast_winding_number(V, F, Q).reshape(len(zz), dims[0], dims[1]).transpose(1, 2, 0)
    del X, Y, XY
    log("winding number done")

    # 3. largest component
    B = W > iso
    del W
    lab, n = ndimage.label(B)
    if n == 0:
        raise RuntimeError("nothing is inside the mesh at this iso; try a lower --iso or smaller --voxel")
    sz = np.bincount(lab.ravel()); sz[0] = 0
    B = lab == sz.argmax()
    del lab
    dropped = int(sz.sum() - sz.max())
    log(f"{n} components, kept largest ({int(sz.max())} voxels), dropped {dropped} voxels")

    # 4. narrow-band signed distance (+inside), sign from winding number
    band = ndimage.binary_dilation(B, iterations=k_band) & ~ndimage.binary_erosion(B, iterations=k_band)
    idx = np.nonzero(band)
    del band
    d = np.sqrt(igl.point_mesh_squared_distance(lo + h * np.c_[idx], V, F)[0])
    d = np.maximum(d, 1e-4 * h)                    # keep the sign on points that sit on the surface
    far = (k_band + 1) * h
    field = np.where(B, far, -far).astype(np.float32)
    field[idx] = np.where(B[idx], d, -d)
    del idx, d, B
    if min_feature > 0:
        field += min_feature
    log("signed distance band done")

    # 5. surface + manifold repair
    v, f, _, _ = marching_cubes(field, 0.0, spacing=(h, h, h))
    del field
    v += lo
    ms = ml.MeshSet()
    ms.add_mesh(ml.Mesh(v, f[:, ::-1]))
    ms.meshing_remove_duplicate_vertices()
    ms.meshing_repair_non_manifold_edges()
    ms.meshing_repair_non_manifold_vertices()
    ms.meshing_remove_connected_component_by_face_number(mincomponentsize=1000)
    log(f"marching cubes: {ms.current_mesh().face_number()} tris")

    # 6. biharmonic fairing of hole patches
    n_faired = 0
    if fair_tau > 0:
        cm = ms.current_mesh()
        pv = cm.vertex_matrix().copy()
        pf = cm.face_matrix().astype(np.int64)
        dist = np.sqrt(igl.point_mesh_squared_distance(pv, V, F)[0]) - min_feature
        free = dist > fair_tau * h
        nv = len(pv)
        e = np.r_[pf[:, [0, 1]], pf[:, [1, 2]], pf[:, [2, 0]]]
        A = sp.coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(nv, nv)).tocsr()
        A = ((A + A.T) > 0).astype(float)
        for _ in range(fair_grow):
            free = free | (A @ free.astype(float) > 0)
        n_faired = int(free.sum())
        if 0 < n_faired < nv:
            deg = np.asarray(A.sum(1)).ravel()
            L = sp.eye(nv) - sp.diags(1 / deg) @ A     # uniform Laplacian
            Mb = (L.T @ L).tocsr()                      # bi-Laplacian
            fi, ci = np.nonzero(free)[0], np.nonzero(~free)[0]
            solve = spl.factorized(Mb[fi][:, fi].tocsc())
            rhs = -(Mb[fi][:, ci] @ pv[ci])
            pv[fi] = np.column_stack([solve(rhs[:, j]) for j in range(3)])
            ms = ml.MeshSet()
            ms.add_mesh(ml.Mesh(pv, pf.astype(np.int32)))
        log(f"faired {n_faired} hole-patch vertices ({100 * n_faired / nv:.1f}%)")

    # 7. decimate + verify
    if ms.current_mesh().face_number() > target_faces:
        ms.meshing_decimation_quadric_edge_collapse(targetfacenum=target_faces, preservetopology=True,
                                                    preserveboundary=True, qualitythr=0.3, planarquadric=True)
    ms.meshing_remove_duplicate_vertices()
    ms.meshing_repair_non_manifold_edges()
    cm = ms.current_mesh()
    r = trimesh.Trimesh(cm.vertex_matrix(), cm.face_matrix(), process=True)
    trimesh.repair.fix_normals(r)
    M = mf.Manifold(mf.Mesh(vert_properties=r.vertices.astype(np.float32), tri_verts=r.faces.astype(np.uint32)))
    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    r.export(dst)
    res = {
        "ok": bool(r.is_watertight and r.is_volume and M.status().name == "NoError"),
        "output": os.path.abspath(dst), "triangles": int(len(r.faces)),
        "watertight": bool(r.is_watertight), "is_volume": bool(r.is_volume),
        "bodies": len(r.split()), "manifold3d": M.status().name,
        "size_mm": [round(float(x), 3) for x in r.extents],
        "volume_mm3": round(float(r.volume), 3),
        "components_dropped_voxels": dropped, "faired_vertices": n_faired,
        "params": dict(scale=scale, voxel=voxel, min_feature=min_feature, iso=iso,
                       target_faces=target_faces, fair_tau=fair_tau, fair_grow=fair_grow),
    }
    log(f"out {res['triangles']} tris, watertight={res['watertight']}, bodies={res['bodies']}, "
        f"manifold3d={res['manifold3d']}, size={res['size_mm']} mm -> {dst}")
    return res


class GridTooLarge(RuntimeError):
    pass


# ---------------------------------------------------------------- entry point
def main(argv=None) -> int:
    p = build_parser()
    a = p.parse_args(argv)
    if not os.path.isfile(a.src):
        p.error(f"input not found: {a.src}")
    if not a.inspect:
        if not a.dst:
            p.error("dst is required unless --inspect is given")
        if os.path.splitext(a.dst)[1].lower() not in OUT_FORMATS:
            p.error(f"unsupported output extension; use one of {', '.join(OUT_FORMATS)}")
    _import_deps()
    log = Log(a.quiet)
    t0 = time.time()

    def emit(obj, code):
        obj["seconds"] = round(time.time() - t0, 1)
        obj["version"] = __version__
        if a.json:
            print(json.dumps(obj, indent=2))
        return code

    try:
        if a.inspect:
            info = inspect(a.src, a.scale)
            if not a.json:
                for k, v in info.items():
                    print(f"{k:28s} {v}")
            return emit({"ok": True, "mode": "inspect", **info}, EXIT_OK)
        res = solidify(a.src, a.dst, scale=a.scale, voxel=a.voxel, min_feature=a.min_feature, iso=a.iso,
                       target_faces=a.target_faces, fair_tau=a.fair_tau, fair_grow=a.fair_grow,
                       max_voxels=a.max_voxels, log=log)
        return emit({"mode": "solidify", **res}, EXIT_OK if res["ok"] else EXIT_VERIFY)
    except GridTooLarge as e:
        print(f"mesh-solidify: {e}", file=sys.stderr)
        return emit({"ok": False, "error": str(e), "code": "grid_too_large"}, EXIT_GRID)
    except Exception as e:  # noqa: BLE001 - report every failure in one shape
        print(f"mesh-solidify: {type(e).__name__}: {e}", file=sys.stderr)
        return emit({"ok": False, "error": f"{type(e).__name__}: {e}", "code": "error"}, EXIT_ERR)


if __name__ == "__main__":
    sys.exit(main())
