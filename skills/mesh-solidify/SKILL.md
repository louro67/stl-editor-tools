---
name: mesh-solidify
description: Use when a mesh (STL, OBJ, 3MF, PLY) is not watertight, has holes, open shells or self-intersections, or a slicer reports floating parts, empty layers or "close holes" failures, and it needs to become one printable solid. Runs the mesh-solidify CLI.
---

# mesh-solidify

Turns a broken mesh into one watertight, manifold solid. It rebuilds the
solid from an inside/outside test instead of patching holes, so it works on
scans, game assets and meshes with huge or non-planar holes, where hole-fill
tools fail.

## Run it

```bash
mesh-solidify IN OUT [options]
```

If it isn't on PATH, download the binary for this platform from
https://github.com/JakeMathews/mesh-solidify/releases/latest, or build it with
`cargo install --git https://github.com/JakeMathews/mesh-solidify mesh-solidify-cli`.
Inputs: .stl, .obj, .ply. Outputs: .stl, .3mf, .obj.

Always pass `--json`. Progress goes to stderr and one JSON object goes to
stdout, so parse stdout only.

## Workflow

1. Inspect first. It's cheap (seconds) and needs no output path:
   `mesh-solidify in.stl --inspect --scale 10% --json`.
   Read `holes`, `components`, `watertight` and `size_mm`. If `watertight` is already true and
   `components` is 1, the mesh probably doesn't need this tool.
2. Pick the print scale. All tolerances are in output mm, so always pass `--scale`
   for the final print size (`0.1` or `10%`).
3. Pick `--voxel`. It's the detail size in output mm. Default 0.125 suits a 0.4 mm nozzle.
   The grid grows with (size / voxel)^3 at 1 byte per voxel; 130M voxels take ~35 s and ~3 GB.
   If you get exit code 4, use the `--voxel` value suggested in the error, or `--voxel auto`.
4. Run: `mesh-solidify in.stl out.3mf --scale 10% --json`.
5. Check the JSON: `ok` true, `watertight` true, `bodies` 1, `size_mm` as expected.

## Options worth changing

| option | default | when to change |
|---|---|---|
| `--scale` | 1 | always set it to the print scale |
| `--voxel` | 0.125 | smaller for finer detail, larger for speed or memory |
| `--min-feature` | 0 | thin parts vanish in the slicer: offset the whole surface outward by this many mm (e.g. 0.2) |
| `--iso` | 0.5 | a cavity behind a big opening should be filled: try 0.4 (lower fills more) |
| `--target-faces` | 400000 | output triangle budget |
| `--fair-tau` | 0.5 | 0 disables smoothing of hole patches |

## Exit codes

0 ok · 1 error · 2 bad arguments · 3 output failed verification · 4 grid too large

## In SolidCraft

The same tool is a SolidCraft plugin (Plugins › Mesh Solidify), but SolidCraft 0.41 only lets a
plugin read bodies of up to 2,730 triangles. For larger meshes run the command line and insert
the result with `file.insert_mesh`.

## Limits

- The output is solid. Hollow parts (an empty tank) come out filled.
- It keeps only the largest connected solid. Separate loose parts are dropped
  (the JSON field `components_dropped_voxels` reports how much).
- Features thinner than about one voxel can disappear. Use a smaller `--voxel` or `--min-feature`.
