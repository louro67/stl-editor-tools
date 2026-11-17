# Changelog

## 0.2.0

Rewritten in Rust, primarily as a SolidCraft plugin.

- **SolidCraft plugin** (`mesh-solidify-plugin.zip`, SolidCraft 0.41+): *Solidify Mesh…* rebuilds
  the selected mesh body as a new watertight body (one undo step); *Inspect Mesh…* reports its
  defects. Its dialog defaults to 10% scale. Within SolidCraft 0.41's plugin limits: input up to
  2,730 triangles, output up to about 30,000.
- **Native command line** with the same flags, JSON output and exit codes as 0.1. No Python
  dependencies. About 2× faster than 0.1 and a quarter of the memory on large grids (978k-triangle
  compressor: 35 s instead of 76 s; 650k-triangle toolbox at 400M voxels: 34 s and 2 GB instead
  of 123 s and 7.8 GB).
- Octree winding-number evaluation: only a thin shell around the surface is evaluated at full
  resolution.
- Marching tetrahedra instead of marching cubes plus repair; decimation keeps the mesh manifold.
- Enclosed voids are filled.
- `--voxel auto` picks the finest spacing that fits `--max-voxels`; the exit-code-4 suggestion is
  now exact.
- New `--threads`. Reads PLY. Writes `.stl`, `.3mf`, `.obj` (0.1 also wrote `.ply`, `.off`, `.glb`).
- The Python implementation moved to `python/` (0.1.1: accepts `4e8`-style counts and rounds the
  voxel suggestion up).

## 0.1.0

First release: Python script and package.
