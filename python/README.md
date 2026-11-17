# mesh-solidify (Python, legacy)

The original Python implementation (0.1.x), kept for reference. The maintained version is the
Rust command line and SolidCraft plugin in the repository root; it gives the same results,
about twice as fast with a quarter of the memory, and has no Python dependencies.

```bash
uv run solidify.py in.stl out.3mf --scale 10%
pipx install "git+https://github.com/JakeMathews/mesh-solidify#subdirectory=python"
```

It needs numpy, scipy, scikit-image, trimesh, libigl, pymeshlab (GPL-3.0), manifold3d, lxml
and rtree; `uv run` installs them from the script header. Flags match the Rust CLI's.
