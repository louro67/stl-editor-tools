use mesh_solidify::{check::verify, inspect, io, solidify, Mesh, NoMonitor, Params};

/// Unit sphere from a subdivided octahedron, outward-oriented.
fn sphere(levels: usize) -> Mesh {
    let mut v: Vec<[f64; 3]> = vec![
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    for _ in 0..levels {
        let mut mid = std::collections::HashMap::new();
        let mut m = |a: u32, b: u32, v: &mut Vec<[f64; 3]>| -> u32 {
            *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let (p, q) = (v[a as usize], v[b as usize]);
                let c = [(p[0] + q[0]) / 2., (p[1] + q[1]) / 2., (p[2] + q[2]) / 2.];
                let l = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
                v.push([c[0] / l, c[1] / l, c[2] / l]);
                (v.len() - 1) as u32
            })
        };
        let mut nf = Vec::new();
        for t in &f {
            let (a, b, c) = (t[0], t[1], t[2]);
            let (ab, bc, ca) = (m(a, b, &mut v), m(b, c, &mut v), m(c, a, &mut v));
            nf.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
        }
        f = nf;
    }
    Mesh { v, f }
}

#[test]
fn closed_sphere_stays_a_sphere() {
    let s = sphere(4);
    assert!(inspect(&s).watertight);
    let p = Params {
        voxel: 0.04,
        target_faces: 3000,
        ..Params::default()
    };
    let (out, r) = solidify(&s, &p, &NoMonitor).unwrap();
    assert!(r.ok && r.watertight && r.bodies == 1, "{r:?}");
    let vol = out.volume();
    let exact = 4.0 / 3.0 * std::f64::consts::PI;
    assert!((vol - exact).abs() < 0.05 * exact, "volume {vol}");
}

#[test]
fn holes_get_closed() {
    let mut s = sphere(4);
    // Cut a big hole (the cap above z = 0.6) and a few small ones.
    let tri = |m: &Mesh, i: usize| m.tri(i);
    let keep: Vec<bool> = (0..s.f.len())
        .map(|i| {
            let c = tri(&s, i).iter().fold([0.0; 3], |a, p| {
                [a[0] + p[0] / 3., a[1] + p[1] / 3., a[2] + p[2] / 3.]
            });
            c[2] < 0.6 && !(c[0] > 0.9 && c[1].abs() < 0.1)
        })
        .collect();
    s.f =
        s.f.iter()
            .zip(&keep)
            .filter(|(_, k)| **k)
            .map(|(t, _)| *t)
            .collect();
    let before = inspect(&s);
    assert!(!before.watertight && before.holes >= 2, "{before:?}");
    let p = Params {
        voxel: 0.04,
        target_faces: 3000,
        ..Params::default()
    };
    let (out, r) = solidify(&s, &p, &NoMonitor).unwrap();
    assert!(r.ok, "{r:?}");
    let v = verify(&out);
    assert!(v.watertight && v.is_volume && v.bodies == 1);
    assert!(out.f.len() <= 3000);
    // The cap is filled, so the volume is near the full sphere's, not a hollow shell's.
    let exact = 4.0 / 3.0 * std::f64::consts::PI;
    assert!(out.volume() > 0.8 * exact, "volume {}", out.volume());
}

#[test]
fn loose_parts_are_dropped_and_scale_applies() {
    let mut s = sphere(3);
    let mut small = sphere(2);
    small.scale(0.2);
    small.translate([3.0, 0.0, 0.0]);
    let off = s.v.len() as u32;
    s.v.extend(small.v);
    s.f.extend(small.f.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    let p = Params {
        scale: 10.0,
        voxel: 0.4,
        target_faces: 2000,
        ..Params::default()
    };
    let (out, r) = solidify(&s, &p, &NoMonitor).unwrap();
    assert!(
        r.ok && r.bodies == 1 && r.components_dropped_voxels > 0,
        "{r:?}"
    );
    let size = out.size();
    assert!((size[0] - 20.0).abs() < 1.0, "{size:?}");
}

#[test]
fn grid_limit_is_reported_with_a_working_suggestion() {
    let s = sphere(2);
    let p = Params {
        voxel: 0.001,
        max_voxels: 1_000_000,
        ..Params::default()
    };
    match solidify(&s, &p, &NoMonitor) {
        Err(mesh_solidify::Error::GridTooLarge { suggest, .. }) => {
            let p2 = Params {
                voxel: suggest,
                max_voxels: 1_000_000,
                target_faces: 1000,
                ..Params::default()
            };
            assert!(
                solidify(&s, &p2, &NoMonitor).is_ok(),
                "suggested voxel {suggest} must fit"
            );
        }
        other => panic!("expected GridTooLarge, got {:?}", other.map(|x| x.1)),
    }
}

#[test]
fn file_formats_round_trip() {
    let s = sphere(2);
    let stl = io::write_stl(&s);
    let back = io::read("x.stl", &stl).unwrap();
    assert_eq!(back.f.len(), s.f.len());
    assert_eq!(back.v.len(), s.v.len());
    let obj = io::write_obj(&s);
    let back = io::read("x.obj", &obj).unwrap();
    assert_eq!(back.f.len(), s.f.len());
    let ply = "ply\nformat ascii 1.0\nelement vertex 4\nproperty float x\nproperty float y\nproperty float z\nelement face 2\nproperty list uchar int vertex_indices\nend_header\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n4 0 1 2 3\n3 0 2 3\n";
    let back = io::read("x.ply", ply.as_bytes()).unwrap();
    assert_eq!(back.f.len(), 3);
}
