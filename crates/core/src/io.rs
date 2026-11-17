//! Mesh file formats without dependencies: STL (binary and ASCII), OBJ and PLY in;
//! binary STL and OBJ out. (3MF output lives in the CLI, which has a zip writer.)

use crate::mesh::{Mesh, V3};

/// Read a mesh, choosing the format from `name`'s extension, falling back to sniffing.
pub fn read(name: &str, data: &[u8]) -> Result<Mesh, String> {
    let ext = name.rsplit('.').next().unwrap_or("").to_ascii_lowercase();
    let mesh = match ext.as_str() {
        "stl" => read_stl(data),
        "obj" => read_obj(data),
        "ply" => read_ply(data),
        _ if data.starts_with(b"ply") => read_ply(data),
        _ => read_stl(data),
    }?;
    if mesh.f.is_empty() {
        return Err("the file has no triangles".into());
    }
    if mesh.v.iter().flatten().any(|x| !x.is_finite()) {
        return Err("the file has non-finite coordinates".into());
    }
    Ok(mesh)
}

pub fn read_stl(data: &[u8]) -> Result<Mesh, String> {
    if data.len() >= 84 {
        let n = u32::from_le_bytes(data[80..84].try_into().unwrap()) as usize;
        if data.len() == 84 + 50 * n || (!looks_ascii(data) && data.len() >= 84 + 50 * n) {
            let mut tris = Vec::with_capacity(n);
            for i in 0..n {
                let rec = &data[84 + 50 * i..84 + 50 * i + 50];
                let mut t = [[0.0; 3]; 3];
                for (k, p) in t.iter_mut().enumerate() {
                    for (j, c) in p.iter_mut().enumerate() {
                        let o = 12 + 12 * k + 4 * j;
                        *c = f32::from_le_bytes(rec[o..o + 4].try_into().unwrap()) as f64;
                    }
                }
                tris.push(t);
            }
            return Ok(Mesh::from_soup(tris));
        }
    }
    let text = std::str::from_utf8(data).map_err(|_| "not a valid STL file")?;
    let mut tris = Vec::new();
    let mut cur: Vec<V3> = Vec::with_capacity(3);
    for line in text.lines() {
        let mut it = line.split_whitespace();
        match it.next() {
            Some("vertex") => {
                let p = parse3(&mut it).ok_or("bad ASCII STL vertex")?;
                cur.push(p);
            }
            Some("endloop") => {
                // Fan-triangulate loops with more than three vertices.
                for k in 1..cur.len().saturating_sub(1) {
                    tris.push([cur[0], cur[k], cur[k + 1]]);
                }
                cur.clear();
            }
            _ => {}
        }
    }
    if tris.is_empty() {
        return Err("not a valid STL file".into());
    }
    Ok(Mesh::from_soup(tris))
}

fn looks_ascii(data: &[u8]) -> bool {
    let head = &data[..data.len().min(512)];
    let t = head
        .iter()
        .skip_while(|b| b.is_ascii_whitespace())
        .take(5)
        .copied()
        .collect::<Vec<_>>();
    t.eq_ignore_ascii_case(b"solid") && head.iter().all(|b| b.is_ascii())
}

fn parse3<'a>(it: &mut impl Iterator<Item = &'a str>) -> Option<V3> {
    Some([
        it.next()?.parse().ok()?,
        it.next()?.parse().ok()?,
        it.next()?.parse().ok()?,
    ])
}

pub fn read_obj(data: &[u8]) -> Result<Mesh, String> {
    let text = String::from_utf8_lossy(data);
    let mut m = Mesh::default();
    for line in text.lines() {
        let mut it = line.split_whitespace();
        match it.next() {
            Some("v") => m.v.push(parse3(&mut it).ok_or("bad OBJ vertex")?),
            Some("f") => {
                let n = m.v.len() as i64;
                let idx: Vec<u32> = it
                    .map(|tok| {
                        let i: i64 = tok
                            .split('/')
                            .next()
                            .unwrap_or("")
                            .parse()
                            .map_err(|_| "bad OBJ face")?;
                        let i = if i < 0 { n + i } else { i - 1 };
                        if i < 0 || i >= n {
                            return Err("OBJ face index out of range");
                        }
                        Ok(i as u32)
                    })
                    .collect::<Result<_, _>>()?;
                for k in 1..idx.len().saturating_sub(1) {
                    m.push_face([idx[0], idx[k], idx[k + 1]]);
                }
            }
            _ => {}
        }
    }
    weld(m)
}

/// Merge vertices with identical coordinates (OBJ and PLY may duplicate them per face).
fn weld(m: Mesh) -> Result<Mesh, String> {
    Ok(Mesh::from_soup((0..m.f.len()).map(|i| m.tri(i))))
}

#[derive(Clone, Copy, PartialEq)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    F32,
    F64,
}

impl Ty {
    fn parse(s: &str) -> Option<Ty> {
        Some(match s {
            "char" | "int8" => Ty::I8,
            "uchar" | "uint8" => Ty::U8,
            "short" | "int16" => Ty::I16,
            "ushort" | "uint16" => Ty::U16,
            "int" | "int32" => Ty::I32,
            "uint" | "uint32" => Ty::U32,
            "float" | "float32" => Ty::F32,
            "double" | "float64" => Ty::F64,
            _ => return None,
        })
    }
    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::F64 => 8,
        }
    }
}

enum Prop {
    Scalar(String, Ty),
    List(String, Ty, Ty),
}

struct Element {
    name: String,
    count: usize,
    props: Vec<Prop>,
}

pub fn read_ply(data: &[u8]) -> Result<Mesh, String> {
    let end = find(data, b"end_header").ok_or("PLY header has no end_header")?;
    let body_start = data[end..]
        .iter()
        .position(|&b| b == b'\n')
        .map(|p| end + p + 1)
        .ok_or("bad PLY header")?;
    let header = std::str::from_utf8(&data[..end]).map_err(|_| "bad PLY header")?;
    let mut format = "";
    let mut elements: Vec<Element> = Vec::new();
    for line in header.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.as_slice() {
            ["format", f, ..] => format = f,
            ["element", name, n] => elements.push(Element {
                name: name.to_string(),
                count: n.parse().map_err(|_| "bad PLY element count")?,
                props: Vec::new(),
            }),
            ["property", "list", c, t, name] => elements
                .last_mut()
                .ok_or("PLY property before element")?
                .props
                .push(Prop::List(
                    name.to_string(),
                    Ty::parse(c).ok_or("bad PLY type")?,
                    Ty::parse(t).ok_or("bad PLY type")?,
                )),
            ["property", t, name] => elements
                .last_mut()
                .ok_or("PLY property before element")?
                .props
                .push(Prop::Scalar(
                    name.to_string(),
                    Ty::parse(t).ok_or("bad PLY type")?,
                )),
            _ => {}
        }
    }
    let mut r = PlyReader::new(&data[body_start..], format)?;
    let mut m = Mesh::default();
    for el in &elements {
        for _ in 0..el.count {
            let mut p = [0.0; 3];
            for prop in &el.props {
                match prop {
                    Prop::Scalar(name, t) => {
                        let x = r.num(*t)?;
                        if el.name == "vertex" {
                            match name.as_str() {
                                "x" => p[0] = x,
                                "y" => p[1] = x,
                                "z" => p[2] = x,
                                _ => {}
                            }
                        }
                    }
                    Prop::List(name, ct, it) => {
                        let n = r.num(*ct)? as usize;
                        let mut idx = Vec::with_capacity(n);
                        for _ in 0..n {
                            idx.push(r.num(*it)? as i64);
                        }
                        if el.name == "face" && (name == "vertex_indices" || name == "vertex_index")
                        {
                            for k in 1..idx.len().saturating_sub(1) {
                                m.push_face([idx[0] as u32, idx[k] as u32, idx[k + 1] as u32]);
                            }
                        }
                    }
                }
            }
            if el.name == "vertex" {
                m.v.push(p);
            }
        }
    }
    let nv = m.v.len() as u32;
    if m.f.iter().flatten().any(|&i| i >= nv) {
        return Err("PLY face index out of range".into());
    }
    weld(m)
}

fn find(h: &[u8], n: &[u8]) -> Option<usize> {
    h.windows(n.len()).position(|w| w == n)
}

struct PlyReader<'a> {
    data: &'a [u8],
    pos: usize,
    mode: u8, // 0 ascii, 1 little endian, 2 big endian
    tokens: Option<std::str::SplitAsciiWhitespace<'a>>,
}

impl<'a> PlyReader<'a> {
    fn new(data: &'a [u8], format: &str) -> Result<Self, String> {
        let mode = match format {
            "ascii" => 0,
            "binary_little_endian" => 1,
            "binary_big_endian" => 2,
            _ => return Err(format!("unsupported PLY format {format:?}")),
        };
        let tokens = if mode == 0 {
            Some(
                std::str::from_utf8(data)
                    .map_err(|_| "bad ASCII PLY")?
                    .split_ascii_whitespace(),
            )
        } else {
            None
        };
        Ok(PlyReader {
            data,
            pos: 0,
            mode,
            tokens,
        })
    }

    fn num(&mut self, t: Ty) -> Result<f64, String> {
        if let Some(tok) = &mut self.tokens {
            return tok
                .next()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| "truncated PLY".into());
        }
        let n = t.size();
        let b = self
            .data
            .get(self.pos..self.pos + n)
            .ok_or("truncated PLY")?;
        self.pos += n;
        let mut a = [0u8; 8];
        a[..n].copy_from_slice(b);
        if self.mode == 2 {
            a[..n].reverse();
        }
        Ok(match t {
            Ty::I8 => a[0] as i8 as f64,
            Ty::U8 => a[0] as f64,
            Ty::I16 => i16::from_le_bytes([a[0], a[1]]) as f64,
            Ty::U16 => u16::from_le_bytes([a[0], a[1]]) as f64,
            Ty::I32 => i32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
            Ty::U32 => u32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
            Ty::F32 => f32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
            Ty::F64 => f64::from_le_bytes(a),
        })
    }
}

/// Binary STL.
pub fn write_stl(m: &Mesh) -> Vec<u8> {
    let mut out = Vec::with_capacity(84 + 50 * m.f.len());
    let mut header = [0u8; 80];
    let tag = b"mesh-solidify";
    header[..tag.len()].copy_from_slice(tag);
    out.extend_from_slice(&header);
    out.extend_from_slice(&(m.f.len() as u32).to_le_bytes());
    for i in 0..m.f.len() {
        let [a, b, c] = m.tri(i);
        let n = crate::mesh::cross(crate::mesh::sub(b, a), crate::mesh::sub(c, a));
        let l = crate::mesh::norm(n).max(1e-30);
        for x in [n[0] / l, n[1] / l, n[2] / l] {
            out.extend_from_slice(&(x as f32).to_le_bytes());
        }
        for p in [a, b, c] {
            for x in p {
                out.extend_from_slice(&(x as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&[0, 0]);
    }
    out
}

pub fn write_obj(m: &Mesh) -> Vec<u8> {
    use std::fmt::Write;
    let mut s = String::with_capacity(m.v.len() * 40 + m.f.len() * 24);
    s.push_str("# mesh-solidify\n");
    for p in &m.v {
        let _ = writeln!(s, "v {} {} {}", p[0] as f32, p[1] as f32, p[2] as f32);
    }
    for t in &m.f {
        let _ = writeln!(s, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1);
    }
    s.into_bytes()
}
