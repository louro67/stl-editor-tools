//! mesh-solidify command line. Progress goes to stderr; with --json, stdout gets exactly one
//! JSON object. Exit codes: 0 ok, 1 error, 2 bad arguments, 3 output failed verification,
//! 4 grid too large.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;
use std::time::Instant;

use clap::Parser;
use mesh_solidify::{inspect, io, solidify, Error, Mesh, Monitor, Params};
use serde_json::{json, Value};

const AFTER: &str = "\
Examples:
  mesh-solidify model.stl fixed.3mf --scale 10%
  mesh-solidify scan.ply fixed.stl --voxel 0.2 --min-feature 0.2
  mesh-solidify model.stl --inspect --json

Exit codes: 0 ok, 1 error, 2 bad arguments, 3 output failed verification, 4 grid too large.";

#[derive(Parser, Debug)]
#[command(name = "mesh-solidify", version, about = "Turn a broken mesh (holes, open shells, self-intersections) into one watertight, printable solid.", after_help = AFTER)]
struct Args {
    /// Input mesh (.stl, .obj, .ply)
    src: PathBuf,
    /// Output mesh (.stl, .3mf, .obj); not needed with --inspect
    dst: Option<PathBuf>,
    /// Uniform scale, as a factor or percentage: 0.1 or 10%
    #[arg(long, default_value = "1", value_parser = scale)]
    scale: f64,
    /// Grid spacing in output mm (detail size), or "auto" to fit --max-voxels
    #[arg(long, default_value = "0.125", value_parser = voxel)]
    voxel: f64,
    /// Offset the surface outward by this many mm so thin parts survive slicing
    #[arg(long, default_value = "0", value_parser = non_negative)]
    min_feature: f64,
    /// Winding-number threshold for inside, between 0 and 1; lower fills more
    #[arg(long, default_value = "0.5", value_parser = unit)]
    iso: f64,
    /// Output triangle budget
    #[arg(long, default_value = "400000", value_parser = count)]
    target_faces: usize,
    /// Fair vertices farther than this many voxels from the input (hole patches); 0 = off
    #[arg(long, default_value = "0.5", value_parser = non_negative)]
    fair_tau: f64,
    /// Extra vertex rings faired around each patch
    #[arg(long, default_value = "2", value_parser = count)]
    fair_grow: usize,
    /// Refuse grids with more nodes than this (about 1 byte of RAM each, plus the surface)
    #[arg(long, default_value = "400000000", value_parser = count)]
    max_voxels: usize,
    /// Worker threads (default: all cores)
    #[arg(long, value_parser = count)]
    threads: Option<usize>,
    /// Only report the input's defects (holes, components, size)
    #[arg(long)]
    inspect: bool,
    /// Print one JSON object on stdout
    #[arg(long)]
    json: bool,
    /// No progress on stderr
    #[arg(short, long)]
    quiet: bool,
}

fn number(s: &str) -> Result<f64, String> {
    let x: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("{s:?} is not a number"))?;
    if !x.is_finite() {
        return Err(format!("{s} must be finite"));
    }
    Ok(x)
}
fn scale(s: &str) -> Result<f64, String> {
    let s = s.trim();
    let x = match s.strip_suffix('%') {
        Some(p) => number(p)? / 100.0,
        None => number(s)?,
    };
    if x <= 0.0 {
        return Err(format!("{s} must be greater than 0"));
    }
    Ok(x)
}
fn voxel(s: &str) -> Result<f64, String> {
    if s.eq_ignore_ascii_case("auto") {
        return Ok(0.0);
    }
    let x = number(s)?;
    if x <= 0.0 {
        return Err(format!("{s} must be greater than 0 (or \"auto\")"));
    }
    Ok(x)
}
fn non_negative(s: &str) -> Result<f64, String> {
    let x = number(s)?;
    if x < 0.0 {
        return Err(format!("{s} must be 0 or more"));
    }
    Ok(x)
}
fn unit(s: &str) -> Result<f64, String> {
    let x = number(s)?;
    if !(x > 0.0 && x < 1.0) {
        return Err(format!("{s} must be between 0 and 1"));
    }
    Ok(x)
}
/// Whole numbers, also written as 4e8 or 400_000_000.
fn count(s: &str) -> Result<usize, String> {
    let x = number(&s.replace('_', ""))?;
    if x < 0.0 || x.fract() != 0.0 || x > 1e15 {
        return Err(format!("{s} must be a whole number of 0 or more"));
    }
    Ok(x as usize)
}

struct Log {
    start: Instant,
    quiet: bool,
    stage: Mutex<String>,
}

impl Log {
    fn line(&self, msg: &str) {
        if !self.quiet {
            eprintln!("[{:7.1}s] {msg}", self.start.elapsed().as_secs_f64());
        }
    }
}

impl Monitor for Log {
    fn progress(&self, _fraction: f32, stage: &str) {
        let mut last = self.stage.lock().unwrap();
        if *last != stage {
            *last = stage.to_string();
            drop(last);
            self.line(&format!("{stage}…"));
        }
    }
    fn log(&self, message: &str) {
        self.line(message);
    }
}

fn write_3mf(m: &Mesh, path: &Path) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let mut model = String::with_capacity(m.v.len() * 60 + m.f.len() * 50);
    model.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<model unit=\"millimeter\" xml:lang=\"en-US\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n<metadata name=\"Application\">mesh-solidify</metadata>\n<resources><object id=\"1\" type=\"model\"><mesh><vertices>\n");
    for p in &m.v {
        let _ = writeln!(
            model,
            "<vertex x=\"{}\" y=\"{}\" z=\"{}\"/>",
            p[0] as f32, p[1] as f32, p[2] as f32
        );
    }
    model.push_str("</vertices><triangles>\n");
    for t in &m.f {
        let _ = writeln!(
            model,
            "<triangle v1=\"{}\" v2=\"{}\" v3=\"{}\"/>",
            t[0], t[1], t[2]
        );
    }
    model.push_str(
        "</triangles></mesh></object></resources><build><item objectid=\"1\"/></build></model>\n",
    );
    let file = std::fs::File::create(path)?;
    let mut z = zip::ZipWriter::new(std::io::BufWriter::new(file));
    let opt = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);
    z.start_file("[Content_Types].xml", opt)?;
    z.write_all(b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/></Types>\n")?;
    z.start_file("_rels/.rels", opt)?;
    z.write_all(b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/></Relationships>\n")?;
    z.start_file("3D/3dmodel.model", opt)?;
    z.write_all(model.as_bytes())?;
    z.finish()?;
    Ok(())
}

fn write_mesh(m: &Mesh, path: &Path) -> Result<(), String> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let r = match ext.as_str() {
        "stl" => std::fs::write(path, io::write_stl(m)),
        "obj" => std::fs::write(path, io::write_obj(m)),
        "3mf" => write_3mf(m, path),
        _ => {
            return Err(format!(
                "unsupported output format {ext:?}; use .stl, .3mf or .obj"
            ))
        }
    };
    r.map_err(|e| format!("cannot write {}: {e}", path.display()))
}

fn emit(json_mode: bool, v: &Value) {
    if json_mode {
        println!("{}", serde_json::to_string_pretty(v).unwrap());
    }
}

fn fail(json_mode: bool, code: u8, kind: &str, msg: &str, start: &Instant) -> ExitCode {
    if json_mode {
        emit(
            true,
            &json!({"ok": false, "error": msg, "code": kind, "seconds": round1(start), "version": mesh_solidify::VERSION}),
        );
    } else {
        eprintln!("mesh-solidify: error: {msg}");
    }
    ExitCode::from(code)
}

fn round1(start: &Instant) -> f64 {
    (start.elapsed().as_secs_f64() * 10.0).round() / 10.0
}

fn main() -> ExitCode {
    let a = Args::parse();
    let start = Instant::now();
    let log = Log {
        start,
        quiet: a.quiet,
        stage: Mutex::new(String::new()),
    };
    if let Some(t) = a.threads {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(t.max(1))
            .build_global();
    }
    if !a.inspect && a.dst.is_none() {
        eprintln!("mesh-solidify: error: an output path is needed (or pass --inspect)");
        return ExitCode::from(2);
    }
    let name = a.src.to_string_lossy().to_string();
    let data = match std::fs::read(&a.src) {
        Ok(d) => d,
        Err(e) => {
            return fail(
                a.json,
                1,
                "read_failed",
                &format!("cannot read {name}: {e}"),
                &start,
            )
        }
    };
    let mut mesh = match io::read(&name, &data) {
        Ok(m) => m,
        Err(e) => return fail(a.json, 1, "read_failed", &format!("{name}: {e}"), &start),
    };
    drop(data);

    if a.inspect {
        mesh.scale(a.scale);
        let r = inspect(&mesh);
        let mut v = serde_json::to_value(&r).unwrap();
        let o = v.as_object_mut().unwrap();
        o.insert("ok".into(), json!(true));
        o.insert("mode".into(), json!("inspect"));
        o.insert("file".into(), json!(name));
        o.insert("scale".into(), json!(a.scale));
        o.insert("seconds".into(), json!(round1(&start)));
        o.insert("version".into(), json!(mesh_solidify::VERSION));
        if a.json {
            emit(true, &v);
        } else {
            println!(
                "{name}: {} triangles, {} components, {} holes, watertight={}, size {:.2} x {:.2} x {:.2} mm",
                r.triangles, r.components, r.holes, r.watertight, r.size_mm[0], r.size_mm[1], r.size_mm[2]
            );
        }
        return ExitCode::SUCCESS;
    }

    let dst = a.dst.clone().unwrap();
    let p = Params {
        scale: a.scale,
        voxel: a.voxel,
        min_feature: a.min_feature,
        iso: a.iso,
        target_faces: a.target_faces,
        fair_tau: a.fair_tau,
        fair_grow: a.fair_grow,
        max_voxels: a.max_voxels,
    };
    log.line(&format!("read {name}: {} triangles", mesh.f.len()));
    let (out, report) = match solidify(&mesh, &p, &log) {
        Ok(x) => x,
        Err(e @ Error::GridTooLarge { .. }) => {
            return fail(a.json, 4, "grid_too_large", &e.to_string(), &start)
        }
        Err(e @ Error::InvalidParam(_)) => {
            return fail(a.json, 2, "invalid_argument", &e.to_string(), &start)
        }
        Err(e) => return fail(a.json, 1, "failed", &e.to_string(), &start),
    };
    if let Err(e) = write_mesh(&out, &dst) {
        return fail(a.json, 1, "write_failed", &e, &start);
    }
    log.line(&format!("wrote {}", dst.display()));
    let mut v = serde_json::to_value(&report).unwrap();
    let o = v.as_object_mut().unwrap();
    o.insert("mode".into(), json!("solidify"));
    o.insert(
        "output".into(),
        json!(std::fs::canonicalize(&dst).unwrap_or(dst).to_string_lossy()),
    );
    o.insert("seconds".into(), json!(round1(&start)));
    o.insert("version".into(), json!(mesh_solidify::VERSION));
    emit(a.json, &v);
    if report.ok {
        ExitCode::SUCCESS
    } else {
        if !a.json {
            eprintln!(
                "mesh-solidify: the output failed verification (watertight={}, bodies={})",
                report.watertight, report.bodies
            );
        }
        ExitCode::from(3)
    }
}
