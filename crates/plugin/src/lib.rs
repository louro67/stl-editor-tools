//! mesh-solidify as a SolidCraft plugin (`craft:plugin@0.1.0`, capability `commands@0.1`).
//!
//! Entries:
//! - `solidify` (edit): reads a mesh body with `model.mesh`, rebuilds it as one watertight
//!   solid and adds the result as a new mesh body with `timeline.add`. One undo step.
//! - `inspect` (read): reports holes, loose parts and size of a mesh body.
//!
//! SolidCraft 0.41 limits what a plugin can exchange: `model.mesh` returns three vertices per
//! triangle and refuses more than 8,192, so a mesh body can have at most 2,730 triangles; and
//! every command's JSON is capped at 1 MiB, so the output is simplified until it fits. Larger
//! scans need the `mesh-solidify` command line for now.

use mesh_solidify::{inspect, simplify, solidify, Error, Mesh, Monitor, Params};
use serde::Deserialize;
use serde_json::{json, Value};

wit_bindgen::generate!({
    path: "wit",
    world: "command-plugin",
});

use craft::plugin::types::{CommandFailure, FieldFailure, PluginError};
use craft::plugin::{commands, control, log};

struct Plugin;
export!(Plugin);

/// Output JSON budget: the host caps every command's arguments at 1 MiB.
const MAX_ARGS: usize = 1_000_000;

#[derive(Deserialize)]
#[serde(default)]
struct SolidifyInput {
    body: String,
    scale: f64,
    voxel: f64,
    min_feature: f64,
    iso: f64,
    target_faces: usize,
    fair_tau: f64,
    fair_grow: usize,
    max_voxels: usize,
}

impl Default for SolidifyInput {
    fn default() -> Self {
        SolidifyInput {
            body: String::new(),
            scale: 0.1,
            voxel: 0.0,
            min_feature: 0.0,
            iso: 0.5,
            target_faces: 20_000,
            fair_tau: 0.5,
            fair_grow: 2,
            max_voxels: 1_500_000,
        }
    }
}

#[derive(Deserialize)]
#[serde(default)]
struct InspectInput {
    body: String,
    scale: f64,
}

impl Default for InspectInput {
    fn default() -> Self {
        InspectInput {
            body: String::new(),
            scale: 1.0,
        }
    }
}

struct Host;
impl Monitor for Host {
    fn progress(&self, fraction: f32, stage: &str) {
        control::progress(fraction.clamp(0.0, 1.0), stage);
    }
    fn log(&self, message: &str) {
        log::log(log::Level::Info, message);
    }
    fn cancelled(&self) -> bool {
        control::cancelled()
    }
}

fn invalid(field: &str, reason: impl Into<String>) -> PluginError {
    PluginError::InvalidInput(FieldFailure {
        field: field.into(),
        reason: reason.into(),
    })
}

fn failed(reason: impl Into<String>) -> PluginError {
    PluginError::CommandFailed(CommandFailure {
        command: "mesh-solidify".into(),
        reason: reason.into(),
    })
}

fn run(command: &str, args: &Value) -> Result<Value, PluginError> {
    let out = commands::run(command, &args.to_string())?;
    serde_json::from_str(&out).map_err(|e| failed(format!("{command} returned invalid JSON: {e}")))
}

/// The body named in the input, else the selected body (or the body of a selected face).
fn pick_body(name: &str) -> Result<String, PluginError> {
    if !name.trim().is_empty() {
        return Ok(name.trim().to_string());
    }
    let sel = run("selection.get", &json!({}))?;
    for item in sel["items"].as_array().into_iter().flatten() {
        if item["kind"] == "body" || item["type"] == "body" {
            if let Some(n) = item["name"].as_str() {
                return Ok(n.to_string());
            }
        }
        if let Some(n) = item["body"].as_str() {
            return Ok(n.to_string());
        }
    }
    Err(invalid(
        "body",
        "select a mesh body first, or type its name",
    ))
}

fn read_body(body: &str) -> Result<Mesh, PluginError> {
    let r = run("model.mesh", &json!({"body": body})).map_err(|e| match e {
        PluginError::CommandFailed(f) if f.reason.contains("budget") => failed(format!(
            "{body} has more than 2,730 triangles, the most SolidCraft 0.41 lets a plugin read. \
             Use the mesh-solidify command line for this mesh, then insert its result."
        )),
        other => other,
    })?;
    let pos = r["positions"]
        .as_array()
        .ok_or_else(|| failed("model.mesh returned no positions"))?;
    let tri = r["triangles"]
        .as_array()
        .ok_or_else(|| failed("model.mesh returned no triangles"))?;
    let p = |v: &Value| -> [f64; 3] {
        let a = v.as_array();
        let c = |i: usize| {
            a.and_then(|a| a.get(i))
                .and_then(Value::as_f64)
                .unwrap_or(0.0)
        };
        [c(0), c(1), c(2)]
    };
    let verts: Vec<[f64; 3]> = pos.iter().map(p).collect();
    let n = verts.len();
    let mut soup = Vec::with_capacity(tri.len());
    for t in tri {
        let a = t.as_array().ok_or_else(|| failed("bad triangle"))?;
        let i = |k: usize| {
            a.get(k)
                .and_then(Value::as_u64)
                .map(|x| x as usize)
                .filter(|&x| x < n)
        };
        let (Some(x), Some(y), Some(z)) = (i(0), i(1), i(2)) else {
            return Err(failed("triangle index out of range"));
        };
        soup.push([verts[x], verts[y], verts[z]]);
    }
    let m = Mesh::from_soup(soup);
    if m.f.is_empty() {
        return Err(failed(format!("{body} has no triangles")));
    }
    Ok(m)
}

/// `timeline.add` arguments for a mesh body, with coordinates rounded to 0.1 µm.
fn mesh_feature(m: &Mesh, name: &str) -> String {
    use std::fmt::Write;
    let mut s = String::with_capacity(m.v.len() * 30 + m.f.len() * 20 + 200);
    let fname = serde_json::to_string(name).unwrap();
    let _ = write!(
        s,
        r#"{{"feature":{{"type":"mesh_import","name":{fname},"file":{fname},"meshes":[{{"name":{fname},"positions":["#
    );
    for (i, p) in m.v.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        for (k, x) in p.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            let r = (x * 10_000.0).round() / 10_000.0;
            let _ = write!(s, "{r}");
        }
    }
    s.push_str(r#"],"triangles":["#);
    for (i, t) in m.f.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(s, "{},{},{}", t[0], t[1], t[2]);
    }
    s.push_str("]}]}}");
    s
}

fn do_solidify(params: &str) -> Result<String, PluginError> {
    let input: SolidifyInput =
        serde_json::from_str(params).map_err(|e| invalid("params", e.to_string()))?;
    let body = pick_body(&input.body)?;
    let mesh = read_body(&body)?;
    let p = Params {
        scale: input.scale,
        voxel: input.voxel,
        min_feature: input.min_feature,
        iso: input.iso,
        target_faces: input.target_faces,
        fair_tau: input.fair_tau,
        fair_grow: input.fair_grow,
        max_voxels: input.max_voxels,
    };
    let (mut out, report) = solidify(&mesh, &p, &Host).map_err(|e| match e {
        Error::Cancelled => PluginError::Cancelled,
        Error::InvalidParam(s) => invalid("params", s),
        e @ Error::GridTooLarge { .. } => invalid("voxel", e.to_string()),
        e => failed(e.to_string()),
    })?;
    let name = format!("{body} solid");
    let mut args = mesh_feature(&out, &name);
    while args.len() > MAX_ARGS {
        let target = (out.f.len() as f64 * 0.8 * MAX_ARGS as f64 / args.len() as f64) as usize;
        log::log(
            log::Level::Info,
            &format!(
                "{} triangles don't fit the 1 MiB command limit; simplifying to {target}",
                out.f.len()
            ),
        );
        simplify(&mut out, target.max(100));
        args = mesh_feature(&out, &name);
    }
    let added = commands::run("timeline.add", &args)?;
    let added: Value = serde_json::from_str(&added).unwrap_or(Value::Null);
    let mut result = serde_json::to_value(&report).unwrap();
    let o = result.as_object_mut().unwrap();
    o.insert("source_body".into(), json!(body));
    o.insert("input_triangles".into(), json!(mesh.f.len()));
    o.insert("triangles".into(), json!(out.f.len()));
    o.insert("feature".into(), added);
    o.insert("version".into(), json!(mesh_solidify::VERSION));
    Ok(result.to_string())
}

fn do_inspect(params: &str) -> Result<String, PluginError> {
    let input: InspectInput =
        serde_json::from_str(params).map_err(|e| invalid("params", e.to_string()))?;
    if !(input.scale.is_finite() && input.scale > 0.0) {
        return Err(invalid("scale", "must be greater than 0"));
    }
    let body = pick_body(&input.body)?;
    let mut mesh = read_body(&body)?;
    mesh.scale(input.scale);
    let mut v = serde_json::to_value(inspect(&mesh)).unwrap();
    v.as_object_mut()
        .unwrap()
        .insert("body".into(), json!(body));
    Ok(v.to_string())
}



#[derive(Deserialize)]
#[serde(default)]
struct BodyInput { body: String }
impl Default for BodyInput { fn default() -> Self { Self { body: String::new() } } }

#[derive(Deserialize)]
#[serde(default)]
struct ScaleInput { body: String, factor: f64, fx: f64, fy: f64, fz: f64, uniform: bool }
impl Default for ScaleInput { fn default() -> Self { Self { body:String::new(), factor:1.0, fx:1.0, fy:1.0, fz:1.0, uniform:true } } }

#[derive(Deserialize)]
#[serde(default)]
struct SplitInput { body:String, axis:String, offset:f64 }
impl Default for SplitInput { fn default() -> Self { Self { body:String::new(), axis:"XY".into(), offset:0.0 } } }

#[derive(Deserialize)]
#[serde(default)]
struct DeleteFaceInput { body:String, x:f64, y:f64, z:f64, heal:bool }
impl Default for DeleteFaceInput { fn default() -> Self { Self { body:String::new(), x:0.0,y:0.0,z:0.0,heal:true } } }

#[derive(Deserialize)]
#[serde(default)]
struct PressPullInput { body:String, x:f64, y:f64, z:f64, distance:f64 }
impl Default for PressPullInput { fn default() -> Self { Self { body:String::new(), x:0.0,y:0.0,z:0.0,distance:1.0 } } }

#[derive(Deserialize)]
#[serde(default)]
struct MoveInput { body:String, dx:f64, dy:f64, dz:f64 }
impl Default for MoveInput { fn default() -> Self { Self { body:String::new(),dx:0.0,dy:0.0,dz:0.0 } } }

fn do_convert(params:&str)->Result<String,PluginError>{
    let input:BodyInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    let body=pick_body(&input.body)?;
    let out=run("mesh.convert",&json!({"body":body.clone(),"tolerance":0.05}))?;
    Ok(json!({"operation":"mesh.convert","body":body,"result":out}).to_string())
}
fn do_scale(params:&str)->Result<String,PluginError>{
    let input:ScaleInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    let body=pick_body(&input.body)?;
    if input.uniform {
        if !input.factor.is_finite() || input.factor<=0.0 { return Err(invalid("factor","must be greater than 0")); }
        let out=run("solid.scale",&json!({"bodies":[body.clone()],"factor":input.factor}))?;
        Ok(json!({"operation":"solid.scale","body":body,"factor":input.factor,"result":out}).to_string())
    } else {
        let fs=[input.fx,input.fy,input.fz];
        if fs.iter().any(|v| !v.is_finite() || *v<=0.0) { return Err(invalid("factors","all axis factors must be greater than 0")); }
        let out=run("solid.scale",&json!({"bodies":[body.clone()],"factors":fs}))?;
        Ok(json!({"operation":"solid.scale","body":body,"factors":fs,"result":out}).to_string())
    }
}
fn do_split(params:&str)->Result<String,PluginError>{
    let input:SplitInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    let body=pick_body(&input.body)?;
    if !input.offset.is_finite(){return Err(invalid("offset","must be finite"));}
    let (normal,origin)=match input.axis.to_ascii_uppercase().as_str(){
        "XY"=>([0.0,0.0,1.0],[0.0,0.0,input.offset]),
        "XZ"=>([0.0,1.0,0.0],[0.0,input.offset,0.0]),
        "YZ"=>([1.0,0.0,0.0],[input.offset,0.0,0.0]),
        _=>return Err(invalid("axis","choose XY, XZ or YZ")),
    };
    let out=run("solid.split_body",&json!({"body":body.clone(),"plane":{"origin":origin,"normal":normal}}))?;
    Ok(json!({"operation":"solid.split_body","body":body,"axis":input.axis,"offset":input.offset,"result":out}).to_string())
}
fn do_delete_face(params:&str)->Result<String,PluginError>{
    let input:DeleteFaceInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    let body=pick_body(&input.body)?;
    let p=[input.x,input.y,input.z];
    if p.iter().any(|v|!v.is_finite()){return Err(invalid("point","coordinates must be finite"));}
    let out=run("face.delete",&json!({"body":body.clone(),"faces":[p],"heal":input.heal}))?;
    Ok(json!({"operation":"face.delete","body":body,"point":p,"result":out}).to_string())
}
fn do_press_pull(params:&str)->Result<String,PluginError>{
    let input:PressPullInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    let body=pick_body(&input.body)?;
    let p=[input.x,input.y,input.z];
    if p.iter().any(|v|!v.is_finite()) || !input.distance.is_finite() || input.distance==0.0 {return Err(invalid("distance","point must be finite and distance must be non-zero"));}
    let out=run("solid.press_pull",&json!({"body":body.clone(),"face":p,"distance":input.distance}))?;
    Ok(json!({"operation":"solid.press_pull","body":body,"point":p,"distance":input.distance,"result":out}).to_string())
}
fn do_move(params:&str)->Result<String,PluginError>{
    let input:MoveInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    let body=pick_body(&input.body)?;
    let v=[input.dx,input.dy,input.dz];
    if v.iter().any(|x|!x.is_finite()){return Err(invalid("translation","coordinates must be finite"));}
    let out=run("solid.move",&json!({"bodies":[body.clone()],"translate":v}))?;
    Ok(json!({"operation":"solid.move","body":body,"translation":v,"result":out}).to_string())
}



#[derive(Deserialize)]
#[serde(default)]
struct CleanInput { body:String, area_epsilon:f64 }
impl Default for CleanInput { fn default()->Self { Self { body:String::new(), area_epsilon:1e-10 } } }

#[derive(Deserialize)]
#[serde(default)]
struct SmoothInput { body:String, iterations:usize, factor:f64, preserve_boundary:bool }
impl Default for SmoothInput { fn default()->Self { Self { body:String::new(), iterations:5, factor:0.35, preserve_boundary:true } } }

#[derive(Deserialize)]
#[serde(default)]
struct DecimateInput { body:String, target_faces:usize }
impl Default for DecimateInput { fn default()->Self { Self { body:String::new(), target_faces:1000 } } }

#[derive(Deserialize)]
#[serde(default)]
struct OffsetInput { body:String, distance:f64 }
impl Default for OffsetInput { fn default()->Self { Self { body:String::new(), distance:0.5 } } }

fn add_mesh_result(mesh:&Mesh, source:&str, suffix:&str)->Result<String,PluginError>{
    if mesh.f.is_empty() { return Err(failed("operation produced an empty mesh")); }
    let name=format!("{source} {suffix}");
    let args=mesh_feature(mesh,&name);
    if args.len()>MAX_ARGS { return Err(failed("result is too large for SolidCraft's 1 MiB command limit; reduce the mesh first")); }
    let added=commands::run("timeline.add",&args)?;
    let feature:Value=serde_json::from_str(&added).unwrap_or(Value::Null);
    Ok(json!({"source_body":source,"new_body":name,"vertices":mesh.v.len(),"triangles":mesh.f.len(),"feature":feature}).to_string())
}

// Native triangle cleanup: removes zero-area and duplicate faces, then unused vertices.
fn do_clean_mesh(params:&str)->Result<String,PluginError>{
    let input:CleanInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    if !input.area_epsilon.is_finite() || input.area_epsilon<0.0 { return Err(invalid("area_epsilon","must be zero or greater")); }
    let body=pick_body(&input.body)?;
    let original=read_body(&body)?;
    let mut out=Mesh { v:original.v.clone(), f:Vec::with_capacity(original.f.len()) };
    let mut seen=std::collections::HashSet::<[u32;3]>::new();
    let mut removed_degenerate=0usize;
    let mut removed_duplicate=0usize;
    for t in &original.f {
        let a=out.v[t[0] as usize]; let b=out.v[t[1] as usize]; let c=out.v[t[2] as usize];
        let area2=mesh_solidify::mesh::norm(mesh_solidify::mesh::cross(mesh_solidify::mesh::sub(b,a),mesh_solidify::mesh::sub(c,a)));
        if !area2.is_finite() || area2*0.5<=input.area_epsilon { removed_degenerate+=1; continue; }
        let mut key=*t; key.sort_unstable();
        if !seen.insert(key) { removed_duplicate+=1; continue; }
        out.f.push(*t);
    }
    out.compact();
    let mut result:Value=serde_json::from_str(&add_mesh_result(&out,&body,"cleaned")?).unwrap();
    let o=result.as_object_mut().unwrap();
    o.insert("removed_degenerate_faces".into(),json!(removed_degenerate));
    o.insert("removed_duplicate_faces".into(),json!(removed_duplicate));
    o.insert("input_triangles".into(),json!(original.f.len()));
    Ok(result.to_string())
}

// Laplacian vertex smoothing performed directly on the triangle mesh.
fn do_smooth_mesh(params:&str)->Result<String,PluginError>{
    let input:SmoothInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    if input.iterations==0 || input.iterations>100 { return Err(invalid("iterations","choose 1 to 100")); }
    if !input.factor.is_finite() || input.factor<=0.0 || input.factor>1.0 { return Err(invalid("factor","must be greater than 0 and at most 1")); }
    let body=pick_body(&input.body)?;
    let mut out=read_body(&body)?;
    let mut boundary=vec![false;out.v.len()];
    if input.preserve_boundary {
        let mut edges=std::collections::HashMap::<(u32,u32),usize>::new();
        for t in &out.f { for k in 0..3 { let a=t[k]; let b=t[(k+1)%3]; let e=if a<b {(a,b)} else {(b,a)}; *edges.entry(e).or_insert(0)+=1; } }
        for ((a,b),count) in edges { if count==1 { boundary[a as usize]=true; boundary[b as usize]=true; } }
    }
    for _ in 0..input.iterations {
        let neighbors=out.vertex_neighbors();
        let old=out.v.clone();
        for i in 0..out.v.len() {
            if boundary[i] || neighbors[i].is_empty() { continue; }
            let mut avg=[0.0;3];
            for &j in &neighbors[i] { for k in 0..3 { avg[k]+=old[j as usize][k]; } }
            for k in 0..3 { avg[k]/=neighbors[i].len() as f64; out.v[i][k]=old[i][k]+(avg[k]-old[i][k])*input.factor; }
        }
    }
    add_mesh_result(&out,&body,"smoothed")
}

// Quadric-error mesh reduction, implemented by this plugin's own mesh library.
fn do_decimate_mesh(params:&str)->Result<String,PluginError>{
    let input:DecimateInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    if input.target_faces<100 { return Err(invalid("target_faces","minimum is 100")); }
    let body=pick_body(&input.body)?;
    let mut out=read_body(&body)?;
    let before=out.f.len();
    if input.target_faces<before { simplify(&mut out,input.target_faces); }
    let mut result:Value=serde_json::from_str(&add_mesh_result(&out,&body,"reduced")?).unwrap();
    result.as_object_mut().unwrap().insert("input_triangles".into(),json!(before));
    Ok(result.to_string())
}

// Vertex-normal offset, calculated from the actual input triangles.
fn do_offset_mesh(params:&str)->Result<String,PluginError>{
    let input:OffsetInput=serde_json::from_str(params).map_err(|e|invalid("params",e.to_string()))?;
    if !input.distance.is_finite() { return Err(invalid("distance","must be finite")); }
    let body=pick_body(&input.body)?;
    let mut out=read_body(&body)?;
    let mut normals=vec![[0.0;3];out.v.len()];
    for t in &out.f {
        let a=out.v[t[0] as usize]; let b=out.v[t[1] as usize]; let c=out.v[t[2] as usize];
        let n=mesh_solidify::mesh::cross(mesh_solidify::mesh::sub(b,a),mesh_solidify::mesh::sub(c,a));
        for &i in t { for k in 0..3 { normals[i as usize][k]+=n[k]; } }
    }
    for i in 0..out.v.len() {
        let len=mesh_solidify::mesh::norm(normals[i]);
        if len>1e-15 { for k in 0..3 { out.v[i][k]+=normals[i][k]/len*input.distance; } }
    }
    add_mesh_result(&out,&body,"offset")
}

impl exports::craft::plugin::entry::Guest for Plugin {
    fn execute(entry: String, params: String) -> Result<String, PluginError> {
        match entry.as_str() {
            "solidify" => do_solidify(&params),
            "inspect" => do_inspect(&params),
            "convert" => do_convert(&params),
            "scale" => do_scale(&params),
            "split" => do_split(&params),
            "delete-face" => do_delete_face(&params),
            "press-pull" => do_press_pull(&params),
            "move" => do_move(&params),
            "clean-mesh" => do_clean_mesh(&params),
            "smooth-mesh" => do_smooth_mesh(&params),
            "decimate-mesh" => do_decimate_mesh(&params),
            "offset-mesh" => do_offset_mesh(&params),
            other => Err(invalid("entry", format!("unknown entry {other:?}"))),
        }
    }
}
