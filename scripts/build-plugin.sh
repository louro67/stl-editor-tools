#!/bin/sh
# Build the SolidCraft plugin package: dist/mesh-solidify/ and dist/mesh-solidify-plugin.zip.
# Needs the wasm32-unknown-unknown Rust target and wasm-tools.
set -eu
cd "$(dirname "$0")/.."
cargo build --release -p solidcraft-mesh-solidify --target wasm32-unknown-unknown
rm -rf dist/mesh-solidify dist/mesh-solidify-plugin.zip
mkdir -p dist/mesh-solidify
wasm-tools component new target/wasm32-unknown-unknown/release/solidcraft_mesh_solidify.wasm \
  -o dist/mesh-solidify/mesh-solidify.wasm
wasm-tools validate dist/mesh-solidify/mesh-solidify.wasm
cp crates/plugin/package/plugin.json dist/mesh-solidify/
cp LICENSE dist/mesh-solidify/LICENSE
(cd dist && zip -qr mesh-solidify-plugin.zip mesh-solidify)
ls -l dist/mesh-solidify dist/mesh-solidify-plugin.zip
