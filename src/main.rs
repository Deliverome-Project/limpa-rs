// SPDX-License-Identifier: GPL-3.0-or-later
#![forbid(unsafe_code)]

use limpa_rs::{Model, Protein, fit};
use rayon::prelude::*;
use std::{
    fs::{self, File},
    io::{self, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    time::{Instant, SystemTime, UNIX_EPOCH},
};
// A successful result appears atomically; failed runs leave no output to mistake
// for a complete analysis. Hard-link publication cannot overwrite an existing file.
struct PendingOutput(PathBuf);
impl PendingOutput {
    fn create(destination: &Path) -> io::Result<(Self, File)> {
        let parent = destination.parent().unwrap_or(Path::new("."));
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(io::Error::other)?
            .as_nanos();
        let path = parent.join(format!(".limpa-rs-{}-{nonce}.tmp", std::process::id()));
        let file = File::options().write(true).create_new(true).open(&path)?;
        Ok((Self(path), file))
    }
    fn publish(&self, destination: &Path) -> io::Result<()> {
        fs::hard_link(&self.0, destination)
    }
}
impl Drop for PendingOutput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

fn uint(r: &mut impl Read) -> io::Result<usize> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b) as usize)
}
fn float(r: &mut impl Read) -> io::Result<f64> {
    let mut b = [0; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}
fn floats(r: &mut impl Read, n: usize) -> io::Result<Vec<f64>> {
    (0..n).map(|_| float(r)).collect()
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() != 3 {
        return Err("Usage: limpa-rs INPUT.bin OUTPUT.bin THREADS (positive integer)".into());
    }
    if a[0] == a[1] {
        return Err("input and output paths must differ".into());
    }
    let threads: usize = a[2].parse()?;
    if threads == 0 {
        return Err("THREADS must be positive".into());
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;
    let t = Instant::now();
    let mut r = BufReader::new(File::open(&a[0])?);
    let mut magic = [0; 8];
    r.read_exact(&mut magic)?;
    if &magic != b"LIMPAR01" {
        return Err("invalid input magic/version".into());
    }
    let np = uint(&mut r)?;
    let n = uint(&mut r)?;
    if np == 0 || !(2..=100_000).contains(&n) {
        return Err("invalid dimensions".into());
    }
    let m = Model {
        intercept: float(&mut r)?,
        slope: float(&mut r)?,
        prior_mean: float(&mut r)?,
        prior_sd: float(&mut r)?,
        prior_logfc: float(&mut r)?,
    };
    let destination = Path::new(&a[1]);
    if destination.exists() {
        return Err("output already exists; refusing to overwrite".into());
    }
    let (pending, file) = PendingOutput::create(destination)?;
    let mut out = BufWriter::new(file);
    out.write_all(b"LIMPAO01")?;
    out.write_all(&(np as u32).to_le_bytes())?;
    out.write_all(&(n as u32).to_le_bytes())?;
    let mut max_iter = 0;
    let mut max_gradient = 0.0_f64;
    for begin in (0..np).step_by(64) {
        if begin > 0 && begin % 1024 == 0 {
            eprintln!(
                "completed_proteins={begin}/{np} elapsed_seconds={:.1}",
                t.elapsed().as_secs_f64()
            );
        }
        let mut batch = Vec::new();
        for _ in begin..(begin + 64).min(np) {
            let p = uint(&mut r)?;
            if p == 0 || p > 100_000 || p.checked_mul(n).is_none_or(|v| v > 100_000_000) {
                return Err("invalid peptide dimensions".into());
            }
            let sigma = float(&mut r)?;
            let y = floats(&mut r, p * n)?;
            let start = floats(&mut r, n + p - 1)?;
            batch.push(Protein {
                samples: n,
                peptides: p,
                sigma,
                y,
                start,
            });
        }
        let fits: Vec<_> = batch.par_iter().map(|p| fit(p, m)).collect();
        for (i, f) in fits.into_iter().enumerate() {
            let f = f.map_err(|e| format!("protein {}: {e}", begin + i + 1))?;
            max_iter = max_iter.max(f.iterations);
            max_gradient = max_gradient.max(f.gradient_max);
            for x in f.expression.iter().chain(f.se.iter()) {
                out.write_all(&x.to_le_bytes())?;
            }
        }
    }
    let mut extra = [0];
    if r.read(&mut extra)? != 0 {
        return Err("trailing input data".into());
    }
    out.flush()?;
    out.get_ref().sync_all()?;
    drop(out);
    pending.publish(destination)?;
    eprintln!(
        "proteins={np} samples={n} threads={threads} elapsed_seconds={:.6} max_iterations={max_iter} max_gradient={max_gradient:.3e}",
        t.elapsed().as_secs_f64()
    );
    Ok(())
}
fn main() {
    if let Err(e) = run() {
        eprintln!("limpa-rs: {e}");
        std::process::exit(1);
    }
}
