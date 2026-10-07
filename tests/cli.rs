// SPDX-License-Identifier: GPL-3.0-or-later
#![forbid(unsafe_code)]

use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Workspace(PathBuf);
impl Workspace {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "limpa-rs-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn run(&self, input: &str, output: &str, threads: &str) -> Output {
        Command::new(env!("CARGO_BIN_EXE_limpa-rs"))
            .arg(self.0.join(input))
            .arg(self.0.join(output))
            .arg(threads)
            .output()
            .unwrap()
    }
}
impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn fixture(proteins: u32) -> Vec<u8> {
    let mut bytes = b"LIMPAR01".to_vec();
    bytes.extend(proteins.to_le_bytes());
    bytes.extend(3u32.to_le_bytes());
    for x in [-11.0f64, 0.75, 18.0, 3.0, 2.0] {
        bytes.extend(x.to_le_bytes());
    }
    for i in 0..proteins {
        let shift = f64::from(i) * 0.01;
        bytes.extend(2u32.to_le_bytes());
        bytes.extend(0.4f64.to_le_bytes());
        for x in [18.0, 19.0, f64::NAN, 20.0, f64::NAN, 21.0] {
            bytes.extend((x + shift).to_le_bytes());
        }
        for x in [
            19.0 + shift,
            18.25 + shift,
            19.7269230769231 + shift,
            -0.507692307692306,
        ] {
            bytes.extend(x.to_le_bytes());
        }
    }
    bytes
}
#[test]
fn worker_count_preserves_bytes_and_order_across_batches() {
    let w = Workspace::new();
    fs::write(w.0.join("input"), fixture(70)).unwrap();
    for (threads, name) in [("1", "serial"), ("4", "parallel")] {
        let out = w.run("input", name, threads);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let a = fs::read(w.0.join("serial")).unwrap();
    let b = fs::read(w.0.join("parallel")).unwrap();
    assert_eq!(a, b);
    assert_eq!(&a[..8], b"LIMPAO01");
    assert_eq!(a.len(), 16 + 70 * 6 * 8);
    let first = f64::from_le_bytes(a[16..24].try_into().unwrap());
    let last = f64::from_le_bytes(a[16 + 69 * 48..24 + 69 * 48].try_into().unwrap());
    assert!(
        last > first + 0.5,
        "fixture must contain distinguishable protein fits"
    );
    for row in a[16..].as_chunks::<48>().0 {
        for (j, cell) in row.as_chunks::<8>().0.iter().enumerate() {
            let value = f64::from_le_bytes(*cell);
            assert!(value.is_finite());
            if j >= 3 {
                assert!(value > 0.0);
            }
        }
    }
}
#[test]
fn malformed_and_invalid_numeric_inputs_fail() {
    let w = Workspace::new();
    let good = fixture(1);
    let mut bad_magic = good.clone();
    bad_magic[0] = 0;
    let mut bad_model = good.clone();
    bad_model[16..24].copy_from_slice(&f64::NAN.to_le_bytes());
    for (i, bytes) in [
        vec![],
        bad_magic,
        good[..good.len() - 1].to_vec(),
        bad_model,
    ]
    .into_iter()
    .enumerate()
    {
        fs::write(w.0.join("input"), bytes).unwrap();
        let out = w.run("input", &format!("out{i}"), "1");
        assert!(!out.status.success());
        assert!(!w.0.join(format!("out{i}")).exists());
        assert!(String::from_utf8_lossy(&out.stderr).contains("limpa-rs:"));
    }
    fs::write(w.0.join("input"), good).unwrap();
    assert!(!w.run("input", "zero-threads", "0").status.success());
}
#[test]
fn existing_output_and_input_are_never_overwritten() {
    let w = Workspace::new();
    let input = fixture(1);
    fs::write(w.0.join("input"), &input).unwrap();
    fs::write(w.0.join("out"), b"keep existing results").unwrap();
    assert!(!w.run("input", "out", "1").status.success());
    assert_eq!(fs::read(w.0.join("out")).unwrap(), b"keep existing results");
    assert!(!w.run("input", "input", "1").status.success());
    assert_eq!(fs::read(w.0.join("input")).unwrap(), input);
}

#[test]
fn late_failure_never_publishes_partial_results() {
    let w = Workspace::new();
    let mut input = fixture(70);
    input.push(42); // detected only after both batches have been computed
    fs::write(w.0.join("input"), input).unwrap();
    assert!(!w.run("input", "out", "4").status.success());
    assert!(!w.0.join("out").exists());
    assert_eq!(fs::read_dir(&w.0).unwrap().count(), 1);
}
