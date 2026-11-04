use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use std::collections::HashMap;

mod pileup;

/// Windows as 4-bit base masks (n x L). Keep the windows the motif matches (or does not match),
/// and return them with their PSSM (4 x L, share of kept windows allowing each base).
/// A window matches when every position's bases are a subset of the motif's allowed bases.
#[pyfunction]
fn filter_and_pssm<'py>(
    py: Python<'py>,
    masks: PyReadonlyArray2<'py, u8>,
    motif_mask: PyReadonlyArray1<'py, u8>,
    keep_matches: bool,
) -> PyResult<(Bound<'py, PyArray2<u8>>, Bound<'py, PyArray2<f64>>)> {
    let masks = masks.as_array();
    let motif = motif_mask.as_slice()?;
    let (n, len) = (masks.nrows(), masks.ncols());
    let not_allowed: Vec<u8> = motif.iter().map(|m| !m).collect();
    let data = masks.as_slice().expect("contiguous windows");
    let mut kept: Vec<u8> = Vec::with_capacity(data.len() / 4);
    // per position, how many kept windows have each of the 16 masks; the PSSM follows at the end
    let mut hist = vec![0u64; 16 * len];
    let mut n_kept = 0usize;
    for i in 0..n {
        let row = &data[i * len..(i + 1) * len];
        let matches = row.iter().zip(&not_allowed).all(|(b, na)| b & na == 0);
        if matches == keep_matches {
            kept.extend_from_slice(row);
            for (j, &b) in row.iter().enumerate() {
                hist[j * 16 + (b & 15) as usize] += 1;
            }
            n_kept += 1;
        }
    }
    let mut pssm = vec![0f64; 4 * len];
    for j in 0..len {
        for mask in 0..16usize {
            let c = hist[j * 16 + mask];
            for bit in 0..4 {
                if mask >> bit & 1 == 1 {
                    pssm[bit * len + j] += c as f64;
                }
            }
        }
    }
    for v in pssm.iter_mut() {
        *v /= n_kept.max(1) as f64;
    }
    let kept = numpy::ndarray::Array2::from_shape_vec((n_kept, len), kept).unwrap();
    let pssm = numpy::ndarray::Array2::from_shape_vec((4, len), pssm).unwrap();
    Ok((kept.into_pyarray_bound(py), pssm.into_pyarray_bound(py)))
}

/// Filter a plain modkit pileup and write it as per (bin, mod type) lz4 Arrow IPC files.
/// Returns {(bin, mod_type): [files]}.
#[pyfunction]
#[pyo3(signature = (path, contig_bin, out_dir, min_coverage, threshold, min_mods=50, min_mod_frequency=0.0001, adjacency_distance=8, threads=1, flush_rows=500_000))]
fn partition_pileup(
    py: Python<'_>, path: &str, contig_bin: HashMap<String, String>, out_dir: &str, min_coverage: i64, threshold: f64,
    min_mods: usize, min_mod_frequency: f64, adjacency_distance: i64, threads: usize, flush_rows: usize,
) -> PyResult<HashMap<(String, String), Vec<String>>> {
    let settings = pileup::Settings { min_coverage, threshold, min_mods, min_mod_frequency, adjacency_distance, flush_rows, max_buffered: 8 * flush_rows };
    py.allow_threads(|| pileup::partition(path, &contig_bin, &settings, out_dir, threads.max(1)))
        .map_err(pyo3::exceptions::PyValueError::new_err)
}

#[pymodule]
fn nanomotif_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(filter_and_pssm, m)?)?;
    m.add_function(wrap_pyfunction!(partition_pileup, m)?)?;
    Ok(())
}
