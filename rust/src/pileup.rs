//! Streaming pass over a modkit pileup (bedMethyl): filter, then write per (bin, mod type) Arrow IPC files.
//!
//! Same filters, in the same order, as nanomotif.dataload with the plain-pileup path of the builder:
//! 1. coverage (Nvalid_cov >= min_coverage) and binned contigs only
//! 2. minimum modification: keep a (contig, mod type) when it has more than `min_mods` positions with
//!    fraction >= threshold and that share of its positions is above `min_mod_frequency`
//! 3. adjacency: a position >= threshold is dropped when another position of the same contig and strand
//!    (any mod type) within `adjacency_distance` has a higher fraction
//! Rows of a contig must be contiguous in the file (modkit output); this is checked.

use arrow_array::{ArrayRef, Float64Array, Int64Array, RecordBatch, StringArray};
use arrow_ipc::writer::{FileWriter, IpcWriteOptions};
use arrow_ipc::CompressionType;
use arrow_schema::{DataType, Field, Schema};
use memchr::memchr;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::{Arc, Mutex};

const CHUNK: usize = 16 << 20;

pub struct Settings {
    pub min_coverage: i64,
    pub threshold: f64,
    pub min_mods: usize,
    pub min_mod_frequency: f64,
    pub adjacency_distance: i64,
    pub flush_rows: usize,
    pub max_buffered: usize,
}

struct Row {
    position: i64,
    mod_type: u16, // index into the thread's mod type list
    minus: bool,
    fraction: f64,
    coverage: i64,
}

/// Buffered rows of one (bin, mod type), flushed to numbered IPC files.
struct Sink {
    contig: Vec<String>,
    position: Vec<i64>,
    strand: Vec<&'static str>,
    fraction: Vec<f64>,
    coverage: Vec<i64>,
}

impl Sink {
    fn new() -> Self {
        Sink { contig: vec![], position: vec![], strand: vec![], fraction: vec![], coverage: vec![] }
    }
}

fn schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("contig", DataType::Utf8, false),
        Field::new("position", DataType::Int64, false),
        Field::new("mod_type", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("fraction_mod", DataType::Float64, false),
        Field::new("Nvalid_cov", DataType::Int64, false),
        Field::new("bin", DataType::Utf8, false),
    ]))
}

fn write_ipc(path: &str, bin: &str, mod_type: &str, sink: &mut Sink) -> std::io::Result<()> {
    let n = sink.position.len();
    let columns: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(std::mem::take(&mut sink.contig))),
        Arc::new(Int64Array::from(std::mem::take(&mut sink.position))),
        Arc::new(StringArray::from(vec![mod_type; n])),
        Arc::new(StringArray::from(std::mem::take(&mut sink.strand))),
        Arc::new(Float64Array::from(std::mem::take(&mut sink.fraction))),
        Arc::new(Int64Array::from(std::mem::take(&mut sink.coverage))),
        Arc::new(StringArray::from(vec![bin; n])),
    ];
    let batch = RecordBatch::try_new(schema(), columns).map_err(std::io::Error::other)?;
    let options = IpcWriteOptions::default()
        .try_with_compression(Some(CompressionType::LZ4_FRAME))
        .map_err(std::io::Error::other)?;
    let mut writer = FileWriter::try_new_with_options(File::create(path)?, &schema(), options).map_err(std::io::Error::other)?;
    writer.write(&batch).map_err(std::io::Error::other)?;
    writer.finish().map_err(std::io::Error::other)?;
    Ok(())
}

fn field(line: &[u8], start: usize) -> (&[u8], usize) {
    match memchr(b'\t', &line[start..]) {
        Some(i) => (&line[start..start + i], start + i + 1),
        None => (&line[start..], line.len()),
    }
}

fn parse_i64(b: &[u8]) -> Option<i64> {
    std::str::from_utf8(b).ok()?.parse().ok()
}

/// Byte ranges that each start at the first row of a contig.
pub fn split_ranges(data: &[u8], n_parts: usize) -> Vec<(usize, usize)> {
    let size = data.len();
    let contig_at = |offset: usize| -> (usize, &[u8]) {
        let start = if offset == 0 { 0 } else {
            match memchr(b'\n', &data[offset..]) { Some(i) => offset + i + 1, None => size }
        };
        let end = memchr(b'\t', &data[start..]).map(|i| start + i).unwrap_or(size);
        (start, &data[start..end])
    };
    let mut cuts = vec![0usize];
    for i in 1..n_parts {
        let (start, name) = contig_at(size * i / n_parts);
        if name.is_empty() || start <= *cuts.last().unwrap() {
            continue;
        }
        // rows of a contig are contiguous: binary search for where it ends
        let (mut lo, mut hi) = (start, size);
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if contig_at(mid).1 == name { lo = mid } else { hi = mid }
        }
        let pos = if hi < size { contig_at(hi).0 } else { size };
        if pos > *cuts.last().unwrap() && pos < size {
            cuts.push(pos);
        }
    }
    cuts.push(size);
    cuts.windows(2).filter(|w| w[1] > w[0]).map(|w| (w[0], w[1])).collect()
}

struct Output {
    dir: String,
    part: usize,
    files: HashMap<(String, String), Vec<String>>,
}

/// One thread: stream a byte range contig by contig.
fn process_range(
    path: &str, range: (usize, usize), contig_bin: &HashMap<String, String>, s: &Settings, out: &Mutex<Output>,
) -> Result<Vec<String>, String> {
    use std::io::{Read, Seek, SeekFrom};
    let mut file = File::open(path).map_err(|e| e.to_string())?;
    file.seek(SeekFrom::Start(range.0 as u64)).map_err(|e| e.to_string())?;
    let mut remaining = range.1 - range.0;
    let mut chunk: Vec<u8> = Vec::new();
    let mut buf = vec![0u8; CHUNK];
    let mut mod_types: Vec<String> = vec![];
    let mut sinks: HashMap<(String, u16), Sink> = HashMap::new();
    let mut contigs_seen: Vec<String> = vec![];
    let mut rows: Vec<Row> = vec![];
    let mut current: Option<Vec<u8>> = None;
    let mut seen: HashSet<Vec<u8>> = HashSet::new();

    let flush = |sinks: &mut HashMap<(String, u16), Sink>, mod_types: &Vec<String>, force: bool| -> Result<(), String> {
        // bounded memory: flush everything once the thread holds max_buffered rows in total
        let force = force || sinks.values().map(|x| x.position.len()).sum::<usize>() >= s.max_buffered;
        for ((bin, m), sink) in sinks.iter_mut() {
            if sink.position.is_empty() || (!force && sink.position.len() < s.flush_rows) {
                continue;
            }
            let path = {
                let mut o = out.lock().unwrap();
                o.part += 1;
                let p = format!("{}/{}.arrow", o.dir, o.part);
                o.files.entry((bin.clone(), mod_types[*m as usize].clone())).or_default().push(p.clone());
                p
            };
            write_ipc(&path, bin, &mod_types[*m as usize], sink).map_err(|e| e.to_string())?;
        }
        Ok(())
    };

    let mut finish_contig = |name: &[u8], rows: &mut Vec<Row>, sinks: &mut HashMap<(String, u16), Sink>, n_mod_types: usize| {
        let contig = String::from_utf8_lossy(name).to_string();
        let bin = &contig_bin[&contig];
        // 2. minimum modification per mod type
        let mut n = vec![0usize; n_mod_types];
        let mut k = vec![0usize; n_mod_types];
        for r in rows.iter() {
            n[r.mod_type as usize] += 1;
            if r.fraction >= s.threshold { k[r.mod_type as usize] += 1; }
        }
        let keep_mod: Vec<bool> = (0..n_mod_types)
            .map(|m| n[m] > 0 && (k[m] as f64 / n[m] as f64) > s.min_mod_frequency && k[m] > s.min_mods)
            .collect();
        rows.retain(|r| keep_mod[r.mod_type as usize]);
        // 3. adjacency among the high rows, per strand, sorted by position
        let mut keep = vec![true; rows.len()];
        for minus in [false, true] {
            let mut high: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].minus == minus && rows[i].fraction >= s.threshold).collect();
            high.sort_by_key(|&i| rows[i].position);
            for a in 0..high.len() {
                let mut b = a + 1;
                while b < high.len() && rows[high[b]].position - rows[high[a]].position <= s.adjacency_distance {
                    let (ra, rb) = (&rows[high[a]], &rows[high[b]]);
                    if rb.fraction > ra.fraction { keep[high[a]] = false; }
                    if ra.fraction > rb.fraction { keep[high[b]] = false; }
                    b += 1;
                }
            }
        }
        for (i, r) in rows.iter().enumerate() {
            if !keep[i] { continue; }
            let sink = sinks.entry((bin.clone(), r.mod_type)).or_insert_with(Sink::new);
            sink.contig.push(contig.clone());
            sink.position.push(r.position);
            sink.strand.push(if r.minus { "-" } else { "+" });
            sink.fraction.push(r.fraction);
            sink.coverage.push(r.coverage);
        }
        rows.clear();
    };

    loop {
    // next chunk of whole lines (the partial last line is carried over)
    if remaining == 0 && chunk.is_empty() { break; }
    let take = remaining.min(CHUNK);
    let n_read = file.read(&mut buf[..take]).map_err(|e| e.to_string())?;
    remaining -= n_read;
    chunk.extend_from_slice(&buf[..n_read]);
    let usable = if remaining == 0 || n_read == 0 { chunk.len() } else {
        match memchr::memrchr(b'\n', &chunk) { Some(i) => i + 1, None => continue }
    };
    let block: Vec<u8> = chunk.drain(..usable).collect();
    if n_read == 0 { remaining = 0; }
    let mut pos = 0usize;
    while pos < block.len() {
        let end = memchr(b'\n', &block[pos..]).map(|i| pos + i).unwrap_or(block.len());
        let line = &block[pos..end];
        pos = end + 1;
        if line.is_empty() || line[0] == b'#' { continue; }
        let (contig, p) = field(line, 0);
        if current.as_deref() != Some(contig) {
            if let Some(prev) = current.take() {
                if contig_bin.contains_key(std::str::from_utf8(&prev).unwrap_or("")) {
                    finish_contig(&prev, &mut rows, &mut sinks, mod_types.len());
                    flush(&mut sinks, &mod_types, false)?;
                }
                seen.insert(prev);
            }
            if seen.contains(contig) {
                return Err(format!("pileup rows of contig {} are not contiguous", String::from_utf8_lossy(contig)));
            }
            current = Some(contig.to_vec());
            contigs_seen.push(String::from_utf8_lossy(contig).to_string());
        }
        if !contig_bin.contains_key(std::str::from_utf8(contig).unwrap_or("")) { continue; }
        let (start, p) = field(line, p);
        let (_, p) = field(line, p);
        let (code, p) = field(line, p);
        let (_, p) = field(line, p);
        let (strand, p) = field(line, p);
        let (_, p) = field(line, p);
        let (_, p) = field(line, p);
        let (_, p) = field(line, p);
        let (n_valid, p) = field(line, p);
        let (percent, _) = field(line, p);
        let coverage = parse_i64(n_valid).ok_or("bad Nvalid_cov")?;
        if coverage < s.min_coverage { continue; }
        let code = std::str::from_utf8(code).map_err(|e| e.to_string())?;
        let m = match mod_types.iter().position(|x| x == code) {
            Some(i) => i,
            None => { mod_types.push(code.to_string()); mod_types.len() - 1 }
        };
        let percent: f64 = std::str::from_utf8(percent).ok().and_then(|x| x.parse().ok()).ok_or("bad fraction")?;
        rows.push(Row {
            position: parse_i64(start).ok_or("bad start")?,
            mod_type: m as u16,
            minus: strand == b"-",
            fraction: percent / 100.0,
            coverage,
        });
    }
    }
    if let Some(prev) = current.take() {
        if contig_bin.contains_key(std::str::from_utf8(&prev).unwrap_or("")) {
            finish_contig(&prev, &mut rows, &mut sinks, mod_types.len());
        }
    }
    flush(&mut sinks, &mod_types, true)?;
    Ok(contigs_seen)
}

/// Returns {(bin, mod_type): [ipc files]}.
pub fn partition(
    path: &str, contig_bin: &HashMap<String, String>, s: &Settings, out_dir: &str, threads: usize,
) -> Result<HashMap<(String, String), Vec<String>>, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let data = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| e.to_string())?;
    let ranges = split_ranges(&data, threads * 4);
    drop(data);
    let out = Mutex::new(Output { dir: out_dir.to_string(), part: 0, files: HashMap::new() });
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().map_err(|e| e.to_string())?;
    let seen: Vec<Vec<String>> = pool.install(|| {
        ranges.par_iter().map(|&r| process_range(path, r, contig_bin, s, &out)).collect::<Result<_, _>>()
    })?;
    // a contig split across ranges means its rows are not contiguous
    let mut all = HashSet::new();
    for c in seen.into_iter().flatten() {
        if !all.insert(c.clone()) {
            return Err(format!("pileup rows of contig {} are not contiguous", c));
        }
    }
    Ok(out.into_inner().unwrap().files)
}
