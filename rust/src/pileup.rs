//! Streaming pass over a modkit pileup (bedMethyl): filter, then write per (bin, mod type) Arrow IPC files.
//!
//! Same filters, in the same order, as nanomotif.dataload with the plain-pileup path of the builder:
//! 1. coverage (Nvalid_cov >= min_coverage) and binned contigs only
//! 2. minimum modification: keep a (contig, mod type) when it has more than `min_mods` positions with
//!    fraction >= threshold and that share of its positions is above `min_mod_frequency`
//! 3. adjacency: a position >= threshold is dropped when another position of the same contig and strand
//!    (any mod type) within `adjacency_distance` has a higher fraction
//! Rows of a contig must be contiguous in the file (modkit output); this is checked. Plain files are
//! read in parallel byte ranges; gzipped (bgzip) files as one decompressed stream.

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
/// Rows hold an index into the sink's contig names (rows of a contig arrive together), so a row takes
/// ~30 bytes; the strings are built only when the file is written.
struct Sink {
    names: Vec<String>,
    contig: Vec<u32>,
    position: Vec<i64>,
    minus: Vec<bool>,
    fraction: Vec<f64>,
    coverage: Vec<i64>,
}

impl Sink {
    fn new() -> Self {
        Sink { names: vec![], contig: vec![], position: vec![], minus: vec![], fraction: vec![], coverage: vec![] }
    }

    fn push(&mut self, contig: &str, position: i64, minus: bool, fraction: f64, coverage: i64) {
        if self.names.last().map(|n| n.as_str()) != Some(contig) {
            self.names.push(contig.to_string());
        }
        self.contig.push((self.names.len() - 1) as u32);
        self.position.push(position);
        self.minus.push(minus);
        self.fraction.push(fraction);
        self.coverage.push(coverage);
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
        Arc::new(StringArray::from_iter_values(sink.contig.iter().map(|&i| sink.names[i as usize].as_str()))),
        Arc::new(Int64Array::from(std::mem::take(&mut sink.position))),
        Arc::new(StringArray::from(vec![mod_type; n])),
        Arc::new(StringArray::from_iter_values(sink.minus.iter().map(|&m| if m { "-" } else { "+" }))),
        Arc::new(Float64Array::from(std::mem::take(&mut sink.fraction))),
        Arc::new(Int64Array::from(std::mem::take(&mut sink.coverage))),
        Arc::new(StringArray::from(vec![bin; n])),
    ];
    sink.contig.clear();
    sink.minus.clear();
    sink.names.clear();
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

fn is_gzip(path: &str) -> bool {
    path.ends_with(".gz") || path.ends_with(".bgz")
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
    use std::io::{Seek, SeekFrom};
    let mut f = File::open(path).map_err(|e| e.to_string())?;
    f.seek(SeekFrom::Start(range.0 as u64)).map_err(|e| e.to_string())?;
    process_stream(f, range.1 - range.0, contig_bin, s, out)
}

/// Process up to `remaining` bytes of whole-contig pileup rows from a reader.
fn process_stream<R: std::io::Read>(
    mut file: R, mut remaining: usize, contig_bin: &HashMap<String, String>, s: &Settings, out: &Mutex<Output>,
) -> Result<Vec<String>, String> {
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
            sink.push(&contig, r.position, r.minus, r.fraction, r.coverage);
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

/// A reader over byte blocks received from a channel (one gzip worker's share of the stream).
struct ChannelReader {
    rx: std::sync::mpsc::Receiver<Vec<u8>>,
    block: Vec<u8>,
    pos: usize,
}

impl std::io::Read for ChannelReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        while self.pos == self.block.len() {
            match self.rx.recv() {
                Ok(block) => { self.block = block; self.pos = 0; }
                Err(_) => return Ok(0), // sender done
            }
        }
        let n = buf.len().min(self.block.len() - self.pos);
        buf[..n].copy_from_slice(&self.block[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Start of the last contig's first row in a block of whole lines (rows of a contig are contiguous).
fn last_contig_start(block: &[u8]) -> usize {
    let name_at = |offset: usize| -> (usize, &[u8]) {
        let start = if offset == 0 { 0 } else {
            match memchr(b'\n', &block[offset - 1..]) { Some(i) => offset + i, None => block.len() }
        };
        let end = memchr(b'\t', &block[start..]).map(|i| start + i).unwrap_or(block.len());
        (start, &block[start..end])
    };
    let last_line = match memchr::memrchr(b'\n', &block[..block.len().saturating_sub(1)]) { Some(i) => i + 1, None => 0 };
    let (_, last) = name_at(last_line);
    let (mut lo, mut hi) = (0usize, last_line);
    if name_at(0).1 == last { return 0; }
    // invariant: the line starting at/after lo is another contig, the one at hi is the last contig
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        let (start, name) = name_at(mid);
        if start >= hi || name == last { hi = mid } else { lo = mid }
    }
    name_at(hi).0
}

/// Gzipped (bgzip) pileup: one thread decompresses and hands whole-contig stretches to `threads` workers
/// (switching worker only at a contig boundary, after `SEGMENT` bytes), each running process_stream.
fn partition_gzip(
    path: &str, contig_bin: &HashMap<String, String>, s: &Settings, out: &Mutex<Output>, threads: usize,
) -> Result<Vec<Vec<String>>, String> {
    use std::io::Read;
    const SEGMENT: usize = 64 << 20;
    let mut decoder = flate2::read::MultiGzDecoder::new(std::io::BufReader::with_capacity(
        CHUNK, File::open(path).map_err(|e| e.to_string())?));
    std::thread::scope(|scope| {
        let mut senders = Vec::new();
        let mut handles = Vec::new();
        for _ in 0..threads {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<u8>>(2);
            senders.push(tx);
            handles.push(scope.spawn(move || {
                process_stream(ChannelReader { rx, block: Vec::new(), pos: 0 }, usize::MAX, contig_bin, s, out)
            }));
        }
        let mut worker = 0usize;
        let mut sent = 0usize;
        let mut carry: Vec<u8> = Vec::new();
        let mut buf = vec![0u8; CHUNK];
        let mut error = None;
        loop {
            let n = match decoder.read(&mut buf) { Ok(n) => n, Err(e) => { error = Some(e.to_string()); break } };
            if n == 0 { break; }
            carry.extend_from_slice(&buf[..n]);
            let whole = match memchr::memrchr(b'\n', &carry) { Some(i) => i + 1, None => continue };
            let mut block: Vec<u8> = carry.drain(..whole).collect();
            if sent >= SEGMENT {
                let cut = last_contig_start(&block);
                if cut > 0 {
                    let rest = block.split_off(cut);
                    if senders[worker].send(block).is_err() { break; }
                    block = rest;
                    worker = (worker + 1) % threads;
                    sent = 0;
                }
            }
            sent += block.len();
            if senders[worker].send(block).is_err() { break; }
        }
        if !carry.is_empty() {
            let _ = senders[worker].send(std::mem::take(&mut carry));
        }
        drop(senders);
        let results: Vec<Result<Vec<String>, String>> = handles.into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err("worker panicked".to_string()))).collect();
        if let Some(e) = error { return Err(e); }
        results.into_iter().collect()
    })
}

/// Returns {(bin, mod_type): [ipc files]}.
pub fn partition(
    path: &str, contig_bin: &HashMap<String, String>, s: &Settings, out_dir: &str, threads: usize,
) -> Result<HashMap<(String, String), Vec<String>>, String> {
    let out_files = Mutex::new(Output { dir: out_dir.to_string(), part: 0, files: HashMap::new() });
    if is_gzip(path) {
        let seen = partition_gzip(path, contig_bin, s, &out_files, threads)?;
        check_contiguous(seen)?;
        return Ok(out_files.into_inner().unwrap().files);
    }
    let ranges = {
        let file = File::open(path).map_err(|e| e.to_string())?;
        let data = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| e.to_string())?;
        split_ranges(&data, threads * 4)
    };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().map_err(|e| e.to_string())?;
    let seen: Vec<Vec<String>> = pool.install(|| {
        ranges.par_iter().map(|&r| process_range(path, r, contig_bin, s, &out_files)).collect::<Result<_, _>>()
    })?;
    check_contiguous(seen)?;
    Ok(out_files.into_inner().unwrap().files)
}

/// A contig seen by two readers (ranges or gzip workers) means its rows are not contiguous.
fn check_contiguous(seen: Vec<Vec<String>>) -> Result<(), String> {
    let mut all = HashSet::new();
    for c in seen.into_iter().flatten() {
        if !all.insert(c.clone()) {
            return Err(format!("pileup rows of contig {} are not contiguous", c));
        }
    }
    Ok(())
}
