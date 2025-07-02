use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    io::Write,
    path::{Path, PathBuf},
    sync::mpsc::{self, Receiver, Sender},
    thread::{self, JoinHandle},
    time::Instant,
};

use log::trace;
use serde::{Deserialize, Serialize};

use crate::{
    math::finite_field::{FFRep, FiniteField as FF},
    matrices::mat_trait::RankMatrix,
};

use super::{
    sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    sparse_vec::SparseVector,
};

fn worker_thread_matrix_loop(
    mut mat: SparseFFMatrix,
    sender: Sender<PivotizeMessage>,
    receiver: Receiver<PivotizeMessage>,
    id: usize,
) -> SparseFFMatrix {
    loop {
        let message = receiver.recv();
        if message.is_err() {
            log::error!("Could not retrieve message, connection with coordinator is closed.");
            panic!()
        }
        let message = message.unwrap();
        match message {
            PivotizeMessage::RowNotFound
            | PivotizeMessage::RowRetrieved(_)
            | PivotizeMessage::NNZ(Some(_))
            | PivotizeMessage::NumRows(Some(_))
            | PivotizeMessage::NextPivotFound(_)
            | PivotizeMessage::Completed
            | PivotizeMessage::VerfifyUpperTriangular(Some(_)) => {
                log::error!("Coordinator should not be sending these messages.");
            }
            PivotizeMessage::Error => {
                log::error!("Coordinator error, worker needs to bail.");
                panic!()
            }
            PivotizeMessage::FindNextPivot(previous_pivot) => {
                let next_pivot = mat.find_next_pivot(previous_pivot);
                sender
                    .send(PivotizeMessage::NextPivotFound(next_pivot))
                    .expect("Could not send next pivot.");
            }
            PivotizeMessage::EliminateAllRows(pivot, pivot_row) => {
                mat.pivotize_with_row(pivot, pivot_row);
                sender
                    .send(PivotizeMessage::Completed)
                    .expect("Could not send elimination confirmation.");
            }
            PivotizeMessage::EliminateAllRowsBelow(pivot, pivot_row) => {
                mat.eliminate_rows_below_no_swap(pivot, &pivot_row);
                sender
                    .send(PivotizeMessage::Completed)
                    .expect("Could not send confirmation");
            }
            PivotizeMessage::GetRow(row_ix) => {
                let row = mat.get_row(row_ix);
                if row.is_zero() {
                    sender
                        .send(PivotizeMessage::RowNotFound)
                        .expect("Could not send back to coordinator.");
                } else {
                    sender
                        .send(PivotizeMessage::RowRetrieved(row))
                        .expect("Could not send back to coordinator.");
                }
            }
            PivotizeMessage::MakePivotRow(row_ix) => {
                if mat.create_pivot_in_row(row_ix) {
                    sender
                        .send(PivotizeMessage::RowRetrieved(mat.get_row(row_ix)))
                        .expect("Could not send back to coordinator.");
                } else {
                    sender
                        .send(PivotizeMessage::RowNotFound)
                        .expect("Could not send back to coordinator.")
                }
            }
            PivotizeMessage::SwapRows(row_1, row_2) => {
                mat.swap_rows(row_1, row_2);
            }
            PivotizeMessage::Quit => {
                break;
            }
            PivotizeMessage::NNZ(None) => {
                sender
                    .send(PivotizeMessage::NNZ(Some(mat.nnz())))
                    .expect("Could not send back to coordinator");
            }
            PivotizeMessage::Cache(path_buf) => {
                if path_buf.is_dir() {
                    let mut e = path_buf.clone();
                    let mut s = String::from("parallel_cache_");
                    s.push_str(&id.to_string()[..]);
                    e.push(&s[..]);
                    let data = &serde_json::to_string(&mat).unwrap();
                    let mut f = std::fs::File::create(e.as_path()).expect("cannot create file");
                    f.write(data.as_bytes()).expect("Cannot write data");
                    sender
                        .send(PivotizeMessage::Completed)
                        .expect("Cannot send to master");
                } else {
                    sender
                        .send(PivotizeMessage::Error)
                        .expect("path_buf was not a directory.");
                }
            }
            PivotizeMessage::VerfifyUpperTriangular(None) => {
                let is_subsection_upper_triangular = Some(mat.verify_upper_triangular());
                sender
                    .send(PivotizeMessage::VerfifyUpperTriangular(
                        is_subsection_upper_triangular,
                    ))
                    .expect("Cannot send result back to coordinator.");
            }
            PivotizeMessage::NumRows(None) => {
                sender
                    .send(PivotizeMessage::NumRows(Some(mat.ix_to_section.len())))
                    .expect("Cannot send result back to coordinator.");
            }
            PivotizeMessage::AddEntries(entries) => {
                mat.insert_entries(entries);
                sender
                    .send(PivotizeMessage::Completed)
                    .expect("Cannot send result back to coordinator.");
            }
        }
    }
    return mat;
}

#[derive(Debug, Clone, PartialEq)]
enum PivotizeMessage {
    FindNextPivot(Option<(usize, usize)>),
    NextPivotFound(Option<(usize, usize)>),
    EliminateAllRows((usize, usize), SparseVector),
    /// Eliminate all rows if the row index is greater than the row of the pivot
    EliminateAllRowsBelow((usize, usize), SparseVector),
    GetRow(usize),
    /// expects a `RowRetrieved` message in return with the new pivot row.
    MakePivotRow(usize),
    SwapRows(usize, usize),
    RowRetrieved(SparseVector),
    RowNotFound,
    /// None indicates a request, Some(n) is a response
    NNZ(Option<usize>),
    Cache(PathBuf),
    Completed,
    Quit,
    Error,
    VerfifyUpperTriangular(Option<bool>),
    /// None indicates a request, Some(n) a response.
    NumRows(Option<usize>),
    AddEntries(Vec<(usize, usize, FFRep)>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PivotCache {
    pivots: Vec<(usize, usize)>,
    num_rows: usize,
}

#[derive(Debug)]
pub struct ParallelFFMatrix {
    thread_handles: Vec<JoinHandle<SparseFFMatrix>>,
    channels: Vec<(Sender<PivotizeMessage>, Receiver<PivotizeMessage>)>,
    max_row_ix: usize,
    pivots: Vec<(usize, usize)>,
    field_mod: FFRep,
}

impl ParallelFFMatrix {
    pub fn new(matrices: Vec<SparseFFMatrix>) -> Self {
        let mut max_row_ix: usize = 0;
        let mut thread_handles = Vec::with_capacity(matrices.len());
        let mut channels = Vec::with_capacity(matrices.len());
        let mut field_mod = None;
        let mut id = 0;
        for mat in matrices {
            if field_mod.is_none() {
                field_mod = Some(mat.field_mod);
            }
            max_row_ix = max_row_ix.max(*mat.ix_to_section.keys().max().unwrap_or(&0));
            let (t1, r1) = mpsc::channel::<PivotizeMessage>();
            let (t2, r2) = mpsc::channel::<PivotizeMessage>();
            // let handle = thread::spawn(move || worker_thread_matrix_loop(mat, t1, r2, id.clone()));
            let stack_size = if mat.n_cols > 10_000_000 {
                64 * 1024
            } else {
                32 * 1024
            };
            let thread_builder = thread::Builder::new()
                .name(format!("worker {:}", id))
                .stack_size(stack_size);
            let handle = thread_builder
                .spawn(move || worker_thread_matrix_loop(mat, t1, r2, id.clone()))
                .expect("Cannot create new threads.");
            id += 1;
            // t2.send(PivotizeMessage::Test).unwrap();
            thread_handles.push(handle);
            channels.push((t2.clone(), r1));
        }
        Self {
            thread_handles,
            channels,
            max_row_ix,
            pivots: Vec::new(),
            field_mod: field_mod.unwrap(),
        }
    }

    pub fn num_threads(&self) -> usize {
        self.channels.len()
    }

    pub fn add_entries(&mut self, entries: Vec<(usize, usize, FFRep)>) {
        let mut hm: HashMap<usize, Vec<(usize, FFRep)>> = HashMap::new();
        for (row_ix, col_ix, entry) in entries {
            let e = hm.entry(row_ix).or_default();
            e.push((col_ix, entry));
        }
        let mut num_rows = HashMap::new();
        for (tx, _rx) in self.channels.iter() {
            tx.send(PivotizeMessage::NumRows(None)).unwrap();
        }
        for ix in 0..self.channels.len() {
            let rec = self.channels[ix].1.recv().unwrap();
            match rec {
                PivotizeMessage::NumRows(Some(n_rows)) => {
                    num_rows.insert(ix, n_rows);
                }
                _ => {
                    panic!("Incorrect response received.")
                }
            }
        }
        for (row_ix, entries) in hm.into_iter() {
            let sender_ix = *num_rows.iter().min_by_key(|(_k, v)| **v).unwrap().0;
            let new_entries = entries.into_iter().map(|(c, e)| (row_ix, c, e)).collect();
            self.channels[sender_ix]
                .0
                .send(PivotizeMessage::AddEntries(new_entries))
                .unwrap();
            if self.channels[sender_ix].1.recv().unwrap() != PivotizeMessage::Completed {
                panic!("Incorrect response received.")
            }
            let e = num_rows.get_mut(&sender_ix).unwrap();
            *e += 1;
            self.max_row_ix = self.max_row_ix.max(row_ix);
        }
    }

    pub fn from_disk(cache_dir: PathBuf, num_threads: usize) -> Option<Self> {
        if cache_dir.is_dir() {
            println!("cache_dir is a dir: {:?}", cache_dir);
            let mut mats: Vec<SparseFFMatrix> = cache_dir
                .read_dir()
                .expect("Directory? ")
                .into_iter()
                .filter_map(|x| x.ok())
                .filter_map(|path| {
                    let p = path.path();
                    let filename = p.file_name();
                    if filename.is_none() {
                        return None;
                    }
                    let filename = String::from(filename.unwrap().to_str().unwrap());
                    if filename.starts_with("parallel_cache_") {
                        let s = std::fs::read_to_string(path.path()).unwrap();
                        serde_json::from_str::<SparseFFMatrix>(&s[..]).ok()
                    } else {
                        None
                    }
                })
                .collect();
            if mats.len() == 0 {
                return None;
            }
            // let num_threads = usize::from(available_parallelism().unwrap());
            if num_threads != mats.len() {
                let mut new_mat =
                    SparseFFMatrix::new(0, 0, mats[0].field_mod, MemoryLayout::RowMajor);
                for mut mat in mats.drain(..) {
                    new_mat.append(&mut mat);
                }
                if mats.len() != 0 {
                    panic!("Attempting to realign matrices with new thread_count but could not drain matrices.")
                }
                let mut row_batches = Vec::with_capacity(num_threads);
                for _ in 0..num_threads {
                    row_batches.push(HashSet::with_capacity(new_mat.n_rows / num_threads));
                }
                let mut batch_ix = 0;
                for row_ix in new_mat.ix_to_section.keys() {
                    let batch = row_batches.get_mut(batch_ix).unwrap();
                    batch.insert(*row_ix);
                    batch_ix += 1;
                    batch_ix %= row_batches.len();
                }
                mats = row_batches
                    .into_iter()
                    .map(|ix_batch: HashSet<usize>| new_mat.split(ix_batch))
                    .collect();
                if mats.len() != num_threads {
                    panic!("Readjusting matrices failed, number of matrices does not match number of rows.")
                }
                for entry in cache_dir.read_dir().unwrap() {
                    if let Ok(e) = entry {
                        let path = e.path();
                        if path.starts_with("parallel_cache_") {
                            std::fs::remove_file(path).expect("Could not remove old cache.");
                        }
                    }
                }
                for ix in 0..mats.len() {
                    let mut e = cache_dir.clone();
                    let mut s = String::from("parallel_cache_");
                    s.push_str(&ix.to_string()[..]);
                    e.push(&s[..]);
                    let data = &serde_json::to_string(&mats[ix]).unwrap();
                    let mut f = std::fs::File::create(e.as_path()).expect("cannot create file");
                    f.write(data.as_bytes()).expect("Cannot write data");
                }
            }
            let mut thread_handles = Vec::with_capacity(mats.len());
            let mut channels = Vec::with_capacity(mats.len());
            let mut field_mod = None;
            let mut id = 0;
            let mut max_row_ix = 0_usize;
            for mat in mats {
                if field_mod.is_none() {
                    field_mod = Some(mat.field_mod);
                }
                max_row_ix = max_row_ix.max(*mat.ix_to_section.keys().max().unwrap());
                let (t1, r1) = mpsc::channel::<PivotizeMessage>();
                let (t2, r2) = mpsc::channel::<PivotizeMessage>();
                let stack_size = if mat.n_cols > 10_000_000 {
                    64 * 1024
                } else {
                    32 * 1024
                };
                let thread_builder = thread::Builder::new()
                    .name(format!("worker {:}", id))
                    .stack_size(stack_size);
                let handle = thread_builder
                    .spawn(move || worker_thread_matrix_loop(mat, t1, r2, id.clone()))
                    .expect("Cannot create new threads.");
                id += 1;
                thread_handles.push(handle);
                channels.push((t2.clone(), r1));
            }
            let mut pivot_cache_path = cache_dir.clone();
            pivot_cache_path.push("pivot_cache_parallel");
            let pivot_cache_file = std::fs::File::open(pivot_cache_path.as_path()).unwrap();
            let cached_data: PivotCache = serde_json::from_reader(pivot_cache_file).unwrap();
            return Some(Self {
                thread_handles,
                channels,
                pivots: cached_data.pivots,
                max_row_ix,
                field_mod: field_mod.unwrap(),
            });
        }
        log::error!("Need to pass in directory to read ParallelFFMatrix from cache.");
        None
    }

    pub fn verify_is_upper_triangular(&self) -> bool {
        for (tx, _rx) in self.channels.iter() {
            tx.send(PivotizeMessage::VerfifyUpperTriangular(None))
                .expect("Cannot send verify message.");
        }
        let mut all_workers_upper_triangular = true;
        for (_, rx) in self.channels.iter() {
            if let Ok(rec) = rx.recv() {
                match rec {
                    PivotizeMessage::VerfifyUpperTriangular(Some(section_is_upper)) => {
                        all_workers_upper_triangular &= section_is_upper;
                    }
                    _ => {
                        panic!("Message received out of order.");
                    }
                }
            }
        }
        all_workers_upper_triangular
    }

    pub fn cache(&self, dir: &Path) {
        for (tx, _rx) in self.channels.iter() {
            tx.send(PivotizeMessage::Cache(PathBuf::from(dir)))
                .expect("Cannot send to channel.");
        }
        for (_tx, rx) in self.channels.iter() {
            let ret = rx.recv().unwrap();
            if ret == PivotizeMessage::Error {
                panic!("Worker failed to cache, dont know what to do.")
            }
        }
        let mut pivot_cache_file = PathBuf::from(dir);
        pivot_cache_file.push("pivot_cache_parallel");
        let pivot_cache = PivotCache {
            pivots: self.pivots.clone(),
            num_rows: self.max_row_ix,
        };
        let s = serde_json::to_string(&pivot_cache).unwrap();
        std::fs::write(pivot_cache_file.as_path(), &s[..]).unwrap();
    }

    /// Number NonZeros
    pub fn nnz(&self) -> usize {
        for (tx, _) in self.channels.iter() {
            tx.send(PivotizeMessage::NNZ(None))
                .expect("Could not send.");
        }
        let mut stats = Vec::new();
        for (_, rx) in self.channels.iter() {
            if let PivotizeMessage::NNZ(Some(x)) = rx.recv().expect("Could not receive nnz reply.")
            {
                stats.push(x);
            }
        }
        // println!(
        //     "min, avg, max: {:?}, {:}, {:?}",
        //     stats.iter().min(),
        //     stats.iter().sum::<usize>() as f64 / stats.len() as f64,
        //     stats.iter().max()
        // );
        stats.iter().fold(0, |mut acc, x| {
            acc += x;
            acc
        })
    }
    fn get_row(&self, row_ix: usize) -> SparseVector {
        for (tx, _) in self.channels.iter() {
            tx.send(PivotizeMessage::GetRow(row_ix))
                .expect("Could not send.");
        }
        let mut ret = SparseVector::new_empty();
        for (_, rx) in self.channels.iter() {
            let r = rx.recv().expect("Could not receive row.");
            match r {
                PivotizeMessage::RowRetrieved(row) => {
                    ret = row;
                }
                PivotizeMessage::RowNotFound => continue,
                _ => log::error!("Why is channel responding without the row?"),
            }
        }
        ret
    }

    fn swap_rows(&mut self, row_1: usize, row_2: usize) {
        if row_1 == row_2 {
            return;
        }
        for (tx, _) in self.channels.iter() {
            tx.send(PivotizeMessage::SwapRows(row_1, row_2))
                .expect("Could not send swap message.");
        }
    }

    pub fn pivotize_row(&mut self, pivot_row_ix: usize) -> Option<usize> {
        let mut pivot_row = self.get_row(pivot_row_ix);
        let col_ix = pivot_row.first_nonzero();
        if col_ix.is_none() {
            return None;
        }
        let pivot_col_ix = col_ix.unwrap().0;
        let entry = FF::new(col_ix.unwrap().1, self.field_mod);
        pivot_row.scale(entry.modular_inverse().0, self.field_mod);
        for (tx, _) in self.channels.iter() {
            tx.send(PivotizeMessage::EliminateAllRows(
                (pivot_row_ix, pivot_col_ix),
                pivot_row.clone(),
            ))
            .expect("Could not send elimination.");
        }
        for (_, rx) in self.channels.iter() {
            let rec = rx
                .recv()
                .expect("Could not receive elimination confirmation.");
            match rec {
                PivotizeMessage::Completed => continue,
                _ => {
                    log::error!("Received the following message in pivotize row: {:?}", rec);
                    panic!("Row could not be eliminated?")
                }
            }
        }
        Some(pivot_col_ix)
    }

    pub fn quit(self) -> SparseFFMatrix {
        for (tx, _) in self.channels {
            tx.send(PivotizeMessage::Quit)
                .expect("Could not quit... hustle too hard.");
        }
        self.thread_handles
            .into_iter()
            .map(|handle| handle.join().expect("Some thread broke"))
            .fold(
                SparseFFMatrix::new(0, 0, self.field_mod, MemoryLayout::RowMajor),
                |mut acc, mut mat| {
                    acc.append(&mut mat);
                    acc
                },
            )
    }

    fn find_next_pivot(&mut self, prev_pivot: Option<(usize, usize)>) -> Option<(usize, usize)> {
        let mut worker_pivots = Vec::with_capacity(self.channels.len());
        for (tx, _) in self.channels.iter() {
            tx.send(PivotizeMessage::FindNextPivot(prev_pivot))
                .expect("Could not send message")
        }
        for (_, rx) in self.channels.iter() {
            match rx.recv().expect("Could not receive found pivot.") {
                PivotizeMessage::NextPivotFound(piv) => {
                    worker_pivots.push(piv);
                }
                error => {
                    log::error!(
                        "Finding pivots routine broken. Recieved this instead: {:?}",
                        error
                    );
                    panic!()
                }
            }
        }
        worker_pivots.sort_by(|a, b| match (a, b) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Greater,
            (Some(_), None) => Ordering::Less,
            (Some((row_a, col_a)), Some((row_b, col_b))) => {
                if col_a == col_b {
                    row_a.cmp(row_b)
                } else {
                    col_a.cmp(col_b)
                }
            }
        });
        if worker_pivots.len() == 0 {
            return None;
        }
        worker_pivots[0]
    }

    fn ensure_pivot_with_swap(
        &mut self,
        previous_pivot: Option<(usize, usize)>,
    ) -> Option<((usize, usize), SparseVector)> {
        let new_pivot = self.find_next_pivot(previous_pivot);
        if new_pivot.is_none() {
            return None;
        }
        let mut new_pivot = new_pivot.unwrap();
        if previous_pivot.is_none() {
            self.swap_rows(0, new_pivot.0);
            new_pivot.0 = 0;
        } else {
            self.swap_rows(previous_pivot.unwrap().0 + 1, new_pivot.0);
            new_pivot.0 = previous_pivot.unwrap().0 + 1;
        }

        // // make sure the pivot entry is 1
        for (tx, _) in self.channels.iter() {
            tx.send(PivotizeMessage::MakePivotRow(new_pivot.0))
                .expect("Could not send message.");
        }
        let mut pivot_row = SparseVector::new_empty();
        for (_, rx) in self.channels.iter() {
            let rec = rx.recv().expect("Could not retrieve message");
            match rec {
                PivotizeMessage::RowRetrieved(row) => {
                    pivot_row = row;
                }
                PivotizeMessage::RowNotFound => {
                    continue;
                }
                _ => {
                    log::error!(
                        "Worker threads out of sync. Expected row retrieval but got {:?}",
                        rec
                    );
                    panic!()
                }
            }
        }
        if pivot_row.is_zero() {
            panic!("Have zero pivot row.")
        }
        Some((new_pivot, pivot_row))
    }

    /// Returns the pivots created during this invocation. Guarantees the
    /// resulting matrix will be in reduced row echelon form
    pub fn rref(
        &mut self,
        cache_dir: Option<&Path>,
        cache_rate: Option<usize>,
        log_rate: Option<usize>,
    ) -> Vec<(usize, usize)> {
        trace!("rref!");
        if self.max_row_ix == 0 {
            return Vec::new();
        }
        let mut current_pivot = self.ensure_pivot_with_swap(self.pivots.last().cloned());
        let mut ret = Vec::new();
        let mut cur_row_ix = if let Some((row_ix, _)) = self.pivots.last() {
            *row_ix
        } else {
            0
        };
        // let log_rate = (self.num_rows / 1000).max(1);
        let mut num_pivots_made = 0;
        let mut time_spent_pivoting = 0.0;

        let mut cache_checkpoints = Vec::new();
        if let Some(cr) = cache_rate {
            let prev_cache_point = cur_row_ix / cr;
            let mut next_cache_point = prev_cache_point + cr;
            while next_cache_point < self.max_row_ix - 1 {
                cache_checkpoints.push(next_cache_point);
                next_cache_point += cr;
            }
        }
        log::trace!(
            "Starting Row Echelon computation. Counter: {:}, cache_rate: {:?}, Cache checkpoints: {:?}, num_rows: {:}, log_rate: {:?}",
            cur_row_ix, cache_rate, cache_checkpoints,self.max_row_ix, log_rate
        );
        while current_pivot.is_some() {
            let elimination_start = Instant::now();
            for (tx, _) in self.channels.iter() {
                let (pivot, row) = current_pivot.clone().unwrap();
                tx.send(PivotizeMessage::EliminateAllRows(pivot, row))
                    .expect("Could not send message.");
            }
            for (_, rx) in self.channels.iter() {
                let rec = rx.recv().expect("Could not retrieve message.");
                if rec != PivotizeMessage::Completed {
                    log::error!("Worker out of sync. Bailing.");
                    panic!()
                }
            }

            self.pivots.push(current_pivot.clone().unwrap().0);
            ret.push(current_pivot.clone().unwrap().0);
            current_pivot = self.ensure_pivot_with_swap(current_pivot.map(|(piv, _)| piv));
            time_spent_pivoting += elimination_start.elapsed().as_secs_f64();
            num_pivots_made += 1;
            cur_row_ix += 1;
            if Some(&cur_row_ix) == cache_checkpoints.first() && cache_dir.is_some() {
                cache_checkpoints.remove(0);
                log::trace!("Parallel Matrix is caching.");
                self.cache(cache_dir.unwrap());
                log::trace!(
                    "Done Caching! Next checkpoint: {:?}",
                    cache_checkpoints.first()
                );
            }
            if let Some(lr) = log_rate {
                if cur_row_ix % lr == 0 {
                    let time_per_pivot = time_spent_pivoting / num_pivots_made as f64;
                    let worst_time_remaining =
                        (self.max_row_ix - cur_row_ix) as f64 * time_per_pivot;
                    log::trace!(
                        "At row: {:} out of {:}. Found {:} pivots and matrix has {:} nonzeros",
                        cur_row_ix,
                        self.max_row_ix,
                        self.pivots.len(),
                        self.nnz()
                    );
                    num_pivots_made = 0;
                    time_spent_pivoting = 0.0;
                    log::trace!(
                        "Estimated time remaining: {:} days, time per pivot: {:}",
                        worst_time_remaining / (3600.0 * 24.0),
                        time_per_pivot
                    );
                }
            }
        }
        if let Some(cd) = cache_dir {
            self.cache(cd);
        }
        ret
    }

    /// Returns the pivots in the row echelon form matrix. Will swap rows to put pivots
    /// as high up in the matrix as possible.
    pub fn row_echelon_form(
        &mut self,
        cache_dir: Option<&Path>,
        cache_rate: Option<usize>,
        log_rate: Option<usize>,
    ) -> Vec<(usize, usize)> {
        // let mut pivots: Vec<(usize, usize)> = Vec::new();
        if self.max_row_ix == 0 {
            return Vec::new();
        }
        println!("In row echelon form");
        let mut current_pivot = self.ensure_pivot_with_swap(self.pivots.last().cloned());
        println!("ensured pivot, nrows x ncols = {:}", self.max_row_ix);
        let mut cur_row_ix = if let Some((row_ix, _)) = self.pivots.last() {
            *row_ix
        } else {
            0
        };
        // let log_rate = (self.num_rows / 1000).max(1);
        let mut num_pivots_made = 0;
        let mut time_spent_pivoting = 0.0;

        let mut cache_checkpoints = Vec::new();
        if let Some(cr) = cache_rate {
            let prev_cache_point = cur_row_ix / cr;
            let mut next_cache_point = prev_cache_point + cr;
            while next_cache_point < self.max_row_ix - 1 {
                cache_checkpoints.push(next_cache_point);
                next_cache_point += cr;
            }
        }
        log::trace!(
            "Starting Row Echelon computation. Counter: {:}, cache_rate: {:?}, Cache checkpoints: {:?}, num_rows: {:}, log_rate: {:?}",
            cur_row_ix, cache_rate, cache_checkpoints,self.max_row_ix, log_rate
        );
        while current_pivot.is_some() {
            let elimination_start = Instant::now();
            for (tx, _) in self.channels.iter() {
                let (pivot, row) = current_pivot.clone().unwrap();
                tx.send(PivotizeMessage::EliminateAllRowsBelow(pivot, row))
                    .expect("Could not send message.");
            }
            for (_, rx) in self.channels.iter() {
                let rec = rx.recv().expect("Could not retrieve message.");
                if rec != PivotizeMessage::Completed {
                    log::error!("Worker out of sync. Bailing.");
                    panic!()
                }
            }

            self.pivots.push(current_pivot.clone().unwrap().0);
            current_pivot = self.ensure_pivot_with_swap(current_pivot.map(|(piv, _)| piv));
            time_spent_pivoting += elimination_start.elapsed().as_secs_f64();
            num_pivots_made += 1;
            cur_row_ix += 1;
            if Some(&cur_row_ix) == cache_checkpoints.first() && cache_dir.is_some() {
                cache_checkpoints.remove(0);
                log::trace!("Parallel Matrix is caching.");
                self.cache(cache_dir.unwrap());
                log::trace!(
                    "Done Caching! Next checkpoint: {:?}",
                    cache_checkpoints.first()
                );
            }
            if let Some(lr) = log_rate {
                if cur_row_ix % lr == 0 {
                    let time_per_pivot = time_spent_pivoting / num_pivots_made as f64;
                    let worst_time_remaining =
                        (self.max_row_ix - cur_row_ix) as f64 * time_per_pivot;
                    log::trace!(
                        "At row: {:} out of {:}. Found {:} pivots and matrix has {:} nonzeros",
                        cur_row_ix,
                        self.max_row_ix,
                        self.pivots.len(),
                        self.nnz()
                    );
                    num_pivots_made = 0;
                    time_spent_pivoting = 0.0;
                    log::trace!(
                        "Estimated time remaining: {:} days, time per pivot: {:}",
                        worst_time_remaining / (3600.0 * 24.0),
                        time_per_pivot
                    );
                }
            }
        }
        if let Some(cd) = cache_dir {
            self.cache(cd);
        }
        self.pivots.clone()
    }
}

#[cfg(test)]
mod test {
    use std::{collections::HashSet, path::PathBuf};

    use simple_logger::SimpleLogger;

    use crate::matrices::{
        mat_trait::RankMatrix,
        parallel_matrix::ParallelFFMatrix,
        sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    };

    #[test]
    fn parallelization() {
        let mut entries = Vec::new();
        for ix in 0..10 {
            entries.push((0 as usize, ix as usize, 1_u32));
        }
        for ix in 1..100 {
            entries.push((ix, 0, 1));
            entries.push((ix, 1 + ix % 7, ix as u32))
        }
        let mut mat = <SparseFFMatrix as RankMatrix>::new(199);
        mat.insert_entries(entries);

        println!("mat before:\n {:}", mat.clone().to_dense());
        // mat.parallel_eliminate_all_rows((0, 0));
        println!("mat after:\n {:}", mat.clone().to_dense());
    }

    #[test]
    fn add_entries() {
        let mut pm = ParallelFFMatrix::new(vec![
            SparseFFMatrix::new_with_entries(1, 2, 7, MemoryLayout::RowMajor, vec![(0, 0, 1)]),
            SparseFFMatrix::new_with_entries(2, 2, 7, MemoryLayout::RowMajor, vec![(1, 1, 2)]),
        ]);
        pm.add_entries(vec![(2, 2, 3), (3, 3, 4)]);
        let s = pm.quit();
        let m = s.to_dense();
        println!("{:}", m);
    }

    #[test]
    fn upper_triangular() {
        let mut mat = SparseFFMatrix::new_random(40, 40, 7, 0.05);
        let mut parallel_mat = mat.clone().split_into_parallel(
            mat.ix_to_section.keys().cloned().collect(),
            std::thread::available_parallelism().unwrap().into(),
        );
        let parallel_pivots = parallel_mat.row_echelon_form(None, None, None);
        let rayon_pivots = mat.row_echelon_form();
        assert_eq!(parallel_pivots, rayon_pivots);
        assert!(parallel_mat.verify_is_upper_triangular());
        let parallel_mat_rref = parallel_mat.quit();
        for ix in 0..parallel_mat_rref.n_rows {
            let parallel_row = parallel_mat_rref.get_row(ix);
            let other_row = mat.get_row(ix);
            println!("row ix: {ix}\nl: {:?}\nr: {:?}", parallel_row, other_row);
        }
        assert_eq!(parallel_mat_rref, mat);
    }
    #[test]
    fn cache_validation() {
        let _logger = SimpleLogger::new().init().unwrap();
        let dir = PathBuf::from("/Users/matt/repos/qec/tmp/par_test");
        for entry in dir.read_dir().expect("Directory? ") {
            let ent = entry.unwrap();
            std::fs::remove_file(ent.path()).expect("Could not remove file.");
        }
        let first_test = ParallelFFMatrix::from_disk(dir.clone(), 2);
        assert!(first_test.is_none());
        let dim = 100;
        let mut mat = SparseFFMatrix::new_random(dim, 2 * dim, 3, 1.0);
        let mut par = mat.split_into_parallel((0..dim).collect(), 2);
        par.row_echelon_form(Some(dir.as_path()), None, None);
        par.quit();
        let from_cache = ParallelFFMatrix::from_disk(dir.clone(), 2);
        assert!(from_cache.is_some());

        dbg!(&from_cache);
        if let Some(par_mat) = from_cache {
            par_mat.quit();
        }
    }

    #[test]
    fn stack_overflow_replication() {
        let mut dense = SparseFFMatrix::new_random(1, 1000, 3, 1.0);
        println!("{:}", dense.clone().to_dense());
        let mut par = dense.split_into_parallel(HashSet::from([0]), 1);
        let _pivots = par.row_echelon_form(None, None, None);
    }
}
