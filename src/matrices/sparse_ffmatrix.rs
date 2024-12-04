use core::{panic, time};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, HashMap, HashSet},
    fs::{self, File},
    io::{Read, Write},
    path::{Path, PathBuf},
    rc::Rc,
    str::FromStr,
    sync::mpsc::{Receiver, Sender},
    thread::{self, available_parallelism, JoinHandle},
    time::Instant,
};

use fxhash::FxHashSet;
use rand::{thread_rng, Rng, SeedableRng};
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use super::{ffmatrix::FFMatrix, sparse_vec::SparseVector};
use crate::math::finite_field::{FFRep, FiniteField as FF};
use crate::matrices::RankMatrix;
use std::sync::mpsc;
/// Used to determine if a memory layout for a matrix
/// is RowMajor or ColMajor, used to determine if the sparse sections
/// retrieved are rows or columns.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MemoryLayout {
    RowMajor,
    ColMajor,
}

pub fn benchmark_rate(dim: usize, num_samples: usize) {
    log::trace!("Generating benchmark matrices.");
    let mut rayon_mats: Vec<SparseFFMatrix> = (0..num_samples)
        .into_iter()
        .map(|_| SparseFFMatrix::new_random(dim, 2 * dim, 3, 1e-7))
        .collect();
    let mut home_roll_mats: Vec<ParallelFFMatrix> = rayon_mats
        .clone()
        .into_iter()
        .map(|mut mat| {
            mat.split_into_parallel((0..100).collect(), available_parallelism().unwrap().into())
        })
        .collect();
    log::trace!("Starting Rayon timing.");
    let rayon_start = Instant::now();
    for mat in rayon_mats.iter_mut() {
        mat.row_echelon_form();
    }
    let rayon_avg_time = rayon_start.elapsed().as_secs_f64() / rayon_mats.len() as f64;
    log::trace!("Starting Home Rolled timing.");
    let home_roll_start = Instant::now();
    for mat in home_roll_mats.iter_mut() {
        mat.row_echelon_form();
    }
    let home_roll_avg_time = home_roll_start.elapsed().as_secs_f64() / home_roll_mats.len() as f64;
    for mat in home_roll_mats {
        mat.quit();
    }
    println!("rayon = {:}", rayon_avg_time);
    println!("home_roll = {:}", home_roll_avg_time);
    println!(
        "rayon / home_roll = {:}",
        rayon_avg_time / home_roll_avg_time
    );
}

impl MemoryLayout {
    pub fn get_row_col(&self, section: usize, position: usize) -> (usize, usize) {
        match self {
            MemoryLayout::RowMajor => (section, position),
            MemoryLayout::ColMajor => (position, section),
        }
    }

    /// Gets the `section` (which memory section to retrieve) and `position` (which position within the memory section to get)
    ///
    /// Returns `(section, position)`.
    pub fn get_section_position(&self, row_ix: usize, col_ix: usize) -> (usize, usize) {
        match self {
            MemoryLayout::RowMajor => (row_ix, col_ix),
            MemoryLayout::ColMajor => (col_ix, row_ix),
        }
    }
}

fn worker_thread_matrix_loop(
    mut mat: SparseFFMatrix,
    sender: Sender<PivotizeMessage>,
    receiver: Receiver<PivotizeMessage>,
) -> SparseFFMatrix {
    let mut err_count = 0;
    loop {
        let message = receiver.recv();
        if message.is_err() {
            log::error!("Could not retrieve message?");
            err_count += 1;
            if err_count > 10 {
                panic!("Too many receiving errors.");
            }
            continue;
        }
        let message = message.unwrap();
        match message {
            PivotizeMessage::RowEliminated
            | PivotizeMessage::RowNotFound
            | PivotizeMessage::RowRetrieved(_)
            | PivotizeMessage::NNZ(Some(_))
            | PivotizeMessage::Capacity(Some(_))
            | PivotizeMessage::NextPivotFound(_)
            | PivotizeMessage::Completed => {
                log::trace!("Coordinator should not be sending these messages.");
            }
            PivotizeMessage::Test => log::trace!("Test received!"),
            PivotizeMessage::FindNextPivot(previous_pivot) => {
                let next_pivot = mat.find_next_pivot(previous_pivot);
                sender
                    .send(PivotizeMessage::NextPivotFound(next_pivot))
                    .expect("Could not send next pivot.");
            }
            PivotizeMessage::EliminateAllRows(pivot, pivot_row) => {
                mat.pivotize_with_row(pivot, pivot_row);
                sender
                    .send(PivotizeMessage::RowEliminated)
                    .expect("Could not send elimination confirmation.");
            }
            PivotizeMessage::EliminateAllRowsBelow(pivot, pivot_row) => {
                mat.eliminate_rows_below_no_swap(pivot, &pivot_row);
                sender
                    .send(PivotizeMessage::Completed)
                    .expect("Could not send confirmation");
            }
            PivotizeMessage::GetRow(row_ix) => {
                let row = mat.row(row_ix);
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
                        .send(PivotizeMessage::RowRetrieved(mat.row(row_ix)))
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
            PivotizeMessage::Capacity(None) => {
                sender
                    .send(PivotizeMessage::NNZ(Some(mat.capacity())))
                    .expect("Could not send back to coordinator");
            }
        }
    }
    return mat;
}

pub struct ParallelFFMatrix {
    thread_handles: Vec<JoinHandle<SparseFFMatrix>>,
    channels: Vec<(Sender<PivotizeMessage>, Receiver<PivotizeMessage>)>,
    field_mod: FFRep,
}

#[derive(Debug, Clone, PartialEq)]
enum PivotizeMessage {
    /// Eliminate all rows if the row index is greater than the row of the pivot
    EliminateAllRowsBelow((usize, usize), SparseVector),
    FindNextPivot(Option<(usize, usize)>),
    NextPivotFound(Option<(usize, usize)>),
    EliminateAllRows((usize, usize), SparseVector),
    RowEliminated,
    GetRow(usize),
    /// expects a `RowRetrieved` message in return with the new pivot row.
    MakePivotRow(usize),
    SwapRows(usize, usize),
    RowRetrieved(SparseVector),
    RowNotFound,
    /// None indicates a request, Some(n) is a response
    NNZ(Option<usize>),
    Capacity(Option<usize>),
    Completed,
    Quit,
    Test,
}

impl ParallelFFMatrix {
    pub fn new(matrices: Vec<SparseFFMatrix>) -> Self {
        // log::info!(
        //     "Number of available threads: {:}",
        //     std::thread::available_parallelism().unwrap()
        // );
        let mut thread_handles = Vec::with_capacity(matrices.len());
        let mut channels = Vec::with_capacity(matrices.len());
        let mut field_mod = None;
        for mat in matrices {
            if field_mod.is_none() {
                field_mod = Some(mat.field_mod);
            }
            let (t1, r1) = mpsc::channel::<PivotizeMessage>();
            let (t2, r2) = mpsc::channel::<PivotizeMessage>();
            let handle = thread::spawn(|| worker_thread_matrix_loop(mat, t1, r2));
            // t2.send(PivotizeMessage::Test).unwrap();
            thread_handles.push(handle);
            channels.push((t2.clone(), r1));
        }
        Self {
            thread_handles,
            channels,
            field_mod: field_mod.unwrap(),
        }
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
        println!(
            "min, avg, max: {:?}, {:}, {:?}",
            stats.iter().min(),
            stats.iter().sum::<usize>() as f64 / stats.len() as f64,
            stats.iter().max()
        );
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
                PivotizeMessage::RowEliminated => continue,
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
    /// Returns the pivots in the row echelon form matrix. Will swap rows to put pivots
    /// as high up in the matrix as possible.
    pub fn row_echelon_form(&mut self) -> Vec<(usize, usize)> {
        let mut pivots: Vec<(usize, usize)> = Vec::new();
        let mut current_pivot = self.ensure_pivot_with_swap(None);
        while current_pivot.is_some() {
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
            pivots.push(current_pivot.clone().unwrap().0);
            current_pivot = self.ensure_pivot_with_swap(current_pivot.map(|(piv, _)| piv));
        }
        pivots
    }
}

// TODO: Should I implement a Compressed Sparse Row version of this?
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SparseFFMatrix {
    ix_to_section: BTreeMap<usize, SparseVector>,
    pub n_rows: usize,
    pub n_cols: usize,
    pub field_mod: FFRep,
    memory_layout: MemoryLayout,
}

impl SparseFFMatrix {
    /// Creates an all zeros sparse matrix with the specified memory layout.
    /// If your matrix is row-sparse then use `MemoryLayout::RowMajor`, if it is
    /// column sparse use `MemoryLayout::ColMajor`.
    pub fn new(n_rows: usize, n_cols: usize, field_mod: FFRep, layout: MemoryLayout) -> Self {
        Self {
            ix_to_section: BTreeMap::new(),
            n_rows,
            n_cols,
            field_mod,
            memory_layout: layout,
        }
    }

    /// Shuffles rows and columns to reduce gaps and start at 1.
    pub fn reorder_row_and_col_indices(&mut self) {
        if self.memory_layout != MemoryLayout::RowMajor {
            panic!("No clean up for Col Major matrices.")
        }
        let mut row_counter = 0;
        let mut col_counter = 0;
        let mut old_to_new_col: HashMap<usize, usize> = HashMap::new();
        let mut new_ix_to_section: BTreeMap<usize, SparseVector> = BTreeMap::new();
        while self.ix_to_section.is_empty() == false {
            let (_, row) = self.ix_to_section.pop_last().unwrap();
            let new_row_ix = row_counter;
            row_counter += 1;
            let entries = row
                .to_vec()
                .into_iter()
                .map(|(ix, entry)| {
                    let new_ix = old_to_new_col.entry(ix).or_insert_with(|| {
                        let ret = col_counter;
                        col_counter += 1;
                        ret
                    });
                    (*new_ix, entry)
                })
                .collect();
            new_ix_to_section.insert(new_row_ix, SparseVector::new_with_entries(entries));
        }
        self.ix_to_section.append(&mut new_ix_to_section);
    }

    /// Number NonZeros
    pub fn nnz(&self) -> usize {
        self.ix_to_section.iter().fold(0, |mut acc, (_, row)| {
            acc += row.nnz();
            acc
        })
    }

    pub fn capacity(&self) -> usize {
        self.ix_to_section.iter().fold(0, |mut acc, (_, row)| {
            acc += row.capacity();
            acc
        })
    }

    pub fn new_with_entries(
        n_rows: usize,
        n_cols: usize,
        field_mod: FFRep,
        layout: MemoryLayout,
        entries: Vec<(usize, usize, FFRep)>,
    ) -> Self {
        let mut new = SparseFFMatrix::new(n_rows, n_cols, field_mod, layout);
        new.insert_entries(entries);
        new
    }

    pub fn new_random(n_rows: usize, n_cols: usize, field_mod: FFRep, avg_sparsity: f64) -> Self {
        let mut entries = Vec::with_capacity(n_rows * n_cols);
        let mut rng = thread_rng();
        for row_ix in 0..n_rows {
            for col_ix in 0..n_cols {
                if rng.gen_bool(avg_sparsity) {
                    entries.push((row_ix, col_ix, rng.gen_range(0..field_mod)));
                }
            }
        }
        Self::new_with_entries(n_rows, n_cols, field_mod, MemoryLayout::RowMajor, entries)
    }

    pub fn append(&mut self, other: &mut Self) {
        for (ix, other_row) in other.ix_to_section.iter_mut() {
            let e = self
                .ix_to_section
                .entry(*ix)
                .or_insert(SparseVector::new_empty());
            e.add_to_self(&other_row, self.field_mod);
            self.n_rows = self.n_rows.max(*ix);
            self.n_cols = self.n_cols.max(other_row.max_index());
        }
    }

    pub fn row(&self, row_ix: usize) -> SparseVector {
        if self.memory_layout == MemoryLayout::RowMajor {
            if let Some(v) = self.ix_to_section.get(&row_ix) {
                v.clone()
            } else {
                SparseVector::new_empty()
            }
        } else {
            panic!("Currently cannot get rows for ColMajor matrices.")
        }
    }

    pub fn normalize_pivot(&mut self, pivot: (usize, usize)) {
        if self.memory_layout == MemoryLayout::RowMajor {
            if let Some(v) = self.ix_to_section.get(&pivot.0) {
                let x = v.query(&pivot.1);
                if x == 0 {
                    return;
                }
                let mod_inv = FF::new(x, self.field_mod).modular_inverse();
                self.scale_section(pivot.0, mod_inv.0);
            }
        } else {
            panic!("Currently do not support pivoting for ColMajor matrices.")
        }
    }

    pub fn to_disk(&self, filename: &Path) {
        let s = serde_json::to_string(self).expect("Could not serialize self.");
        let mut file = File::create(filename).expect("Could not create file for write");
        file.write_all(s.as_bytes()).expect("Could not write");
    }

    pub fn from_disk(filename: &Path) -> Self {
        let mut s = String::new();
        let mut file = File::open(filename).expect("Could not open file for matrix read.");
        file.read_to_string(&mut s)
            .expect("Could not read matrix from disk");
        serde_json::from_str::<SparseFFMatrix>(&s).expect("Could not deserialize file.")
    }

    /// Assumes that the sections are appropriately ordered. Assumes that the
    /// largest seen index in any section is the appropriate number of rows/
    /// cols.
    pub fn new_with_sections(
        field_mod: FFRep,
        layout: MemoryLayout,
        sections: &Vec<SparseVector>,
    ) -> Self {
        let mut section_len = 0;
        let mut ix_to_section = BTreeMap::new();
        for ix in 0..sections.len() {
            section_len = section_len.max(sections[ix].max_index());
            ix_to_section.insert(ix, sections[ix].clone());
        }
        let (n_rows, n_cols) = layout.get_row_col(sections.len(), section_len);
        SparseFFMatrix {
            ix_to_section,
            n_rows,
            n_cols,
            field_mod,
            memory_layout: layout,
        }
    }

    pub fn transpose(&mut self) {
        match self.memory_layout {
            MemoryLayout::RowMajor => {
                self.memory_layout = MemoryLayout::ColMajor;
            }
            MemoryLayout::ColMajor => {
                self.memory_layout = MemoryLayout::RowMajor;
            }
        }
    }

    pub fn to_row_layout(&mut self) {
        if self.memory_layout == MemoryLayout::ColMajor {
            self.swap_layout();
        }
    }

    pub fn to_col_layout(&mut self) {
        if self.memory_layout == MemoryLayout::RowMajor {
            self.swap_layout();
        }
    }

    pub fn swap_layout(&mut self) {
        let mut ix_to_sections: BTreeMap<usize, SparseVector> = BTreeMap::new();
        for col_ix in 0..self.n_cols {
            for row_ix in 0..self.n_rows {
                let entry = self.get(row_ix, col_ix);
                match self.memory_layout {
                    MemoryLayout::RowMajor => {
                        let section = ix_to_sections
                            .entry(col_ix)
                            .or_insert(SparseVector::new_empty());
                        section.insert(row_ix, entry);
                    }
                    MemoryLayout::ColMajor => {
                        let section = ix_to_sections
                            .entry(row_ix)
                            .or_insert(SparseVector::new_empty());
                        section.insert(col_ix, entry);
                    }
                }
            }
        }
        self.ix_to_section = ix_to_sections;
        self.memory_layout = match self.memory_layout {
            MemoryLayout::RowMajor => MemoryLayout::ColMajor,
            MemoryLayout::ColMajor => MemoryLayout::RowMajor,
        };
    }

    pub fn is_square(&self) -> bool {
        self.n_rows == self.n_cols
    }

    /// Gets the row and column indices associated with the given `section` and `position`
    pub fn get_row_col(&self, section: usize, position: usize) -> (usize, usize) {
        match self.memory_layout {
            MemoryLayout::RowMajor => (section, position),
            MemoryLayout::ColMajor => (position, section),
        }
    }

    /// Adds a collection of (`row_ix`, `col_ix`, `entry`) to the matrix.
    pub fn insert_entries(&mut self, entries: Vec<(usize, usize, FFRep)>) {
        for (row_ix, col_ix, entry) in entries.into_iter() {
            self.insert(row_ix, col_ix, entry);
        }
    }

    /// Returns the first row with a nonzero entry at the given column with the guarantee that the row returned is larger than the provided `minimum_row`
    fn find_first_nonzero_row(&self, col_ix: usize, minimum_row: usize) -> Option<usize> {
        if minimum_row > self.n_rows - 1 {
            return None;
        }
        if col_ix > self.n_cols - 1 {
            return None;
        }
        let mut ret = None;
        for row_ix in minimum_row..self.n_rows {
            if self.get(row_ix, col_ix) != 0 {
                ret = Some(row_ix);
                break;
            }
        }
        ret
    }

    pub fn swap_rows(&mut self, row_1: usize, row_2: usize) {
        if row_1 == row_2 {
            return;
        }
        match self.memory_layout {
            MemoryLayout::RowMajor => self.swap_sections(row_1, row_2),
            MemoryLayout::ColMajor => {
                for (_, col) in self.ix_to_section.iter_mut() {
                    col.swap(row_1, row_2);
                }
            }
        }
    }

    pub fn swap_cols(&mut self, col_1: usize, col_2: usize) {
        match self.memory_layout {
            MemoryLayout::ColMajor => self.swap_sections(col_1, col_2),
            MemoryLayout::RowMajor => {
                for (_, row) in self.ix_to_section.iter_mut() {
                    row.swap(col_1, col_2);
                }
            }
        }
    }

    fn swap_sections(&mut self, section_1: usize, section_2: usize) {
        if section_1 == section_2 {
            return;
        }
        match self.memory_layout {
            MemoryLayout::RowMajor => {
                if section_1 >= self.n_rows || section_2 >= self.n_rows {
                    panic!("Attempting to swap rows that are outside the current matrix.")
                }
            }
            MemoryLayout::ColMajor => {
                if section_1 >= self.n_cols || section_2 >= self.n_cols {
                    panic!("Attempting to swap cols that are outside the current matrix.")
                }
            }
        }
        let tmp_1 = self.ix_to_section.remove_entry(&section_1);
        let tmp_2 = self.ix_to_section.remove_entry(&section_2);
        if let Some((_, entries_1)) = tmp_1 {
            self.ix_to_section.insert(section_2, entries_1);
        }
        if let Some((_, entries_2)) = tmp_2 {
            self.ix_to_section.insert(section_1, entries_2);
        }
    }

    /// adds `source * scalar` to `target`
    fn add_multiple_of_section_to_other(&mut self, source: usize, target: usize, scalar: FFRep) {
        if let Some(source_section) = self.ix_to_section.get(&source) {
            let mut summand = source_section.clone();
            summand.scale(scalar, self.field_mod);
            if let Some(target_section) = self.ix_to_section.get_mut(&target) {
                target_section.add_to_self(&summand, self.field_mod);
            } else {
                self.ix_to_section.insert(target, summand);
            }
        }
    }

    /// Finds the next pivot entry given the previous one. Will swap rows
    /// to prep matrix for next round.
    fn find_next_pivot_with_swap(
        &mut self,
        previous_pivot: (usize, usize),
    ) -> Option<(usize, usize)> {
        if previous_pivot.0 == self.n_rows - 1 {
            return None;
        }
        if previous_pivot.1 == self.n_cols - 1 {
            return None;
        }
        let pivot_row = previous_pivot.0 + 1;
        for col_ix in (previous_pivot.1 + 1)..self.n_cols {
            if let Some(nonzero_row) = self.find_first_nonzero_row(col_ix, pivot_row) {
                self.swap_sections(pivot_row, nonzero_row);
                if self.get(pivot_row, col_ix) != 0 {
                    return Some((pivot_row, col_ix));
                }
            }
        }
        None
    }

    fn reduce_column_from_pivot(&mut self, pivot: (usize, usize)) {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        let pivot_entry = self.get(pivot.0, pivot.1);
        if pivot_entry == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry != 1 {
            let pivot_inv = FF::from((pivot_entry, self.field_mod)).modular_inverse();
            self.scale_section(pivot.0, pivot_inv.0);
        }
        for row_ix in 0..self.n_rows {
            if row_ix == pivot.0 {
                continue;
            }
            let entry = FF::from((self.get(row_ix, pivot.1), self.field_mod));
            let scalar = -1 * entry;
            if scalar.0 != 0 {
                self.add_multiple_of_section_to_other(pivot.0, row_ix, scalar.0);
            }
        }
    }

    pub fn rref(&mut self) {
        if self.memory_layout != MemoryLayout::RowMajor {
            panic!("Trying to put sparse matrix in reduced row echelon form when matrix is in column major order.")
        }
        // Find the pivot column of the first row, then use helpers for the rest
        let mut first_col_ix = None;
        for col_ix in 0..self.n_cols {
            if let Some(first_nonzero_row) = self.find_first_nonzero_row(0, col_ix) {
                self.swap_sections(0, first_nonzero_row);
                first_col_ix = Some(col_ix);
                break;
            }
        }
        if first_col_ix.is_none() {
            println!("Could not find pivot for first row. Doing nothing");
            return;
        }
        let mut pivot = (0, first_col_ix.unwrap());
        self.reduce_column_from_pivot(pivot);
        while let Some(new_pivot) = self.find_next_pivot_with_swap(pivot) {
            self.reduce_column_from_pivot(new_pivot);
            pivot = new_pivot;
        }
    }

    /// Returns the created pivots. WARNING: Currently eliminates
    /// zero rows and columns and reindexes the matrix. This is bad
    /// but
    pub fn row_echelon_form(&mut self) -> Vec<(usize, usize)> {
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("Rank calculations only available for RowMajor matrices.")
        }
        let mut pivots = Vec::new();
        let mut current_pivot: Option<(usize, usize)> = self.ensure_pivot_with_swap(None);
        while current_pivot.is_some() {
            let pivot_row = self
                .ix_to_section
                .get(&current_pivot.unwrap().0)
                .unwrap()
                .clone();
            self.eliminate_rows_below_no_swap(current_pivot.unwrap(), &pivot_row);
            pivots.push(current_pivot.unwrap());
            current_pivot = self.ensure_pivot_with_swap(current_pivot);
        }
        pivots
    }

    pub fn is_col_zero(&self, col_ix: usize) -> bool {
        for (_, row) in self.ix_to_section.iter() {
            if row.query(&col_ix) != 0 {
                return false;
            }
        }
        true
    }

    /// Creates a new pivot, with row swapping, given the provided previous pivot. Row swapping is
    /// done to guarantee that the next pivot is one row below the previous pivot. If another
    /// pivot cannot be made then this returns `None`.
    pub fn ensure_pivot_with_swap(
        &mut self,
        previous_pivot: Option<(usize, usize)>,
    ) -> Option<(usize, usize)> {
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

        // make sure the pivot entry is 1
        let row = self
            .ix_to_section
            .get_mut(&new_pivot.0)
            .expect("Just pivotized?");
        let first_nonzero = row.first_nonzero().expect("Row should not be zero");
        let scalar = FF::new(first_nonzero.1, self.field_mod).modular_inverse();
        row.scale(scalar.0, scalar.1);
        Some(new_pivot)
    }

    /// Does not swap or create the pivot. If the provided input is None then
    /// this will find the first pivot (logic is slightly different)
    pub fn find_next_pivot(
        &self,
        previous_pivot: Option<(usize, usize)>,
    ) -> Option<(usize, usize)> {
        let mut most_left_row = None;
        let ix_range = match previous_pivot {
            Some((prev_row, _)) => prev_row + 1..,
            None => 0..,
        };
        for (row_ix, row) in self.ix_to_section.range(ix_range) {
            let first_nonzero = row.first_nonzero();
            if first_nonzero.is_none() {
                continue;
            }
            let (col_ix, _) = first_nonzero.unwrap();
            if let Some((_, prev_col)) = &previous_pivot {
                if col_ix <= *prev_col {
                    dbg!(previous_pivot);
                    dbg!(col_ix);
                    self.dense_print();
                    panic!("Found a row that has not been properly reduced. Pivot invariant has been broken.")
                } else if col_ix == prev_col + 1 {
                    // found it.
                    most_left_row = Some((*row_ix, col_ix));
                    break;
                }
            } else {
                if col_ix == 0 {
                    most_left_row = Some((*row_ix, col_ix));
                    break;
                }
            }

            if most_left_row.is_none() {
                most_left_row = Some((*row_ix, col_ix));
            } else if col_ix < most_left_row.unwrap().1 {
                most_left_row = Some((*row_ix, col_ix));
            }
        }
        if most_left_row.is_none() {
            // could not find a pivot to put in this row.
            return None;
        }
        most_left_row
    }

    /// Assumes pivot_row is already `1` at pivot col. Does not swap rows or columns.
    pub fn eliminate_rows_below_no_swap(
        &mut self,
        pivot: (usize, usize),
        pivot_row: &SparseVector,
    ) {
        self.ix_to_section.par_iter_mut().for_each(|(row_ix, row)| {
            if *row_ix <= pivot.0 {
                return;
            }
            row.add_scaled_row_to_self(
                -1 * FF::new(row.query(&pivot.1), self.field_mod),
                pivot_row,
            );
        });
    }

    /// Attempts to create a pivot at the provided row and column but only eliminates all rows
    /// below `row_ix`. This is done to reduce the growth of the Grassmanian. Will swap rows around
    /// to make sure that the pivot is created at `col_ix`.
    pub fn pivotize_all_rows_below(&mut self, row_ix: usize, col_ix: usize) -> bool {
        if row_ix >= self.n_rows {
            return false;
        }
        if let Some(nonzero_row) = self.find_first_nonzero_row(col_ix, row_ix) {
            self.swap_sections(row_ix, nonzero_row);
        } else {
            return false;
        }
        let pivotizing_row = self.ix_to_section.get_mut(&row_ix).unwrap();
        let inv = FF::new(pivotizing_row.query(&col_ix), self.field_mod).modular_inverse();
        pivotizing_row.scale(inv.0, inv.1);
        let pivot_row = pivotizing_row.clone();
        self.ix_to_section.par_iter_mut().for_each(|(ix, row)| {
            if *ix <= row_ix {
                return;
            }
            let current_entry = row.query(&col_ix);
            if current_entry == 0 {
                return;
            }
            let scalar = -1 * FF::new(current_entry, self.field_mod);
            row.add_scaled_row_to_self(scalar, &pivot_row);
        });
        true
    }

    pub fn rank(&self) -> usize {
        let mut new = self.clone();
        new.rref();
        // TODO: could probably just check the non-zero diagonal elements,
        // but then the pivot may be further off the diagonal. So you'd have to
        // keep a running column index and sweep it to the right off the diagonal and once it hits the end it
        // will just stay there. but thats kind of complicated for me rn.
        let mut num_zero_rows = 0;
        for row_ix in 0..new.n_rows {
            if new.check_row_is_zero(row_ix) {
                num_zero_rows += 1;
            }
        }
        new.n_rows - num_zero_rows
    }

    fn check_row_is_zero(&self, row_ix: usize) -> bool {
        match self.memory_layout {
            MemoryLayout::RowMajor => {
                if let Some(row) = self.ix_to_section.get(&row_ix) {
                    row.is_zero()
                } else {
                    false
                }
            }
            MemoryLayout::ColMajor => {
                let mut are_all_zero = true;
                for col_ix in 0..self.n_cols {
                    if self.get(row_ix, col_ix) != 0 {
                        are_all_zero = false;
                        break;
                    }
                }
                are_all_zero
            }
        }
    }

    fn scale_section(&mut self, section: usize, scalar: FFRep) {
        if let Some(memory_section) = self.ix_to_section.get_mut(&section) {
            memory_section.scale(scalar, self.field_mod);
        }
    }

    pub fn pivotize_with_row(&mut self, pivot: (usize, usize), pivot_row: SparseVector) {
        for (ix, row) in self.ix_to_section.iter_mut() {
            if *ix == pivot.0 {
                let entry = row.query(&pivot.1);
                if entry == 0 {
                    panic!("Matrix found row I am supposed to be pivoting on, but the pivot entry is zero.");
                }
                let scalar = FF::new(entry, self.field_mod).modular_inverse().0;
                row.scale(scalar, self.field_mod);
            }
            let current_entry = row.query(&pivot.1);
            if current_entry == 0 {
                continue;
            }
            let scalar = -1 * FF::new(current_entry, self.field_mod);
            row.add_scaled_row_to_self(scalar, &pivot_row);
        }
    }

    /// shrinks `n_rows` and `n_cols` to what is actually being used
    fn shrink_to_fit(&mut self) {
        let mut max_row_ix = usize::MIN;
        let mut max_col_ix = usize::MIN;
        for section in self.ix_to_section.iter() {
            let (row, col) = self
                .memory_layout
                .get_row_col(*section.0, section.1.max_index());
            max_row_ix = max_row_ix.max(row);
            max_col_ix = max_col_ix.max(col);
        }
        self.n_rows = max_row_ix + 1;
        self.n_cols = max_col_ix + 1;
    }

    /// Finds the first row that has a nonzero entry in the specified column
    /// within the provided row index range.
    pub fn find_nonzero_entry_among_rows(
        &self,
        col_ix: usize,
        row_ix_range: Vec<usize>,
    ) -> Option<usize> {
        if row_ix_range.len() == 0 {
            log::debug!("find_nonzero_entry_among_rows: Empty range provided.");
            return None;
        } else if row_ix_range.len() == 1 {
            if self.get(row_ix_range[0], col_ix) != 0 {
                return Some(row_ix_range[0]);
            }
        }
        if col_ix > self.n_cols - 1 {
            return None;
        }
        for row_ix in row_ix_range {
            if self.get(row_ix, col_ix) != 0 {
                return Some(row_ix);
            }
        }
        None
    }

    /// scales the given row by the multiplicative inverse of the first nonzero entry.
    /// Assumes the user knows what they are doing!
    pub fn create_pivot_in_row(&mut self, row_ix: usize) -> bool {
        if self.ix_to_section.contains_key(&row_ix) == false {
            return false;
        }
        let entry = FF::new(
            self.ix_to_section
                .get(&row_ix)
                .unwrap()
                .first_nonzero()
                .unwrap()
                .1,
            self.field_mod,
        );
        self.scale_section(row_ix, entry.modular_inverse().0);
        true
    }

    /// Returns the FFMatrix corresponding to the block given by the two
    /// corners `corner1` and `corner2`. If the corners are the same
    /// then it returns a 1x1 matrix corresponding to the entry at that position.
    /// Note corners are inclusive, so `mat.get_block((0,0), (n_rows -1, n_cols - 1))` should be equal to the original matrix.
    /// corners should be given as (row_1, col_1) and (row_2, col_2)
    /// corners are inclusive! meaning `corner1 = (row, col)` and `corner2 = (row, col)` will give
    /// a single entry
    pub fn clone_block(&self, corner1: (usize, usize), corner2: (usize, usize)) -> SparseFFMatrix {
        if corner1 == corner2 {
            let entry = self.get(corner1.0, corner1.1);
            return SparseFFMatrix::new_with_entries(
                1,
                1,
                self.field_mod,
                self.memory_layout.clone(),
                vec![(0, 0, entry)],
            );
        }
        let top_row = usize::min(corner1.0, corner2.0);
        let bot_row = usize::max(corner1.0, corner2.0);
        let left_col = usize::min(corner1.1, corner2.1);
        let right_col = usize::max(corner1.1, corner2.1);
        if bot_row > self.n_rows - 1 {
            panic!("Row indexing out of bounds for getting block of a matrix.");
        }
        if right_col > self.n_cols - 1 {
            println!(
                "Columns improperly indexed. column 1: {:}, column 2: {:}, number columns: {:}",
                left_col, right_col, self.n_cols
            );
            panic!("Column indexing out of bounds for getting block of a matrix.");
        }
        // let num_rows = bot_row - top_row + 1;
        // let num_cols = right_col - left_col + 1;
        let mut entries = Vec::new();
        for row_ix in top_row..=bot_row {
            for col_ix in left_col..right_col {
                let q = self.get(row_ix, col_ix);
                entries.push((row_ix, col_ix, q));
            }
        }
        SparseFFMatrix::new_with_entries(
            bot_row - top_row + 1,
            right_col - left_col + 1,
            self.field_mod,
            self.memory_layout.clone(),
            entries,
        )
    }

    pub fn dense_print(&self) {
        let mut s = String::new();
        for row_ix in 0..self.n_rows {
            for col_ix in 0..self.n_cols {
                let q = self.get(row_ix, col_ix);
                s.push(q.to_string().chars().nth(0).unwrap());
            }
            s.push('\n');
        }
        println!("{s}");
    }

    pub fn to_dense(mut self) -> FFMatrix {
        self.shrink_to_fit();
        let zeros = (0..self.n_rows * self.n_cols)
            .into_iter()
            .map(|_| FF::new(0, self.field_mod))
            .collect();
        let mut ret = FFMatrix::new(zeros, self.n_rows, self.n_cols);
        let memory_layout = self.memory_layout.clone();
        for (section_ix, section) in self.ix_to_section.into_iter() {
            for (position_ix, entry) in section.0.into_iter() {
                let (row_ix, col_ix) = memory_layout.get_row_col(section_ix, position_ix);
                ret.set_entry(row_ix, col_ix, FF::new(entry, self.field_mod));
            }
        }
        ret
    }
}

impl RankMatrix for SparseFFMatrix {
    fn new(field_mod: u32) -> Self {
        SparseFFMatrix::new(0, 0, field_mod, MemoryLayout::RowMajor)
    }

    fn split(&mut self, row_ixs: HashSet<usize>) -> Self {
        let mut new_ix_to_section = BTreeMap::new();
        for (ix, row) in self.ix_to_section.iter() {
            if row_ixs.contains(ix) {
                new_ix_to_section.insert(*ix, row.clone());
            }
        }
        for ix in row_ixs {
            self.ix_to_section.remove(&ix);
        }
        Self {
            ix_to_section: new_ix_to_section,
            n_rows: self.n_rows,
            n_cols: self.n_cols,
            field_mod: self.field_mod,
            memory_layout: self.memory_layout.clone(),
        }
    }

    /// interleaves the rows across separate threads to avoid "chunking"
    /// where all upper rows are done with work and the remaining threads have a bunch of work
    /// to do.
    fn split_into_parallel(
        &mut self,
        row_ixs: HashSet<usize>,
        num_threads: usize,
    ) -> ParallelFFMatrix {
        let mut row_batches: Vec<HashSet<usize>> = Vec::new();
        for _ in 0..num_threads {
            row_batches.push(HashSet::with_capacity(row_ixs.len() / num_threads));
        }
        let mut batch_ix = 0;
        for row_ix in row_ixs {
            let batch = row_batches.get_mut(batch_ix).unwrap();
            batch.insert(row_ix);
            batch_ix += 1;
            batch_ix %= row_batches.len();
        }
        let mats: Vec<SparseFFMatrix> = row_batches
            .into_iter()
            .map(|ix_batch: HashSet<usize>| self.split(ix_batch))
            .collect();
        ParallelFFMatrix::new(mats)
    }

    /// Creates a new pivot across all matrix rows at the returned column, or the row is linearly
    /// dependent on the other rows in the matrix.
    fn pivotize_row(&mut self, pivot_row_ix: usize) -> Option<usize> {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("eliminating rows with column major matrix is not supported.");
        }
        let pivot_row = self.ix_to_section.get_mut(&pivot_row_ix);
        if pivot_row.is_none() {
            return None;
        }
        let pivot_row = pivot_row.unwrap();
        let col_ix = pivot_row.first_nonzero();
        if col_ix.is_none() {
            return None;
        }
        let pivot_col_ix = col_ix.unwrap().0;
        self.eliminate_all_rows((pivot_row_ix, pivot_col_ix));
        Some(pivot_col_ix)
    }

    /// Returns the column of the newly created pivot if one is possible, returns None if the row
    /// is zero.
    fn pivotize_row_within_range(
        &mut self,
        pivot_row_ix: usize,
        range: impl AsRef<[usize]>,
    ) -> Option<usize> {
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("Eliminating rows with column major matrix not supported.")
        }
        let pivot_row = self.ix_to_section.get_mut(&pivot_row_ix);
        if pivot_row.is_none() {
            return None;
        }
        let pivot_row = pivot_row.unwrap();
        let nonzero = pivot_row.first_nonzero();
        if nonzero.is_none() {
            return None;
        }
        let (pivot_col_ix, entry) = nonzero.unwrap();
        let mod_inv = FF::new(entry, self.field_mod).modular_inverse();
        pivot_row.scale(mod_inv.0, mod_inv.1);
        let pivot_row = pivot_row.clone();
        let mut zero_rows = Vec::new();
        for ix in range.as_ref().iter() {
            let row = self.ix_to_section.get_mut(&ix);
            if row.is_none() {
                continue;
            }
            let row = row.unwrap();
            if *ix == pivot_row_ix {
                continue;
            }
            let current_entry = row.query(&pivot_col_ix);
            if current_entry == 0 {
                continue;
            }
            let scalar = -1 * FF::new(current_entry, self.field_mod);
            row.add_scaled_row_to_self(scalar, &pivot_row);
            if row.is_zero() {
                zero_rows.push(ix);
            }
        }
        for zero_row_ix in zero_rows {
            self.ix_to_section.remove(&zero_row_ix);
        }
        Some(pivot_col_ix)
    }

    fn eliminate_rows(&mut self, pivot: (usize, usize), rows_to_eliminate: impl AsRef<[usize]>) {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        let pivot_entry = self.get(pivot.0, pivot.1);
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("eliminating rows with column major matrix is not supported.");
        }
        if pivot_entry == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry != 1 {
            let pivot_inv = FF::from((pivot_entry, self.field_mod)).modular_inverse();
            self.scale_section(pivot.0, pivot_inv.0);
        }
        let pivot_row = self.ix_to_section.get(&pivot.0).unwrap().clone();
        for row_ix in rows_to_eliminate.as_ref().iter() {
            if *row_ix == pivot.0 {
                continue;
            }
            let row = self.ix_to_section.get_mut(&row_ix);
            if row.is_none() {
                continue;
            }
            let row = row.unwrap();
            let entry = row.query(&pivot.1);
            let scalar = -1 * FF::new(entry, self.field_mod);
            row.add_scaled_row_to_self(scalar, &pivot_row);
        }
    }

    fn eliminate_all_rows(&mut self, pivot: (usize, usize)) {
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("Eliminating all rows with a column-major matrix is not implemented.")
        }
        let pivot_entry = FF::from((self.get(pivot.0, pivot.1), self.field_mod));
        if pivot_entry.0 == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry.0 != 1 {
            let pivot_inv = pivot_entry.modular_inverse();
            self.scale_section(pivot.0, pivot_inv.0);
        }
        let pivot_row = self.ix_to_section.get(&pivot.0).unwrap().clone();
        for (ix, row) in self.ix_to_section.iter_mut() {
            if *ix == pivot.0 {
                continue;
            }
            let current_entry = row.query(&pivot.1);
            if current_entry == 0 {
                continue;
            }
            let scalar = -1 * FF::new(current_entry, self.field_mod);
            row.add_scaled_row_to_self(scalar, &pivot_row);
        }
    }

    fn insert(&mut self, row_ix: usize, col_ix: usize, entry: FFRep) {
        let (section, position) = self.memory_layout.get_section_position(row_ix, col_ix);
        if let Some(memory_section) = self.ix_to_section.get_mut(&section) {
            memory_section.insert(position, entry % self.field_mod);
        } else {
            let new_section =
                SparseVector::new_with_entries(vec![(position, entry % self.field_mod)]);
            self.ix_to_section.insert(section, new_section);
        }
        self.n_rows = self.n_rows.max(row_ix + 1);
        self.n_cols = self.n_cols.max(col_ix + 1);
    }

    fn get(&self, row_ix: usize, col_ix: usize) -> FFRep {
        let (section, position) = self.memory_layout.get_section_position(row_ix, col_ix);
        if let Some(memory_section) = self.ix_to_section.get(&section) {
            memory_section.query(&position)
        } else {
            0
        }
    }

    fn to_disk(self, filename: impl AsRef<Path>) -> Result<(), ()> {
        SparseFFMatrix::to_disk(&self, filename.as_ref());
        Ok(())
    }

    fn from_disk(path: impl AsRef<Path>) -> Option<std::rc::Rc<Self>> {
        Some(Rc::new(Self::from_disk(path.as_ref())))
    }

    fn n_rows(&self) -> usize {
        self.n_rows
    }

    fn n_cols(&self) -> usize {
        self.n_cols
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use rand::thread_rng;

    use crate::matrices::{mat_trait::RankMatrix, sparse_ffmatrix::ParallelFFMatrix};

    use super::SparseFFMatrix;

    #[test]
    fn test_rref() {
        let p = 7;
        let mut mat = crate::matrices::sparse_ffmatrix::SparseFFMatrix::new(
            3,
            3,
            p,
            super::MemoryLayout::RowMajor,
        );
        let entries = vec![
            (0, 0, 1),
            (0, 1, 2),
            (0, 2, 3),
            (1, 0, 4),
            (1, 1, 5),
            (1, 2, 6),
            (2, 0, 7),
            (2, 1, 8),
            (2, 2, 9),
        ];
        mat.insert_entries(entries);
        mat.rref();
        let mut old_stuff = crate::matrices::ffmatrix::FFMatrix::new(
            (1..=9)
                .into_iter()
                .map(|x| crate::math::finite_field::FiniteField::new(x, p))
                .collect(),
            3,
            3,
        );
        old_stuff.rref();
        let dense = mat.to_dense();
        println!("new: {:}", dense);
        println!("old: {:}", old_stuff);
    }

    #[test]
    fn various_pivoting() {
        let p = 7;
        let mut mat = crate::matrices::sparse_ffmatrix::SparseFFMatrix::new(
            3,
            3,
            p,
            super::MemoryLayout::RowMajor,
        );
        let entries = vec![
            (0, 0, 1),
            (0, 1, 2),
            (0, 2, 3),
            (1, 0, 4),
            (1, 1, 5),
            (1, 2, 6),
            (2, 0, 7),
            (2, 1, 8),
            (2, 2, 9),
        ];
        mat.insert_entries(entries);
        println!("pre: {:}", mat.clone().to_dense());
        let col = mat.pivotize_row_within_range(0, vec![0, 1]);
        println!("mat: {:}", mat.to_dense());
        dbg!(col);
    }

    #[test]
    fn upper_triangular() {
        let mut mat = SparseFFMatrix::new_random(40, 40, 7, 0.05);
        // mat.insert_entries(vec![(0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)]);
        mat.dense_print();
        let mut parallel_mat = mat.clone().split_into_parallel(
            mat.ix_to_section.keys().cloned().collect(),
            std::thread::available_parallelism().unwrap().into(),
        );
        let parallel_pivots = parallel_mat.row_echelon_form();
        let rayon_pivots = mat.row_echelon_form();
        assert_eq!(parallel_pivots.len(), rayon_pivots.len());
        for ix in 0..parallel_pivots.len() {
            println!("{:?} == {:?}", parallel_pivots[ix], rayon_pivots[ix]);
        }
        parallel_mat.quit().dense_print();
        mat.dense_print();
    }

    #[test]
    fn disk_loading_and_retrieving() {
        let mut entries: Vec<(usize, usize, u32)> = Vec::new();
        for ix in 0..10000 {
            entries.push((ix % 100, ix, (ix % 11) as u32))
        }
        let mat = crate::matrices::sparse_ffmatrix::SparseFFMatrix::new_with_entries(
            100,
            10000,
            11,
            super::MemoryLayout::RowMajor,
            entries,
        );
        let filename = PathBuf::from("/Users/matt/repos/qec/tmp/big_mat.txt");
        mat.clone().to_disk(&filename).expect("ERROR!");
        let loaded = SparseFFMatrix::from_disk(&filename);
        assert_eq!(loaded, mat);
    }

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
}
