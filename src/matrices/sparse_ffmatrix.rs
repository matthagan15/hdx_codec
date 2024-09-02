use core::panic;
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{Read, Write},
    path::Path,
};

use serde::{Deserialize, Serialize};

use super::ffmatrix::FFMatrix;
use crate::math::finite_field::{FFRep, FiniteField as FF};

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
/// Meant to be used as a row or column of a sparse matrix. Stores the index and
/// value of the nonzero entries in the vector. The dimension of the vector
/// is not stored, although that probably wouldn't be too hard to add.
pub struct SparseVector(Vec<(usize, FFRep)>);
impl SparseVector {
    pub fn new_empty() -> Self {
        SparseVector(Vec::new())
    }

    pub fn first_nonzero(&self) -> Option<(usize, FFRep)> {
        self.0.first().cloned()
    }

    pub fn to_vec(self) -> Vec<(usize, FFRep)> {
        self.0
    }

    pub fn new_with_entries(entries: Vec<(usize, FFRep)>) -> Self {
        let mut new_entries = entries;
        new_entries.sort_by(|a, b| a.0.cmp(&b.0));
        SparseVector(new_entries)
    }

    pub fn is_zero(&self) -> bool {
        self.0.len() == 0
    }

    pub fn max_index(&self) -> usize {
        self.0.iter().fold(0, |acc, x| acc.max(x.0))
    }

    pub fn clone_and_scale(&self, scalar: FF) -> Self {
        let v: Vec<(usize, FFRep)> = self
            .0
            .iter()
            .map(|(ix, entry)| (*ix, (&scalar * *entry).0))
            .collect();
        SparseVector(v)
    }

    /// Inserts the given `entry` at position `ix`, overwriting existing
    /// entries if they exist.
    pub fn insert(&mut self, ix: usize, entry: FFRep) {
        let insertion_point = self.0.partition_point(|&x| x.0 < ix);
        if insertion_point >= self.0.len() {
            self.0.insert(insertion_point, (ix, entry));
        } else if self.0[insertion_point].0 == ix {
            self.0[insertion_point].1 = entry;
        } else {
            self.0.insert(insertion_point, (ix, entry));
        }
    }
    pub fn query(&self, ix: &usize) -> FFRep {
        if let Ok(position) = self.0.binary_search_by_key(ix, |&(ix, _)| ix) {
            self.0[position].1
        } else {
            0
        }
    }

    pub fn sparsity(&self) -> usize {
        self.0.len()
    }

    pub fn scale(&mut self, scalar: FFRep, field_mod: FFRep) {
        for (_, entry) in self.0.iter_mut() {
            *entry *= scalar;
            *entry %= field_mod;
        }
    }

    pub fn add_to_self(&mut self, rhs: &SparseVector, field_mod: FFRep) {
        let mut hm: HashMap<usize, u32> = self.0.clone().into_iter().collect();
        for (ix, entry) in rhs.0.iter() {
            let e = hm.entry(*ix).or_default();
            *e += *entry;
            *e %= field_mod;
        }
        let mut v: Vec<(usize, FFRep)> =
            hm.into_iter().filter(|(_ix, entry)| *entry != 0).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        self.0 = v;
    }

    pub fn add_scaled_row_to_self(&mut self, scalar: FF, rhs: &SparseVector) {
        let mut self_ix = 0;
        let mut rhs_ix = 0;
        let mut new_vec = Vec::new();
        loop {
            if self_ix == self.0.len() && rhs_ix < rhs.0.len() {
                let (rhs_pos, rhs_entry) = rhs.0[rhs_ix];
                new_vec.push((rhs_pos, (&scalar * rhs_entry).0));
                rhs_ix += 1;
            } else if self_ix < self.0.len() && rhs_ix == rhs.0.len() {
                new_vec.push(self.0[self_ix]);
                self_ix += 1;
            } else if self_ix == self.0.len() && rhs_ix == rhs.0.len() {
                break;
            } else {
                let (self_pos, self_entry) = self.0[self_ix];
                let (rhs_pos, rhs_entry) = rhs.0[rhs_ix];
                if self_pos == rhs_pos {
                    let new_entry = &(&scalar * rhs_entry) + self_entry;
                    if new_entry.0 != 0 {
                        new_vec.push((self_pos, new_entry.0));
                    }
                    self_ix += 1;
                    rhs_ix += 1;
                } else if self_pos < rhs_pos {
                    new_vec.push((self_pos, self_entry));
                    self_ix += 1;
                } else {
                    new_vec.push((rhs_pos, (&scalar * rhs_entry).0));
                    rhs_ix += 1;
                }
            }
        }
        self.0 = new_vec;
    }
}

/// Used to determine if a memory layout for a matrix (posibly extend to tensor?)
/// is RowMajor or ColMajor, used to determine if the sparse sections
/// retrieved are rows or columns.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MemoryLayout {
    RowMajor,
    ColMajor,
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

    /// Takes the given rows in `row_ixs` from `self` and returns them in a new matrix.
    /// Row indices are kept the same across both.
    pub fn split(&mut self, row_ixs: impl AsRef<[usize]>) -> Self {
        let ret_ix_to_section = BTreeMap::new();
        self.ix_to_section
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
                        section.insert(row_ix, entry.0);
                    }
                    MemoryLayout::ColMajor => {
                        let section = ix_to_sections
                            .entry(row_ix)
                            .or_insert(SparseVector::new_empty());
                        section.insert(col_ix, entry.0);
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

    fn get_section_position(&self, row_ix: usize, col_ix: usize) -> (usize, usize) {
        match self.memory_layout {
            MemoryLayout::RowMajor => (row_ix, col_ix),
            MemoryLayout::ColMajor => (col_ix, row_ix),
        }
    }

    /// Gets the row and column indices associated with the given `section` and `position`
    pub fn get_row_col(&self, section: usize, position: usize) -> (usize, usize) {
        match self.memory_layout {
            MemoryLayout::RowMajor => (section, position),
            MemoryLayout::ColMajor => (position, section),
        }
    }

    pub fn insert(&mut self, row_ix: usize, col_ix: usize, entry: FFRep) {
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

    /// Adds a collection of (`row_ix`, `col_ix`, `entry`) to the matrix.
    pub fn insert_entries(&mut self, entries: Vec<(usize, usize, FFRep)>) {
        for (row_ix, col_ix, entry) in entries.into_iter() {
            self.insert(row_ix, col_ix, entry);
        }
    }

    pub fn get(&self, row_ix: usize, col_ix: usize) -> FF {
        let (section, position) = self.memory_layout.get_section_position(row_ix, col_ix);
        if let Some(memory_section) = self.ix_to_section.get(&section) {
            (memory_section.query(&position), self.field_mod).into()
        } else {
            (0, self.field_mod).into()
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
            if self.get(row_ix, col_ix).0 != 0 {
                ret = Some(row_ix);
                break;
            }
        }
        ret
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
        match (tmp_1.is_some(), tmp_2.is_some()) {
            (true, true) => {
                self.ix_to_section.insert(section_2, tmp_1.unwrap().1);
                self.ix_to_section.insert(section_1, tmp_2.unwrap().1);
            }
            (true, false) => {
                self.ix_to_section.insert(section_2, tmp_1.unwrap().1);
            }
            (false, true) => {
                self.ix_to_section.insert(section_1, tmp_2.unwrap().1);
            }
            (false, false) => {}
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
    fn find_next_pivot(&mut self, previous_pivot: (usize, usize)) -> Option<(usize, usize)> {
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
                if self.get(pivot_row, col_ix).0 != 0 {
                    return Some((pivot_row, col_ix));
                }
            }
        }
        None
    }

    fn reduce_column_from_pivot(&mut self, pivot: (usize, usize)) {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        let pivot_entry = self.get(pivot.0, pivot.1);
        if pivot_entry.0 == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry.0 != 1 {
            let pivot_inv = pivot_entry.modular_inverse();
            self.scale_section(pivot.0, pivot_inv.0);
        }
        for row_ix in 0..self.n_rows {
            if row_ix == pivot.0 {
                continue;
            }
            let entry = self.get(row_ix, pivot.1);
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
        while let Some(new_pivot) = self.find_next_pivot(pivot) {
            self.reduce_column_from_pivot(new_pivot);
            pivot = new_pivot;
        }
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
                    if self.get(row_ix, col_ix).0 != 0 {
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
            if self.get(row_ix_range[0], col_ix).0 != 0 {
                return Some(row_ix_range[0]);
            }
        }
        if col_ix > self.n_cols - 1 {
            return None;
        }
        for row_ix in row_ix_range {
            if self.get(row_ix, col_ix).0 != 0 {
                return Some(row_ix);
            }
        }
        None
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
                vec![(0, 0, entry.0)],
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
                entries.push((row_ix, col_ix, q.0));
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
                s.push(q.0.to_string().chars().nth(0).unwrap());
            }
            s.push('\n');
        }
        println!("{s}");
    }

    /// Returns the column of the newly created pivot if one is possible, returns None if the row
    /// is zero.
    pub fn pivotize_row_within_range(
        &mut self,
        pivot_row_ix: usize,
        row_ix_range: impl AsRef<[usize]>,
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
        for ix in row_ix_range.as_ref().iter() {
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

    pub fn eliminate_rows(&mut self, pivot: (usize, usize), row_ix_range: Vec<usize>) {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        let pivot_entry = self.get(pivot.0, pivot.1);
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("eliminating rows with column major matrix is not supported.");
        }
        if pivot_entry.0 == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry.0 != 1 {
            let pivot_inv = pivot_entry.modular_inverse();
            self.scale_section(pivot.0, pivot_inv.0);
        }
        let pivot_row = self.ix_to_section.get(&pivot.0).unwrap().clone();
        for row_ix in row_ix_range {
            if row_ix == pivot.0 {
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

    /// Creates a new pivot across all matrix rows at the returned column, or the row is linearly
    /// dependent on the other rows in the matrix.
    pub fn pivotize_row(&mut self, pivot_row_ix: usize) -> Option<usize> {
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

    /// *panics* if the provided entry is zero.
    pub fn eliminate_all_rows(&mut self, pivot: (usize, usize)) {
        if self.memory_layout == MemoryLayout::ColMajor {
            panic!("Eliminating all rows with a column-major matrix is not implemented.")
        }
        let pivot_entry = self.get(pivot.0, pivot.1);
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
            let mut adder_row = pivot_row.clone();
            adder_row.scale(current_entry, self.field_mod);
            row.add_to_self(&adder_row, self.field_mod);
        }
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

mod tests {
    use std::path::{Path, PathBuf};

    use super::SparseFFMatrix;

    #[test]
    fn test_sparse_section() {
        let v = vec![(0, 10), (200, 4), (10, 30), (1000, 2)];
        let mut ss = crate::matrices::sparse_ffmatrix::SparseVector::new_with_entries(v);
        assert_eq!(ss.query(&7), 0);
        assert_eq!(ss.query(&1000), 2);
        ss.insert(4, 17);
        assert_eq!(ss.query(&4), 17);
    }

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
    fn sparse_vector_additions() {
        let p = 7_u32;
        let mut v = Vec::new();
        for ix in 0..10 {
            v.push((ix as usize, 1));
        }
        let mut u = Vec::new();
        for ix in 10..20 {
            u.push((ix as usize, 1));
        }

        let adder = crate::matrices::sparse_ffmatrix::SparseVector::new_with_entries(vec![
            (3, 4),
            (4, 6),
            (9, 3),
        ]);
        let mut sparse_v = crate::matrices::sparse_ffmatrix::SparseVector::new_with_entries(v);
        let mut sparse_u = crate::matrices::sparse_ffmatrix::SparseVector::new_with_entries(u);
        sparse_v.add_scaled_row_to_self(crate::math::finite_field::FiniteField::new(1, p), &adder);
        sparse_u
            .add_scaled_row_to_self(crate::math::finite_field::FiniteField::new(1, p), &sparse_v);
        dbg!(sparse_u);
        dbg!(sparse_v);
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
        mat.to_disk(&filename);
        let loaded = SparseFFMatrix::from_disk(&filename);
        assert_eq!(loaded, mat);
    }
}
