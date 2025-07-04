use std::{
    fmt::{Display, Write},
    ops::{Index, Mul},
};

use serde::{ser::SerializeStruct, Deserialize, Serialize};

use crate::math::finite_field::{FFRep, FiniteField as FF};

/// Constructs the Vandermonde matrix V such that V_{i,j} = elements[j]^i.
pub fn vandermonde(elements: &Vec<FFRep>, n_rows: usize, field_mod: FFRep) -> FFMatrix {
    let mut entries = Vec::new();
    for k in 0..n_rows {
        let mut tmp = elements
            .clone()
            .into_iter()
            .map(|x| x.pow(k as u32) % field_mod)
            .collect();
        entries.append(&mut tmp);
    }
    FFMatrix::new(entries, n_rows, elements.len(), field_mod)
}

// TODO: Why the hell did this use the entire FiniteField struct instead of just storing the
// value????????????
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FFMatrix {
    pub entries: Vec<FFRep>,
    pub n_rows: usize,
    pub n_cols: usize,
    pub field_mod: u32,
}

impl Serialize for FFMatrix {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut s = String::new();
        if self.entries.len() == 0 {
            s.push_str("[]");
        } else {
            s.push_str("[");
            for node in self.entries.iter() {
                s.push_str(&format!("{:?},", node));
            }
            if s.ends_with(',') {
                s.remove(s.len() - 1);
            }
            s.push_str("]");
        }

        let mut ser = serializer.serialize_struct("Mat", 4)?;
        ser.serialize_field("entries", &s[..])?;
        ser.serialize_field("n_rows", &self.n_rows)?;
        ser.serialize_field("n_cols", &self.n_cols)?;
        ser.serialize_field("field_mod", &self.field_mod)?;
        ser.end()
    }
}

impl<'de> Deserialize<'de> for FFMatrix {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let mut data = <serde_json::Value>::deserialize(deserializer)?;
        let mut entries_string = String::from(data["entries"].take().as_str().unwrap());
        let n_rows = data["n_rows"].take().as_u64().unwrap() as usize;
        let n_cols = data["n_cols"].take().as_u64().unwrap() as usize;
        let field_mod = data["field_mod"].take().as_u64().unwrap() as u32;

        if entries_string.starts_with("[") {
            entries_string.remove(0);
        }
        if entries_string.ends_with("]") {
            entries_string.remove(entries_string.len() - 1);
        }
        if entries_string.contains(",") {
            let v: Vec<FFRep> = entries_string
                .split(',')
                .filter_map(|x| -> Option<FFRep> { x.parse().ok() })
                .collect();
            Ok(FFMatrix::new(v, n_rows, n_cols, field_mod))
        } else {
            if let Ok(n) = entries_string.parse::<u32>() {
                Ok(FFMatrix::new(vec![n], n_rows, n_cols, field_mod))
            } else {
                if entries_string.len() == 0 {
                    Ok(FFMatrix::new(vec![], n_rows, n_cols, field_mod))
                } else {
                    println!("vec_data: {:?}", entries_string);
                    panic!("Could not parse single input.");
                }
            }
        }
    }
}

impl FFMatrix {
    /// Creates a new matrix given entries in row major order. For example, the
    /// 2x2 matrix
    /// A = [[a, b], [c, d]], where the first row is a, b and the second c, d,
    ///  would be created with an entries vec `entries = vec![a, b, c, d];`
    pub fn new(entries: Vec<FFRep>, n_rows: usize, n_cols: usize, field_mod: u32) -> Self {
        if entries.len() == 0 {
            panic!("Why did you try to make an empty matrix?")
        }
        if entries.len() != n_rows * n_cols {
            panic!("Error: tried to create matrices with improper number of rows and columns.")
        }
        Self {
            entries,
            n_rows,
            n_cols,
            field_mod,
        }
    }

    pub fn new_from_entries(entries: Vec<(usize, usize, FFRep)>, field_mod: FFRep) -> Self {
        let mut n_rows = 0;
        let mut n_cols = 0;
        for (row_ix, col_ix, _entry) in entries.iter() {
            n_rows = n_rows.max(row_ix + 1);
            n_cols = n_cols.max(col_ix + 1);
        }
        let mut ret = FFMatrix::zero(n_rows, n_cols, field_mod);
        for (row_ix, col_ix, entry) in entries {
            let ix = ret.convert_indices(row_ix, col_ix);
            let new_entry = ret.entries[ix] + entry;
            ret.entries[ix] = new_entry % field_mod;
        }
        ret
    }

    pub fn is_square(&self) -> bool {
        self.n_rows == self.n_cols
    }

    pub fn rref(&mut self) {
        // Find the pivot column of the first row, then use helpers for the rest
        let mut first_col_ix = None;
        for col_ix in 0..self.n_cols {
            if let Some(first_nonzero_row) = self.find_first_nonzero_row(0, col_ix) {
                self.swap_rows(0, first_nonzero_row);
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

    pub fn row_echelon_form(&mut self) {
        // Find the pivot column of the first row, then use helpers for the rest
        let mut first_col_ix = None;
        for col_ix in 0..self.n_cols {
            if let Some(first_nonzero_row) = self.find_first_nonzero_row(0, col_ix) {
                self.swap_rows(0, first_nonzero_row);
                first_col_ix = Some(col_ix);
                break;
            }
        }
        if first_col_ix.is_none() {
            println!("Could not find pivot for first column. Doing nothing");
            return;
        }
        let mut pivot = (0, first_col_ix.unwrap());
        self.reduce_column_below_pivot(pivot);
        while let Some(new_pivot) = self.find_next_pivot(pivot) {
            self.reduce_column_below_pivot(new_pivot);
            pivot = new_pivot;
        }
    }

    /// Computes the matrix that acts on the row space of self to
    /// put self in rref
    pub fn rref_left_inverse(&self) -> FFMatrix {
        let mut self_with_identity_entries = Vec::new();
        let id = FFMatrix::id(self.n_rows, self.field_mod);
        for row_ix in 0..self.n_rows {
            self_with_identity_entries.append(&mut self.get_row(row_ix));
            self_with_identity_entries.append(&mut id.get_row(row_ix));
        }
        let mut self_with_identity = FFMatrix::new(
            self_with_identity_entries,
            self.n_rows,
            self.n_cols + id.n_cols,
            self.field_mod,
        );
        self_with_identity.rref();
        self_with_identity.clone_block(
            (0, self.n_cols),
            (self.n_rows - 1, self_with_identity.n_cols - 1),
        )
    }

    /// WARNING: Does no checks for bounds or for entry field_mod !
    pub fn set_entry(&mut self, row_ix: usize, col_ix: usize, entry: FF) {
        let ix = self.convert_indices(row_ix, col_ix);
        self.entries[ix] = entry.0;
    }

    pub fn rank(&self) -> usize {
        let mut new = self.clone();
        new.rref();
        // TODO: could probably just check the non-zero diagonal elements,
        // but then the pivot may be further off the diagonal. So you'd have to
        // keep a running column index and sweep it to the right off the diagonal and once it hits the end it
        // will just stay there. but thats kind of complicated for me rn.
        new.n_rows - new.num_zero_rows()
    }

    pub fn num_zero_rows(&self) -> usize {
        let mut num_zero_rows = 0;
        for row_ix in 0..self.n_rows {
            if self.check_row_is_zero(row_ix) {
                num_zero_rows += 1;
            }
        }
        num_zero_rows as usize
    }

    fn check_row_is_zero(&self, row_ix: usize) -> bool {
        let ix = self.convert_indices(row_ix, 0);
        let mut are_all_zero = true;
        for col_ix in 0..self.n_cols {
            if self.entries[ix + col_ix] != 0 {
                are_all_zero = false;
                break;
            }
        }
        are_all_zero
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
                self.swap_rows(pivot_row, nonzero_row);
                if self[[pivot_row, col_ix]] != 0 {
                    return Some((pivot_row, col_ix));
                }
            }
        }
        None
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
            if self[[row_ix, col_ix]] != 0 {
                ret = Some(row_ix);
                break;
            }
        }
        ret
    }

    fn reduce_column_from_pivot(&mut self, pivot: (usize, usize)) {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        let pivot_entry = self.entries[self.convert_indices(pivot.0, pivot.1)];
        if pivot_entry == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry != 1 {
            let pivot_inv = FF::new(pivot_entry, self.field_mod).modular_inverse();
            self.scale_row(pivot.0, pivot_inv);
        }
        for row_ix in 0..self.n_rows {
            if row_ix == pivot.0 {
                continue;
            }
            let entry = self[[row_ix, pivot.1]];
            let scalar = self.field_mod - entry;
            self.add_multiple_of_row_to_other(pivot.0, row_ix, scalar);
        }
    }

    fn reduce_column_below_pivot(&mut self, pivot: (usize, usize)) {
        // First check that the pivot is 1, if not rescale. Panic if 0.
        let pivot_entry = self.entries[self.convert_indices(pivot.0, pivot.1)];
        if pivot_entry == 0 {
            panic!("Cannot reduce a column on a non-zero row.")
        } else if pivot_entry != 1 {
            let pivot_inv = FF::new(pivot_entry, self.field_mod).modular_inverse();
            self.scale_row(pivot.0, pivot_inv);
        }
        for row_ix in pivot.0 + 1..self.n_rows {
            if row_ix == pivot.0 {
                continue;
            }
            let entry = self[[row_ix, pivot.1]];
            let scalar = self.field_mod - entry;
            self.add_multiple_of_row_to_other(pivot.0, row_ix, scalar);
        }
    }

    fn scale_row(&mut self, row_ix: usize, scalar: FF) {
        let start_ix = self.convert_indices(row_ix, 0);
        for col_ix in 0..self.n_cols {
            self.entries[start_ix + col_ix] =
                (self.entries[start_ix + col_ix] * scalar.0) % self.field_mod;
        }
    }

    pub fn scale(&mut self, scalar: FF) {
        for entry in self.entries.iter_mut() {
            *entry = (*entry * scalar.0) % self.field_mod;
        }
    }

    fn add_multiple_of_row_to_other(
        &mut self,
        source_row: usize,
        target_row: usize,
        scalar: FFRep,
    ) {
        let source_start_ix = self.convert_indices(source_row, 0);
        let target_start_ix = self.convert_indices(target_row, 0);
        for col_ix in 0..self.n_cols {
            let new_entry = self.entries[target_start_ix + col_ix]
                + self.entries[source_start_ix + col_ix] * scalar;
            self.entries[target_start_ix + col_ix] = new_entry % self.field_mod;
        }
    }

    fn swap_rows(&mut self, row_1: usize, row_2: usize) {
        if row_1 == row_2 {
            return;
        }
        let r1 = row_1 % self.n_rows;
        let r2 = row_2 % self.n_rows;
        let r1_start_ix = self.convert_indices(r1, 0);
        let r2_start_ix = self.convert_indices(r2, 0);
        for k in 0..self.n_cols {
            let tmp = self.entries[r1_start_ix + k].clone();
            self.entries[r1_start_ix + k] = self.entries[r2_start_ix + k].clone();
            self.entries[r2_start_ix + k] = tmp;
        }
    }

    pub fn transpose(&mut self) {
        // WARNING: In place transpose is difficult so don't try and do it.
        // TODO: Could possibly still store matrices in row-major order,
        // just add a flag if the matrix is transposed or not and just change the indexing but that seems like a bad idea.
        let mut new_entries = Vec::new();
        for col_ix in 0..self.n_cols {
            for row_ix in 0..self.n_rows {
                let ix = self.convert_indices(row_ix, col_ix);
                let x = self.entries[ix];
                new_entries.push(x);
            }
        }
        self.entries = new_entries;
        let n = self.n_rows;
        let m = self.n_cols;
        self.n_rows = m;
        self.n_cols = n;
    }

    pub fn zero(n_rows: usize, n_cols: usize, field_mod: u32) -> Self {
        let mut entries = Vec::with_capacity(n_rows * n_cols);
        for _ in 0..(n_rows * n_cols) {
            entries.push(0);
        }
        Self {
            entries,
            n_rows,
            n_cols,
            field_mod,
        }
    }

    pub fn id(dim: usize, field_mod: u32) -> FFMatrix {
        let mut entries = Vec::with_capacity(dim * dim);
        for row_ix in 0..dim {
            for col_ix in 0..dim {
                entries.push(if row_ix == col_ix { 1 } else { 0 });
            }
        }
        FFMatrix {
            entries,
            n_cols: dim,
            n_rows: dim,
            field_mod,
        }
    }

    /// ix is the row index (starts at 0) and jx is the col index (also
    /// starts at 0)
    pub fn convert_indices(&self, ix: usize, jx: usize) -> usize {
        ((ix % self.n_rows) * self.n_cols) + (jx % self.n_cols)
    }
    /// Clones the specified row of the matrix as `FiniteFieldPolynomials`, does not include information about the quotient polynomial.
    pub fn get_row(&self, row_ix: usize) -> Vec<FFRep> {
        let row_ix_start = self.convert_indices(row_ix, 0);
        self.entries[row_ix_start..(row_ix_start + self.n_cols)]
            .iter()
            .cloned()
            .collect()
    }

    /// Clones the specified column of the matrix as `FiniteFieldPolynomials`, does not include information about the quotient polynomial.
    pub fn get_col(&self, col_ix: usize) -> Vec<FFRep> {
        let mut ret = Vec::with_capacity(self.n_rows);
        for ix in 0..self.n_rows {
            ret.push(self.entries[self.convert_indices(ix, col_ix)].clone());
        }
        ret
    }

    /// Returns the FFMatrix corresponding to the block given by the two
    /// corners `corner1` and `corner2`. If the corners are the same
    /// then it returns a 1x1 matrix corresponding to the entry at that position.
    /// Note corners are inclusive, so `mat.get_block((0,0), (n_rows -1, n_cols - 1))` should be equal to the original matrix.
    /// corners should be given as (row_1, col_1) and (row_2, col_2)
    pub fn clone_block(&self, corner1: (usize, usize), corner2: (usize, usize)) -> FFMatrix {
        if corner1 == corner2 {
            return FFMatrix {
                entries: vec![self[[corner1.0, corner1.1]]],
                n_rows: 1,
                n_cols: 1,
                field_mod: self.field_mod,
            };
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
        let num_rows = bot_row - top_row + 1;
        let num_cols = right_col - left_col + 1;
        let mut entries = Vec::new();
        for row_ix in top_row..=bot_row {
            let start_ix = self.convert_indices(row_ix, left_col);
            let mut row_slice = Vec::new();
            self.entries[start_ix..start_ix + num_cols].clone_into(&mut row_slice);
            entries.append(&mut row_slice);
        }
        FFMatrix {
            entries,
            n_rows: num_rows,
            n_cols: num_cols,
            field_mod: self.field_mod,
        }
    }
}

impl PartialOrd for FFMatrix {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.n_cols != other.n_cols {
            return None;
        }
        if self.n_rows != other.n_rows {
            return None;
        }

        for ix in 0..self.entries.len() {
            if self.entries[ix] == other.entries[ix] {
                continue;
            } else {
                return self.entries[ix].partial_cmp(&other.entries[ix]);
            }
        }
        Some(std::cmp::Ordering::Equal)
    }
}

impl Ord for FFMatrix {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other)
            .expect("Cannot order non-equal shape matrices")
    }
}

impl Mul<&Vec<FF>> for &FFMatrix {
    type Output = Vec<FF>;

    fn mul(self, rhs: &Vec<FF>) -> Self::Output {
        if rhs.len() != self.n_cols {
            println!("rhs = {:}", rhs.len());
            println!("self = {:}, {:}", self.n_rows, self.n_cols);
            panic!("Mismatched shapes in matrix-vector multiplication.")
        }
        let mut ret = Vec::new();
        let p = self.field_mod;

        for ix in 0..self.n_rows {
            let mut tmp = FF::new(0, p);
            for jx in 0..rhs.len() {
                tmp += &rhs[jx] * self.entries[self.convert_indices(ix, jx)]
            }
            ret.push(tmp);
        }
        ret
    }
}

impl Mul for &FFMatrix {
    type Output = FFMatrix;

    fn mul(self, rhs: &FFMatrix) -> Self::Output {
        if self.n_cols != rhs.n_rows {
            panic!("Tried to multiply incompatible matrices.")
        }
        let mut entries = Vec::with_capacity(self.entries.len());
        let field_mod = self.field_mod;
        for ix in 0..self.n_rows {
            for jx in 0..rhs.n_cols {
                let mut entry = 0;
                for kx in 0..self.n_cols {
                    entry += self[[ix, kx]] * rhs[[kx, jx]];
                }
                entries.push(entry % self.field_mod);
            }
        }
        FFMatrix {
            entries,
            n_rows: self.n_rows,
            n_cols: rhs.n_cols,
            field_mod,
        }
    }
}

impl Index<[usize; 2]> for FFMatrix {
    type Output = FFRep;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        self.entries
            .get(self.convert_indices(index[0], index[1]))
            .expect("Matrix indexing out of bounds.")
    }
}

impl Display for FFMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut max_string_len = 0;
        let strings: Vec<String> = self
            .entries
            .iter()
            .map(|q| {
                let s = q.to_string();
                max_string_len = max_string_len.max(s.len());
                s
            })
            .collect();
        let mut rows = Vec::new();
        let mut col_counter = 0;
        let mut row = String::from("| ");
        let mut row_len = 0;
        for mut string in strings {
            let diff_len = max_string_len - string.len();
            string.push_str(&" ".repeat(diff_len));
            col_counter += 1;
            col_counter %= self.n_cols;
            if col_counter == 0 {
                row.push_str(&string);
                row.push_str(" |\n");
                rows.push(row.clone());
                row_len = row.len();
                row.clear();
                row.push_str("| ");
            } else {
                row.push_str(&string);
                // row.push_str(" ; ");
            }
        }
        f.write_char('\n')?;
        f.write_str(&"_".repeat(row_len - 1))?;
        f.write_char('\n')?;
        for row in rows {
            f.write_str(&row)?;
        }
        f.write_str(&"-".repeat(row_len - 1))?;
        f.write_str(&format!(" modulo F_{:}", self.field_mod))
    }
}

#[cfg(test)]
mod tests {
    use crate::math::finite_field::{FFRep, FiniteField};

    use super::FFMatrix;

    fn basic_matrix() -> FFMatrix {
        let p = 9_u32;
        let entries: Vec<FFRep> = (0..12).collect();
        FFMatrix::new(entries, 3, 4, p)
    }
    #[test]
    fn test_transpose() {
        let mut m = basic_matrix();
        println!("m:{:}", m);
        m.transpose();
        println!("m transposed:{:}", m);
    }

    #[test]
    fn test_rref() {
        let mut m = basic_matrix();
        println!("m:{:}", m);
        m.rref();
        println!("m rref:{:}", m);
    }

    #[test]
    fn test_block_access() {
        let m = basic_matrix();
        let block = m.clone_block((2, 2), (0, 0));
        println!("block:{:?}", block);
    }

    #[test]
    fn test_multiplication() {
        let mat = FFMatrix::new(vec![1, 1, 0, 1, 0, 1], 2, 3, 3);
        let vec = vec![
            FiniteField::new(0, 3),
            FiniteField::new(0, 3),
            FiniteField::new(0, 3),
        ];
        println!("{:?}", &mat * &vec);
    }

    #[test]
    fn test_convert_indices() {
        let mat = FFMatrix::id(4, 7);
        for ix in 0..mat.n_rows {
            for jx in 0..mat.n_cols {
                let ix = mat.convert_indices(ix, jx);
                dbg!(&mat.entries[ix]);
            }
        }
    }

    #[test]
    fn test_rank() {
        let id = FFMatrix::id(100, 2);
        dbg!(id.rank());
    }

    #[test]
    fn serialization() {
        let m = FFMatrix::id(3, 7);
        let s = serde_json::to_string_pretty(&m).unwrap();
        println!("s: {:}", s);
        let m: Result<FFMatrix, _> = serde_json::from_str(&s[..]);
        println!("{:}", m.unwrap());
    }
}
