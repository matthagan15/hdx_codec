use std::{collections::HashSet, path::Path, rc::Rc};

use crate::math::finite_field::FFRep;

pub trait RankMatrix {
    fn new(field_mod: FFRep) -> Self;
    fn split(&mut self, row_ixs: HashSet<usize>) -> Self;
    /// creates a pivot at the smallest nonzero entry within the given row
    /// returns None if the row is zero
    fn pivotize_row(&mut self, row_ix: usize) -> Option<usize>;
    fn pivotize_row_within_range(
        &mut self,
        row_ix: usize,
        range: impl AsRef<[usize]>,
    ) -> Option<usize>;
    fn eliminate_rows(&mut self, pivot: (usize, usize), rows_to_eliminate: impl AsRef<[usize]>);
    fn eliminate_all_rows(&mut self, pivot: (usize, usize));
    /// Inserts the given entry at (`row_ix`, `col_ix`) = entry + whatever was already there % field_mod
    fn insert(&mut self, row_ix: usize, col_ix: usize, entry: FFRep);
    fn get(&self, row_ix: usize, col_ix: usize) -> FFRep;
    fn to_disk(self) -> Result<(), ()>;
    fn from_disk(path: impl AsRef<Path>) -> Option<Rc<Self>>;
}