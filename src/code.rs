use crate::math::finite_field::FiniteField as FF;

pub trait Code {
    /// Takes ownership of a logical message and outputs the physical bits.
    fn encode(&self, message: &Vec<FF>) -> Vec<FF>;

    /// Takes in a logical string and returns the closest codeword,
    /// ties broken arbitrarily.
    fn decode(&self, encrypted: &Vec<FF>) -> Vec<FF>;

    /// Returns `true` if `word` is in the codespace, `false` otherwise.
    fn code_check(&self, word: &Vec<FF>) -> bool;

    /// Returns the parity check of the provided message
    fn parity_check(&self, message: &Vec<FF>) -> Vec<FF>;
}