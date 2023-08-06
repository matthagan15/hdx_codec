mod classical;
mod quantum;


mod Pauli {
    use std::ops::Mul;

    #[derive(Debug, Clone, Copy)]
    enum Pauli {
        I = 0,
        X = 1,
        Y = 2,
        Z = 3,
    }

    enum Phase {
        zero,
    }

    impl Mul<&Pauli> for Pauli {
        type Output = Self;

        fn mul(self, rhs: &Self) -> Self::Output {
            match (self, rhs) {
                (Pauli::I, Pauli::I) => Pauli::I,
                (Pauli::I, Pauli::X) => Pauli::X,
                (Pauli::I, Pauli::Y) => Pauli::Y,
                (Pauli::I, Pauli::Z) => Pauli::Z,
                (Pauli::X, Pauli::I) => Pauli::X,
                (Pauli::X, Pauli::X) => Pauli::I,
                (Pauli::X, Pauli::Y) => Pauli::X,
                (Pauli::X, Pauli::Z) => Pauli::X,
                (Pauli::Y, Pauli::I) => Pauli::X,
                (Pauli::Y, Pauli::X) => Pauli::X,
                (Pauli::Y, Pauli::Y) => Pauli::X,
                (Pauli::Y, Pauli::Z) => Pauli::X,
                (Pauli::Z, Pauli::I) => Pauli::X,
                (Pauli::Z, Pauli::X) => Pauli::X,
                (Pauli::Z, Pauli::Y) => Pauli::X,
                (Pauli::Z, Pauli::Z) => Pauli::X,
            }
        }
    }
    
    struct PauliString {
        data: Vec<Pauli>,
    }

    impl PauliString {
        fn mul(&self, rhs: &Self) -> Self {
            self.data.iter().zip(rhs.data.iter()).map(|(a,b)| *a * b);
            todo!()
        }
    }
}



#[cfg(test)]
mod tests {
}
