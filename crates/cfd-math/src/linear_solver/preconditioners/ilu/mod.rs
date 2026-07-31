//! Incomplete LU factorization preconditioner
//!
//! Modular implementation of ILU(k) factorization following SOLID/CUPID principles.
//! Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). SIAM, §10.4.

mod csr_search;
mod ilu0;
mod iluk;
mod triangular;
mod types;

pub use types::IncompleteLU;

#[cfg(test)]
mod tests;
