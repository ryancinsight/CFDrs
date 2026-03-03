//! Shared disjoint-set union (DSU) for arrangement phases.
//!
//! Provides a single union-find implementation with union-by-size and path
//! compression for coplanar grouping and fragment-component labeling.

/// Disjoint-set union over contiguous indices `[0, n)`.
///
/// # Algorithm
///
/// - `find`: path compression (parent-halving).
/// - `union`: union-by-size (attach smaller root under larger root).
///
/// # Theorem — Amortized Near-Constant Operations
///
/// Union-find with path compression and union-by-size has amortized
/// `O(α(n))` complexity per operation, where `α` is the inverse Ackermann
/// function (Tarjan). For all practical mesh sizes this is effectively
/// constant-time. ∎
#[derive(Debug, Clone)]
pub(crate) struct DisjointSet {
    parent: Vec<usize>,
    size: Vec<usize>,
}

impl DisjointSet {
    #[inline]
    pub(crate) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }

    #[inline]
    pub(crate) fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            let p = self.parent[x];
            let gp = self.parent[p];
            self.parent[x] = gp;
            x = gp;
        }
        x
    }

    /// Union the sets containing `a` and `b`.
    ///
    /// Returns `true` iff a merge happened.
    #[inline]
    pub(crate) fn union(&mut self, a: usize, b: usize) -> bool {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return false;
        }
        if self.size[ra] < self.size[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        self.size[ra] += self.size[rb];
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adversarial_chain_connectivity_collapses_to_one_root() {
        let mut dsu = DisjointSet::new(8);
        for i in 0..7 {
            dsu.union(i, i + 1);
        }
        let r0 = dsu.find(0);
        for i in 1..8 {
            assert_eq!(dsu.find(i), r0);
        }
    }

    #[test]
    fn union_order_invariance_of_partition() {
        let edges = [(0usize, 1usize), (1, 2), (4, 5), (2, 3), (5, 6), (6, 7)];

        let mut a = DisjointSet::new(10);
        for &(u, v) in &edges {
            a.union(u, v);
        }

        let mut b = DisjointSet::new(10);
        for &(u, v) in edges.iter().rev() {
            b.union(u, v);
        }

        for i in 0..10 {
            for j in 0..10 {
                let same_a = a.find(i) == a.find(j);
                let same_b = b.find(i) == b.find(j);
                assert_eq!(same_a, same_b);
            }
        }
    }
}
