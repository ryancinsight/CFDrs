//! Element-to-node indices in one contiguous buffer.
//!
//! The spectral meshes and the global assembly each mapped an element to its
//! node indices with a `Vec<usize>` per element: `E + 1` heap blocks for `E`
//! elements, and a pointer hop before any element's nodes could be read.
//! [`ElementNodes`] holds every element's run back to back behind an offset
//! table, so the map is one index allocation plus one offset allocation and
//! assembly walks element nodes in memory order.

/// Element-to-node indices in one contiguous buffer with an offset table.
///
/// Element `e` owns `values[offsets[e]..offsets[e + 1]]`; elements may differ
/// in length, which is why the table exists rather than a stride.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ElementNodes {
    values: Vec<usize>,
    offsets: Vec<usize>,
}

impl ElementNodes {
    /// Assemble from an offset table and its backing buffer.
    ///
    /// # Panics
    ///
    /// Panics when `offsets` is empty, does not start at zero, does not end at
    /// `values.len()`, or is not non-decreasing: each is a constructor
    /// invariant rather than a caller input.
    fn from_parts(values: Vec<usize>, offsets: Vec<usize>) -> Self {
        assert!(
            !offsets.is_empty(),
            "invariant: an element offset table carries its leading zero"
        );
        assert_eq!(
            offsets.first().copied(),
            Some(0),
            "invariant: the first element starts at offset zero"
        );
        assert_eq!(
            offsets.last().copied(),
            Some(values.len()),
            "invariant: the offset table ends at the index buffer's length"
        );
        assert!(
            offsets.windows(2).all(|pair| pair[0] <= pair[1]),
            "invariant: element offsets ascend"
        );
        Self { values, offsets }
    }

    /// Flatten a per-element index list into one buffer.
    pub(super) fn from_elements(elements: Vec<Vec<usize>>) -> Self {
        let total: usize = elements.iter().map(Vec::len).sum();
        let mut values = Vec::with_capacity(total);
        let mut offsets = Vec::with_capacity(elements.len() + 1);
        offsets.push(0);
        for element in elements {
            values.extend_from_slice(&element);
            offsets.push(values.len());
        }
        Self::from_parts(values, offsets)
    }

    /// Build `element_count` equally long runs from a node-index formula.
    ///
    /// `index_of(element, local_node)` is the global node index. Both spectral
    /// meshes share the layout `element * (nodes_per_element - 1) +
    /// local_node`, the shape that shares the boundary node between neighbours.
    pub(super) fn from_uniform(
        element_count: usize,
        nodes_per_element: usize,
        index_of: impl Fn(usize, usize) -> usize,
    ) -> Self {
        let mut values = Vec::with_capacity(element_count * nodes_per_element);
        let mut offsets = Vec::with_capacity(element_count + 1);
        offsets.push(0);
        for element in 0..element_count {
            for local_node in 0..nodes_per_element {
                values.push(index_of(element, local_node));
            }
            offsets.push(values.len());
        }
        Self::from_parts(values, offsets)
    }

    /// The number of elements mapped.
    pub(super) fn element_count(&self) -> usize {
        self.offsets.len() - 1
    }

    /// The node indices of element `element`.
    ///
    /// # Panics
    ///
    /// Panics when `element >= self.element_count()`.
    pub(super) fn element(&self, element: usize) -> &[usize] {
        &self.values[self.offsets[element]..self.offsets[element + 1]]
    }
}

#[cfg(test)]
mod tests {
    use super::ElementNodes;

    #[test]
    fn flattening_preserves_every_element_run() {
        let elements = vec![vec![0, 1, 2], vec![2, 3], vec![], vec![3, 4, 5]];
        let nodes = ElementNodes::from_elements(elements.clone());
        assert_eq!(nodes.element_count(), elements.len());
        for (index, expected) in elements.iter().enumerate() {
            assert_eq!(nodes.element(index), expected.as_slice());
        }
    }

    #[test]
    fn the_uniform_layout_shares_the_boundary_node() {
        let nodes = ElementNodes::from_uniform(3, 3, |element, local| element * 2 + local);
        assert_eq!(nodes.element_count(), 3);
        assert_eq!(nodes.element(0), [0, 1, 2]);
        assert_eq!(nodes.element(1), [2, 3, 4]);
        assert_eq!(nodes.element(2), [4, 5, 6]);
        assert_eq!(nodes.element(1)[0], nodes.element(0)[2]);
    }

    #[test]
    fn an_empty_map_is_one_offset() {
        let nodes = ElementNodes::from_elements(Vec::new());
        assert_eq!(nodes.element_count(), 0);
    }

    #[test]
    fn cloned_maps_share_no_storage() {
        let nodes = ElementNodes::from_uniform(2, 2, |element, local| element + local);
        let clone = nodes.clone();
        assert_eq!(clone, nodes);
        assert_eq!(clone.element(1), nodes.element(1));
    }
}
