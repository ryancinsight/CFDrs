//! Strongly-typed indices for mesh elements.
//!
//! Each element kind (vertex, face, edge, half-edge, region) gets its own newtype
//! so they cannot be accidentally interchanged. All are `Copy`, `Ord`, `Hash`.

use std::fmt;

macro_rules! define_index {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub u32);

        impl $name {
            /// Create a new index from a `u32`.
            #[inline]
            pub const fn new(raw: u32) -> Self {
                Self(raw)
            }

            /// Create from `usize`, panicking if out of range.
            #[inline]
            pub fn from_usize(v: usize) -> Self {
                Self(u32::try_from(v).expect("index overflow"))
            }

            /// Convert to `usize`.
            #[inline]
            pub const fn as_usize(self) -> usize {
                self.0 as usize
            }

            /// Raw `u32` value.
            #[inline]
            pub const fn raw(self) -> u32 {
                self.0
            }

            /// Sentinel value indicating "no element".
            pub const INVALID: Self = Self(u32::MAX);

            /// Check if this is the sentinel value.
            #[inline]
            pub const fn is_valid(self) -> bool {
                self.0 != u32::MAX
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}({})", stringify!($name), self.0)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}", self.0)
            }
        }

        impl From<u32> for $name {
            #[inline]
            fn from(v: u32) -> Self { Self(v) }
        }

        impl From<usize> for $name {
            #[inline]
            fn from(v: usize) -> Self { Self::from_usize(v) }
        }
    };
}

define_index! {
    /// Index into the vertex pool.
    VertexId
}

define_index! {
    /// Index into the face store.
    FaceId
}

define_index! {
    /// Index into the edge store.
    EdgeId
}

define_index! {
    /// Index into the half-edge store.
    HalfEdgeId
}

define_index! {
    /// Region identifier (channel, junction, boundary, etc.).
    RegionId
}
