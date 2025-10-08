//! Windowed operations for iterators

use std::collections::VecDeque;

/// Iterator adapter for windowed operations
pub struct WindowIterator<I, T> {
    iter: I,
    window: VecDeque<T>,
    size: usize,
}

impl<I, T> WindowIterator<I, T>
where
    I: Iterator<Item = T>,
    T: Clone,
{
    /// Create a new windowed iterator
    pub fn new(iter: I, size: usize) -> Self {
        Self {
            iter,
            window: VecDeque::with_capacity(size),
            size,
        }
    }
}

impl<I, T> Iterator for WindowIterator<I, T>
where
    I: Iterator<Item = T>,
    T: Clone,
{
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        // Fill initial window
        while self.window.len() < self.size {
            match self.iter.next() {
                Some(val) => self.window.push_back(val),
                None => {
                    if self.window.is_empty() {
                        return None;
                    }
                    break;
                }
            }
        }

        if self.window.len() == self.size {
            let result = self.window.iter().cloned().collect();

            // Slide window
            if let Some(next) = self.iter.next() {
                self.window.pop_front();
                self.window.push_back(next);
            } else {
                self.window.clear(); // Mark as exhausted
            }

            Some(result)
        } else {
            None
        }
    }
}

/// Strided window iterator for efficient access patterns
pub struct StridedWindowIterator<I, T> {
    iter: I,
    window_size: usize,
    stride: usize,
    buffer: Vec<T>,
}

impl<I, T> StridedWindowIterator<I, T>
where
    I: Iterator<Item = T>,
    T: Clone,
{
    /// Create a new strided window iterator
    pub fn new(iter: I, window_size: usize, stride: usize) -> Self {
        Self {
            iter,
            window_size,
            stride,
            buffer: Vec::new(),
        }
    }
}

impl<I, T> Iterator for StridedWindowIterator<I, T>
where
    I: Iterator<Item = T>,
    T: Clone,
{
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        // Initial fill
        if self.buffer.is_empty() {
            self.buffer
                .extend(self.iter.by_ref().take(self.window_size));
            if self.buffer.len() < self.window_size {
                return None;
            }
            return Some(self.buffer.clone());
        }

        // Advance by stride
        for _ in 0..self.stride {
            if let Some(val) = self.iter.next() {
                self.buffer.push(val);
                self.buffer.remove(0);
            } else {
                return None;
            }
        }

        Some(self.buffer.clone())
    }
}
