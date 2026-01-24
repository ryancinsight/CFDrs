use cfd_math::iterators::StridedWindowIterator;

#[test]
fn test_strided_window_basic() {
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let iter = StridedWindowIterator::new(data.into_iter(), 3, 2);
    let result: Vec<Vec<i32>> = iter.collect();

    let expected = vec![
        vec![1, 2, 3],
        vec![3, 4, 5],
        vec![5, 6, 7],
        vec![7, 8, 9],
    ];

    assert_eq!(result, expected);
}

#[test]
fn test_strided_window_overlap() {
    let data = vec![1, 2, 3, 4, 5];
    // Window 3, stride 1 (standard sliding window)
    let iter = StridedWindowIterator::new(data.into_iter(), 3, 1);
    let result: Vec<Vec<i32>> = iter.collect();

    let expected = vec![
        vec![1, 2, 3],
        vec![2, 3, 4],
        vec![3, 4, 5],
    ];

    assert_eq!(result, expected);
}

#[test]
fn test_strided_window_gap() {
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    // Window 2, stride 3
    let iter = StridedWindowIterator::new(data.into_iter(), 2, 3);
    let result: Vec<Vec<i32>> = iter.collect();

    let expected = vec![
        vec![1, 2],
        vec![4, 5],
        vec![7, 8],
    ];

    assert_eq!(result, expected);
}
