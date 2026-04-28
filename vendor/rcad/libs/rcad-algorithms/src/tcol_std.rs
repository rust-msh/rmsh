//! TColStd-style collection types.
//!
//! Provides standard collection types analogous to OCCT's TColStd package.
//! These are 1-based indexed arrays, sequences, maps, and lists for compatibility
//! with OCCT algorithms.

use std::collections::{HashMap, HashSet};
use std::hash::Hash;

// ============================================================================
// Array Types
// ============================================================================

/// A generic 1-dimensional array with 1-based indexing (OCCT-style).
///
/// In OCCT, arrays use a lower and upper bound, with indices ranging from
/// `lower` to `upper` inclusive. This implementation mirrors that behavior.
#[derive(Debug, Clone, PartialEq)]
pub struct Array1<T> {
    data: Vec<T>,
    lower: i32,
}

impl<T: Clone> Array1<T> {
    /// Creates a new array with elements initialized to their default.
    ///
    /// # Arguments
    /// * `lower` - Lower bound index (inclusive)
    /// * `upper` - Upper bound index (inclusive)
    ///
    /// # Panics
    /// Panics if `upper < lower`.
    pub fn new(lower: i32, upper: i32) -> Self
    where
        T: Default,
    {
        if upper < lower {
            panic!("upper bound ({}) must be >= lower bound ({})", upper, lower);
        }
        let len = (upper - lower + 1) as usize;
        Self {
            data: vec![T::default(); len],
            lower,
        }
    }

    /// Creates a new array from a Vec, with the specified lower bound.
    ///
    /// The upper bound is computed as `lower + data.len() - 1`.
    pub fn from_vec(data: Vec<T>, lower: i32) -> Self {
        let lower = if data.is_empty() { 1 } else { lower };
        Self { data, lower }
    }

    /// Returns the element at the given index.
    ///
    /// Returns `None` if the index is out of bounds.
    pub fn get(&self, index: i32) -> Option<&T> {
        if index < self.lower {
            return None;
        }
        self.data.get((index - self.lower) as usize)
    }

    /// Returns a mutable reference to the element at the given index.
    ///
    /// Returns `None` if the index is out of bounds.
    pub fn get_mut(&mut self, index: i32) -> Option<&mut T> {
        if index < self.lower {
            return None;
        }
        self.data.get_mut((index - self.lower) as usize)
    }

    /// Sets the element at the given index.
    ///
    /// Does nothing if the index is out of bounds.
    pub fn set(&mut self, index: i32, value: T) {
        if let Some(elem) = self.get_mut(index) {
            *elem = value;
        }
    }

    /// Returns the number of elements in the array.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the array is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns the lower bound of the array.
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// Returns the upper bound of the array.
    pub fn upper(&self) -> i32 {
        self.lower + self.data.len() as i32 - 1
    }

    /// Returns an iterator over the elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Returns a mutable iterator over the elements.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }

    /// Initializes all elements to the given value.
    pub fn init(&mut self, value: T) {
        for elem in &mut self.data {
            *elem = value.clone();
        }
    }
}

/// A 1-dimensional array of real numbers (f64).
pub type Array1OfReal = Array1<f64>;

/// A 1-dimensional array of integers (i32).
pub type Array1OfInteger = Array1<i32>;

impl Array1OfReal {
    /// Creates a new array of reals with all elements initialized to 0.0.
    pub fn new_real(lower: i32, upper: i32) -> Self {
        Self::new(lower, upper)
    }

    /// Creates a new array from a Vec of f64, with lower bound 1.
    pub fn from_vec_real(data: Vec<f64>) -> Self {
        Self::from_vec(data, 1)
    }
}

impl Array1OfInteger {
    /// Creates a new array of integers with all elements initialized to 0.
    pub fn new_integer(lower: i32, upper: i32) -> Self {
        Self::new(lower, upper)
    }

    /// Creates a new array from a Vec of i32, with lower bound 1.
    pub fn from_vec_integer(data: Vec<i32>) -> Self {
        Self::from_vec(data, 1)
    }
}

// ============================================================================
// 2D Array Types
// ============================================================================

/// A generic 2-dimensional array with 1-based indexing (OCCT-style).
#[derive(Debug, Clone, PartialEq)]
pub struct Array2<T> {
    data: Vec<Vec<T>>,
    lower_row: i32,
    upper_row: i32,
    lower_col: i32,
    upper_col: i32,
}

impl<T: Clone> Array2<T> {
    /// Creates a new 2D array with elements initialized to their default.
    ///
    /// # Arguments
    /// * `lower_row` - Lower row bound (inclusive)
    /// * `upper_row` - Upper row bound (inclusive)
    /// * `lower_col` - Lower column bound (inclusive)
    /// * `upper_col` - Upper column bound (inclusive)
    ///
    /// # Panics
    /// Panics if `upper_row < lower_row` or `upper_col < lower_col`.
    pub fn new(lower_row: i32, upper_row: i32, lower_col: i32, upper_col: i32) -> Self
    where
        T: Default,
    {
        if upper_row < lower_row {
            panic!("upper_row ({}) must be >= lower_row ({})", upper_row, lower_row);
        }
        if upper_col < lower_col {
            panic!("upper_col ({}) must be >= lower_col ({})", upper_col, lower_col);
        }
        let rows = (upper_row - lower_row + 1) as usize;
        let cols = (upper_col - lower_col + 1) as usize;
        Self {
            data: vec![vec![T::default(); cols]; rows],
            lower_row,
            upper_row,
            lower_col,
            upper_col,
        }
    }

    /// Returns the element at the given (row, col) index.
    pub fn get(&self, row: i32, col: i32) -> Option<&T> {
        if row < self.lower_row || row > self.upper_row {
            return None;
        }
        if col < self.lower_col || col > self.upper_col {
            return None;
        }
        self.data
            .get((row - self.lower_row) as usize)
            .and_then(|r| r.get((col - self.lower_col) as usize))
    }

    /// Returns a mutable reference to the element at the given (row, col) index.
    pub fn get_mut(&mut self, row: i32, col: i32) -> Option<&mut T> {
        if row < self.lower_row || row > self.upper_row {
            return None;
        }
        if col < self.lower_col || col > self.upper_col {
            return None;
        }
        self.data
            .get_mut((row - self.lower_row) as usize)
            .and_then(|r| r.get_mut((col - self.lower_col) as usize))
    }

    /// Sets the element at the given (row, col) index.
    pub fn set(&mut self, row: i32, col: i32, value: T) {
        if let Some(elem) = self.get_mut(row, col) {
            *elem = value;
        }
    }

    /// Returns the number of rows.
    pub fn row_count(&self) -> usize {
        (self.upper_row - self.lower_row + 1) as usize
    }

    /// Returns the number of columns.
    pub fn col_count(&self) -> usize {
        (self.upper_col - self.lower_col + 1) as usize
    }

    /// Returns the total number of elements.
    pub fn len(&self) -> usize {
        self.row_count() * self.col_count()
    }

    /// Returns true if the array is empty.
    pub fn is_empty(&self) -> bool {
        self.row_count() == 0 || self.col_count() == 0
    }

    /// Returns the lower row bound.
    pub fn lower_row(&self) -> i32 {
        self.lower_row
    }

    /// Returns the upper row bound.
    pub fn upper_row(&self) -> i32 {
        self.upper_row
    }

    /// Returns the lower column bound.
    pub fn lower_col(&self) -> i32 {
        self.lower_col
    }

    /// Returns the upper column bound.
    pub fn upper_col(&self) -> i32 {
        self.upper_col
    }
}

/// A 2-dimensional array of real numbers.
pub type Array2OfReal = Array2<f64>;

/// A 2-dimensional array of integers.
pub type Array2OfInteger = Array2<i32>;

// ============================================================================
// Sequence Types
// ============================================================================

/// A sequence of elements (doubly-linked list semantics).
///
/// Analogous to OCCT's `TColStd_SequenceOfReal` and `TColStd_SequenceOfInteger`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Sequence<T> {
    data: Vec<T>,
}

impl<T> Sequence<T> {
    /// Creates a new empty sequence.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Returns the number of elements in the sequence.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Appends an element to the end of the sequence.
    pub fn append(&mut self, value: T) {
        self.data.push(value);
    }

    /// Prepends an element to the beginning of the sequence.
    pub fn prepend(&mut self, value: T) {
        self.data.insert(0, value);
    }

    /// Inserts an element at the given 1-based index.
    ///
    /// Does nothing if the index is out of bounds (must be 1 to len+1).
    pub fn insert(&mut self, index: usize, value: T) {
        if index == 0 || index > self.data.len() + 1 {
            return;
        }
        self.data.insert(index - 1, value);
    }

    /// Removes the element at the given 1-based index.
    ///
    /// Does nothing if the index is out of bounds.
    pub fn remove(&mut self, index: usize) {
        if index == 0 || index > self.data.len() {
            return;
        }
        self.data.remove(index - 1);
    }

    /// Returns a reference to the first element.
    pub fn first(&self) -> Option<&T> {
        self.data.first()
    }

    /// Returns a reference to the last element.
    pub fn last(&self) -> Option<&T> {
        self.data.last()
    }

    /// Returns a reference to the element at the given 1-based index.
    pub fn get(&self, index: usize) -> Option<&T> {
        if index == 0 {
            return None;
        }
        self.data.get(index - 1)
    }

    /// Returns a mutable reference to the element at the given 1-based index.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut T> {
        if index == 0 {
            return None;
        }
        self.data.get_mut(index - 1)
    }

    /// Sets the element at the given 1-based index.
    pub fn set(&mut self, index: usize, value: T) {
        if let Some(elem) = self.get_mut(index) {
            *elem = value;
        }
    }

    /// Returns an iterator over the elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Clears all elements from the sequence.
    pub fn clear(&mut self) {
        self.data.clear();
    }

    /// Reverses the sequence in place.
    pub fn reverse(&mut self) {
        self.data.reverse();
    }

    /// Returns the underlying Vec.
    pub fn into_vec(self) -> Vec<T> {
        self.data
    }
}

/// A sequence of real numbers.
pub type SequenceOfReal = Sequence<f64>;

/// A sequence of integers.
pub type SequenceOfInteger = Sequence<i32>;

// ============================================================================
// Map Types
// ============================================================================

/// A set of unique keys (analogous to OCCT's `TColStd_MapOfInteger`).
#[derive(Debug, Clone, Default)]
pub struct MapOfInteger {
    set: HashSet<i32>,
}

impl MapOfInteger {
    /// Creates a new empty map.
    pub fn new() -> Self {
        Self {
            set: HashSet::new(),
        }
    }

    /// Creates a new empty map with the given capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            set: HashSet::with_capacity(capacity),
        }
    }

    /// Inserts a key into the map.
    ///
    /// Returns `true` if the key was newly inserted, `false` if it already existed.
    pub fn insert(&mut self, key: i32) -> bool {
        self.set.insert(key)
    }

    /// Returns `true` if the map contains the key.
    pub fn contains(&self, key: i32) -> bool {
        self.set.contains(&key)
    }

    /// Removes a key from the map.
    ///
    /// Returns `true` if the key was present, `false` otherwise.
    pub fn remove(&mut self, key: i32) -> bool {
        self.set.remove(&key)
    }

    /// Returns the number of elements in the map.
    pub fn len(&self) -> usize {
        self.set.len()
    }

    /// Returns true if the map is empty.
    pub fn is_empty(&self) -> bool {
        self.set.is_empty()
    }

    /// Clears all elements from the map.
    pub fn clear(&mut self) {
        self.set.clear();
    }

    /// Returns an iterator over the keys.
    pub fn iter(&self) -> impl Iterator<Item = i32> + '_ {
        self.set.iter().copied()
    }
}

/// A set of unique real numbers.
///
/// Note: Uses direct comparison of f64 values. For approximate comparison,
/// consider using a different approach.
#[derive(Debug, Clone, Default)]
pub struct MapOfReal {
    set: HashSet<OrderedFloat>,
}

/// Wrapper for f64 that implements Hash and Eq for use in HashSet/HashMap.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderedFloat(f64);

impl Eq for OrderedFloat {}

impl Hash for OrderedFloat {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // Use bit representation for hashing
        self.0.to_bits().hash(state);
    }
}

impl MapOfReal {
    /// Creates a new empty map.
    pub fn new() -> Self {
        Self {
            set: HashSet::new(),
        }
    }

    /// Creates a new empty map with the given capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            set: HashSet::with_capacity(capacity),
        }
    }

    /// Inserts a key into the map.
    ///
    /// Returns `true` if the key was newly inserted, `false` if it already existed.
    pub fn insert(&mut self, key: f64) -> bool {
        self.set.insert(OrderedFloat(key))
    }

    /// Returns `true` if the map contains the key.
    pub fn contains(&self, key: f64) -> bool {
        self.set.contains(&OrderedFloat(key))
    }

    /// Removes a key from the map.
    ///
    /// Returns `true` if the key was present, `false` otherwise.
    pub fn remove(&mut self, key: f64) -> bool {
        self.set.remove(&OrderedFloat(key))
    }

    /// Returns the number of elements in the map.
    pub fn len(&self) -> usize {
        self.set.len()
    }

    /// Returns true if the map is empty.
    pub fn is_empty(&self) -> bool {
        self.set.is_empty()
    }

    /// Clears all elements from the map.
    pub fn clear(&mut self) {
        self.set.clear();
    }

    /// Returns an iterator over the keys.
    pub fn iter(&self) -> impl Iterator<Item = f64> + '_ {
        self.set.iter().map(|f| f.0)
    }
}

// ============================================================================
// Indexed Map Types
// ============================================================================

/// A map that maintains insertion order and allows indexed access.
///
/// Analogous to OCCT's `TColStd_IndexedMapOfInteger`.
#[derive(Debug, Clone, Default)]
pub struct IndexedMap<K> {
    keys: Vec<K>,
    indices: HashMap<K, usize>,
}

impl<K: Clone + Eq + Hash> IndexedMap<K> {
    /// Creates a new empty indexed map.
    pub fn new() -> Self {
        Self {
            keys: Vec::new(),
            indices: HashMap::new(),
        }
    }

    /// Adds a key to the map.
    ///
    /// Returns the index of the key (1-based). If the key already exists,
    /// returns its existing index.
    pub fn add(&mut self, key: K) -> usize {
        if let Some(&index) = self.indices.get(&key) {
            return index + 1; // 1-based
        }
        let index = self.keys.len();
        self.keys.push(key.clone());
        self.indices.insert(key, index);
        index + 1 // 1-based
    }

    /// Returns true if the map contains the key.
    pub fn contains(&self, key: &K) -> bool {
        self.indices.contains_key(key)
    }

    /// Returns the 1-based index of the key, or 0 if not found.
    pub fn find_index(&self, key: &K) -> usize {
        self.indices.get(key).map(|&i| i + 1).unwrap_or(0)
    }

    /// Returns the key at the given 1-based index.
    pub fn find_key(&self, index: usize) -> Option<&K> {
        if index == 0 {
            return None;
        }
        self.keys.get(index - 1)
    }

    /// Removes a key from the map.
    ///
    /// Note: This operation is O(n) because it needs to rebuild the index map.
    pub fn remove(&mut self, key: &K) -> bool {
        if let Some(index) = self.indices.remove(key) {
            self.keys.remove(index);
            // Rebuild indices
            self.indices.clear();
            for (i, k) in self.keys.iter().enumerate() {
                self.indices.insert(k.clone(), i);
            }
            true
        } else {
            false
        }
    }

    /// Returns the number of elements in the map.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Returns true if the map is empty.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Clears all elements from the map.
    pub fn clear(&mut self) {
        self.keys.clear();
        self.indices.clear();
    }

    /// Returns an iterator over the keys in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = &K> {
        self.keys.iter()
    }
}

/// An indexed map of integers.
pub type IndexedMapOfInteger = IndexedMap<i32>;

/// An indexed map of real numbers.
pub type IndexedMapOfReal = IndexedMap<OrderedFloat>;

// ============================================================================
// List Types
// ============================================================================

/// A doubly-linked list of elements.
///
/// Analogous to OCCT's `TColStd_ListOfReal` and `TColStd_ListOfInteger`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct List<T> {
    data: Vec<T>,
}

impl<T> List<T> {
    /// Creates a new empty list.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Returns the number of elements in the list.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the list is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Appends an element to the end of the list.
    pub fn append(&mut self, value: T) {
        self.data.push(value);
    }

    /// Prepends an element to the beginning of the list.
    pub fn prepend(&mut self, value: T) {
        self.data.insert(0, value);
    }

    /// Removes and returns the first element.
    pub fn pop_front(&mut self) -> Option<T> {
        if self.data.is_empty() {
            None
        } else {
            Some(self.data.remove(0))
        }
    }

    /// Removes and returns the last element.
    pub fn pop_back(&mut self) -> Option<T> {
        self.data.pop()
    }

    /// Returns a reference to the first element.
    pub fn first(&self) -> Option<&T> {
        self.data.first()
    }

    /// Returns a reference to the last element.
    pub fn last(&self) -> Option<&T> {
        self.data.last()
    }

    /// Returns an iterator over the elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Clears all elements from the list.
    pub fn clear(&mut self) {
        self.data.clear();
    }

    /// Reverses the list in place.
    pub fn reverse(&mut self) {
        self.data.reverse();
    }

    /// Returns the underlying Vec.
    pub fn into_vec(self) -> Vec<T> {
        self.data
    }
}

/// A list of real numbers.
pub type ListOfReal = List<f64>;

/// A list of integers.
pub type ListOfInteger = List<i32>;

// ============================================================================
// Data Map Types (Key-Value Maps)
// ============================================================================

/// A hash map from integer keys to values.
///
/// Analogous to OCCT's `TColStd_DataMapOfIntegerReal`.
pub type DataMapOfIntegerReal = HashMap<i32, f64>;

/// A hash map from integer keys to integers.
pub type DataMapOfIntegerInteger = HashMap<i32, i32>;

// ============================================================================
// Conversion Utilities
// ============================================================================

/// Converts a Vec to an Array1 with lower bound 1.
pub fn vec_to_array1<T: Clone>(vec: Vec<T>) -> Array1<T> {
    Array1::from_vec(vec, 1)
}

/// Converts a Vec to an Array1 with a custom lower bound.
pub fn vec_to_array1_with_lower<T: Clone>(vec: Vec<T>, lower: i32) -> Array1<T> {
    Array1::from_vec(vec, lower)
}

/// Converts an Array1 to a Vec.
pub fn array1_to_vec<T: Clone>(array: &Array1<T>) -> Vec<T> {
    array.data.clone()
}

// ============================================================================
// HArray1 Types (Handle to Array1 - Reference Counted)
// ============================================================================

use std::rc::Rc;

/// A reference-counted handle to an Array1.
///
/// Analogous to OCCT's `TColStd_HArray1OfReal`.
#[derive(Debug, Clone)]
pub struct HArray1<T>(Rc<Array1<T>>);

impl<T: Clone> HArray1<T> {
    /// Creates a new handle to an Array1.
    pub fn new(lower: i32, upper: i32) -> Self
    where
        T: Default,
    {
        Self(Rc::new(Array1::new(lower, upper)))
    }

    /// Creates a handle from an existing Array1.
    pub fn from_array(array: Array1<T>) -> Self {
        Self(Rc::new(array))
    }

    /// Returns a reference to the underlying Array1.
    pub fn array(&self) -> &Array1<T> {
        &self.0
    }

    /// Returns the lower bound.
    pub fn lower(&self) -> i32 {
        self.0.lower()
    }

    /// Returns the upper bound.
    pub fn upper(&self) -> i32 {
        self.0.upper()
    }

    /// Returns the length.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns true if the array is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Gets an element by index.
    pub fn get(&self, index: i32) -> Option<&T> {
        self.0.get(index)
    }

    /// Returns an iterator over elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.0.iter()
    }
}

/// A reference-counted handle to an Array1OfReal.
pub type HArray1OfReal = HArray1<f64>;

/// A reference-counted handle to an Array1OfInteger.
pub type HArray1OfInteger = HArray1<i32>;

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Array1 Tests
    #[test]
    fn array1_new_creates_correct_size() {
        let arr: Array1OfReal = Array1::new(1, 5);
        assert_eq!(arr.len(), 5);
        assert_eq!(arr.lower(), 1);
        assert_eq!(arr.upper(), 5);
    }

    #[test]
    fn array1_get_set_works() {
        let mut arr: Array1OfReal = Array1::new(1, 3);
        arr.set(1, 10.0);
        arr.set(2, 20.0);
        arr.set(3, 30.0);

        assert_eq!(arr.get(1), Some(&10.0));
        assert_eq!(arr.get(2), Some(&20.0));
        assert_eq!(arr.get(3), Some(&30.0));
        assert_eq!(arr.get(0), None);
        assert_eq!(arr.get(4), None);
    }

    #[test]
    fn array1_from_vec_works() {
        let arr = Array1OfReal::from_vec_real(vec![1.0, 2.0, 3.0]);
        assert_eq!(arr.len(), 3);
        assert_eq!(arr.lower(), 1);
        assert_eq!(arr.get(1), Some(&1.0));
        assert_eq!(arr.get(2), Some(&2.0));
        assert_eq!(arr.get(3), Some(&3.0));
    }

    #[test]
    fn array1_init_sets_all_elements() {
        let mut arr: Array1OfInteger = Array1::new(1, 5);
        arr.init(42);
        for i in 1..=5 {
            assert_eq!(arr.get(i), Some(&42));
        }
    }

    #[test]
    fn array1_iter_works() {
        let arr = Array1OfReal::from_vec_real(vec![1.0, 2.0, 3.0]);
        let sum: f64 = arr.iter().sum();
        assert!((sum - 6.0).abs() < 1e-9);
    }

    #[test]
    #[should_panic]
    fn array1_invalid_bounds_panics() {
        let _: Array1OfReal = Array1::new(5, 1);
    }

    // Array2 Tests
    #[test]
    fn array2_new_creates_correct_size() {
        let arr: Array2OfReal = Array2::new(1, 3, 1, 4);
        assert_eq!(arr.row_count(), 3);
        assert_eq!(arr.col_count(), 4);
        assert_eq!(arr.len(), 12);
    }

    #[test]
    fn array2_get_set_works() {
        let mut arr: Array2OfReal = Array2::new(1, 2, 1, 2);
        arr.set(1, 1, 1.0);
        arr.set(1, 2, 2.0);
        arr.set(2, 1, 3.0);
        arr.set(2, 2, 4.0);

        assert_eq!(arr.get(1, 1), Some(&1.0));
        assert_eq!(arr.get(1, 2), Some(&2.0));
        assert_eq!(arr.get(2, 1), Some(&3.0));
        assert_eq!(arr.get(2, 2), Some(&4.0));
        assert_eq!(arr.get(0, 1), None);
        assert_eq!(arr.get(3, 1), None);
    }

    // Sequence Tests
    #[test]
    fn sequence_append_works() {
        let mut seq = SequenceOfReal::new();
        seq.append(1.0);
        seq.append(2.0);
        seq.append(3.0);

        assert_eq!(seq.len(), 3);
        assert_eq!(seq.first(), Some(&1.0));
        assert_eq!(seq.last(), Some(&3.0));
    }

    #[test]
    fn sequence_prepend_works() {
        let mut seq = SequenceOfReal::new();
        seq.append(2.0);
        seq.prepend(1.0);

        assert_eq!(seq.len(), 2);
        assert_eq!(seq.first(), Some(&1.0));
        assert_eq!(seq.last(), Some(&2.0));
    }

    #[test]
    fn sequence_insert_works() {
        let mut seq = SequenceOfReal::new();
        seq.append(1.0);
        seq.append(3.0);
        seq.insert(2, 2.0); // Insert at 1-based index 2

        assert_eq!(seq.len(), 3);
        assert_eq!(seq.get(1), Some(&1.0));
        assert_eq!(seq.get(2), Some(&2.0));
        assert_eq!(seq.get(3), Some(&3.0));
    }

    #[test]
    fn sequence_remove_works() {
        let mut seq = SequenceOfReal::new();
        seq.append(1.0);
        seq.append(2.0);
        seq.append(3.0);
        seq.remove(2); // Remove at 1-based index 2

        assert_eq!(seq.len(), 2);
        assert_eq!(seq.get(1), Some(&1.0));
        assert_eq!(seq.get(2), Some(&3.0));
    }

    #[test]
    fn sequence_reverse_works() {
        let mut seq = SequenceOfReal::new();
        seq.append(1.0);
        seq.append(2.0);
        seq.append(3.0);
        seq.reverse();

        assert_eq!(seq.first(), Some(&3.0));
        assert_eq!(seq.last(), Some(&1.0));
    }

    // MapOfInteger Tests
    #[test]
    fn map_of_integer_basic_operations() {
        let mut map = MapOfInteger::new();
        assert!(map.is_empty());

        assert!(map.insert(1));
        assert!(map.insert(2));
        assert!(!map.insert(1)); // Already exists

        assert_eq!(map.len(), 2);
        assert!(map.contains(1));
        assert!(map.contains(2));
        assert!(!map.contains(3));

        assert!(map.remove(1));
        assert!(!map.contains(1));
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn map_of_integer_clear() {
        let mut map = MapOfInteger::new();
        map.insert(1);
        map.insert(2);
        map.insert(3);
        map.clear();

        assert!(map.is_empty());
        assert_eq!(map.len(), 0);
    }

    // MapOfReal Tests
    #[test]
    fn map_of_real_basic_operations() {
        let mut map = MapOfReal::new();
        assert!(map.is_empty());

        assert!(map.insert(1.0));
        assert!(map.insert(2.0));
        assert!(!map.insert(1.0)); // Already exists

        assert_eq!(map.len(), 2);
        assert!(map.contains(1.0));
        assert!(map.contains(2.0));
        assert!(!map.contains(3.0));

        assert!(map.remove(1.0));
        assert!(!map.contains(1.0));
        assert_eq!(map.len(), 1);
    }

    // IndexedMap Tests
    #[test]
    fn indexed_map_add_and_find() {
        let mut map = IndexedMapOfInteger::new();

        let idx1 = map.add(10);
        let idx2 = map.add(20);
        let idx3 = map.add(10); // Duplicate, should return same index

        assert_eq!(idx1, 1);
        assert_eq!(idx2, 2);
        assert_eq!(idx3, 1); // Duplicate returns existing index

        assert_eq!(map.len(), 2);
        assert!(map.contains(&10));
        assert!(map.contains(&20));

        assert_eq!(map.find_index(&10), 1);
        assert_eq!(map.find_index(&20), 2);
        assert_eq!(map.find_index(&30), 0); // Not found

        assert_eq!(map.find_key(1), Some(&10));
        assert_eq!(map.find_key(2), Some(&20));
        assert_eq!(map.find_key(3), None);
    }

    #[test]
    fn indexed_map_iter_preserves_order() {
        let mut map = IndexedMapOfInteger::new();
        map.add(1);
        map.add(2);
        map.add(3);

        let values: Vec<i32> = map.iter().copied().collect();
        assert_eq!(values, vec![1, 2, 3]);
    }

    // List Tests
    #[test]
    fn list_basic_operations() {
        let mut list = ListOfReal::new();
        assert!(list.is_empty());

        list.append(1.0);
        list.append(2.0);
        list.prepend(0.0);

        assert_eq!(list.len(), 3);
        assert_eq!(list.first(), Some(&0.0));
        assert_eq!(list.last(), Some(&2.0));

        assert_eq!(list.pop_front(), Some(0.0));
        assert_eq!(list.len(), 2);
        assert_eq!(list.pop_back(), Some(2.0));
        assert_eq!(list.len(), 1);
    }

    #[test]
    fn list_reverse() {
        let mut list = ListOfReal::new();
        list.append(1.0);
        list.append(2.0);
        list.append(3.0);
        list.reverse();

        assert_eq!(list.first(), Some(&3.0));
        assert_eq!(list.last(), Some(&1.0));
    }

    // Conversion Utilities Tests
    #[test]
    fn vec_to_array1_conversion() {
        let vec = vec![1.0, 2.0, 3.0];
        let arr = vec_to_array1(vec);

        assert_eq!(arr.len(), 3);
        assert_eq!(arr.lower(), 1);
        assert_eq!(arr.get(1), Some(&1.0));
    }

    #[test]
    fn vec_to_array1_with_custom_lower() {
        let vec = vec![1.0, 2.0, 3.0];
        let arr = vec_to_array1_with_lower(vec, 0);

        assert_eq!(arr.len(), 3);
        assert_eq!(arr.lower(), 0);
        assert_eq!(arr.upper(), 2);
        assert_eq!(arr.get(0), Some(&1.0));
    }

    #[test]
    fn array1_to_vec_conversion() {
        let arr = Array1OfReal::from_vec_real(vec![1.0, 2.0, 3.0]);
        let vec = array1_to_vec(&arr);

        assert_eq!(vec, vec![1.0, 2.0, 3.0]);
    }

    // HArray1 Tests
    #[test]
    fn harray1_basic_operations() {
        let harr = HArray1OfReal::new(1, 3);
        assert_eq!(harr.len(), 3);
        assert_eq!(harr.lower(), 1);
        assert_eq!(harr.upper(), 3);
    }

    #[test]
    fn harray1_cloning_shares_data() {
        let harr1 = HArray1OfReal::from_array(Array1OfReal::from_vec_real(vec![1.0, 2.0, 3.0]));
        let harr2 = harr1.clone();

        // Both point to the same underlying data
        assert_eq!(harr1.get(1), harr2.get(1));
    }

    // DataMap Tests
    #[test]
    fn data_map_of_integer_real() {
        let mut map = DataMapOfIntegerReal::new();
        map.insert(1, 1.5);
        map.insert(2, 2.5);

        assert_eq!(map.get(&1), Some(&1.5));
        assert_eq!(map.get(&2), Some(&2.5));
        assert_eq!(map.get(&3), None);
    }

    #[test]
    fn data_map_of_integer_integer() {
        let mut map = DataMapOfIntegerInteger::new();
        map.insert(1, 100);
        map.insert(2, 200);

        assert_eq!(map.get(&1), Some(&100));
        assert_eq!(map.get(&2), Some(&200));
    }
}
