use std::marker::PhantomData;

#[derive(Debug)]
pub struct Mask<T> {
    selected: Vec<bool>,
    _indexes: PhantomData<T>,
}

impl<T> Mask<T> {
    pub fn from_bools(selected: Vec<bool>) -> Self {
        Self { selected, _indexes: PhantomData }
    }

    pub fn from_indices(len: usize, indices: impl IntoIterator<Item = usize>) -> Self {
        let mut selected = vec![false; len];
        for index in indices {
            selected[index] = true;
        }
        Self { selected, _indexes: PhantomData }
    }

    pub fn and(&self, other: &Self) -> Self {
        assert_eq!(self.selected.len(), other.selected.len(), "masks have different length");
        let selected = self.selected.iter().zip(&other.selected).map(|(&a, &b)| a && b).collect();
        Self { selected, _indexes: PhantomData }
    }

    pub fn or(&self, other: &Self) -> Self {
        assert_eq!(self.selected.len(), other.selected.len(), "masks have different length");
        let selected = self.selected.iter().zip(&other.selected).map(|(&a, &b)| a || b).collect();
        Self { selected, _indexes: PhantomData }
    }

    pub fn xor_assign(&mut self, other: &Self) {
        assert_eq!(self.selected.len(), other.selected.len(), "masks have different length");
        for (selected, &other) in self.selected.iter_mut().zip(&other.selected) {
            *selected ^= other;
        }
    }

    pub fn not(&self) -> Self {
        let selected = self.selected.iter().map(|selected| !selected).collect();
        Self { selected, _indexes: PhantomData }
    }

    pub fn indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.selected.iter().enumerate().filter(|(_, &selected)| selected).map(|(index, _)| index)
    }

    pub fn contains(&self, index: usize) -> bool {
        self.selected[index]
    }

    pub fn len(&self) -> usize {
        self.selected.len()
    }
}

// Written by hand because derive would require T: PartialEq, and T is only a tag.
impl<T> PartialEq for Mask<T> {
    fn eq(&self, other: &Self) -> bool {
        self.selected == other.selected
    }
}

impl<T> Clone for Mask<T> {
    fn clone(&self) -> Self {
        Self::from_bools(self.selected.clone())
    }
}

pub struct XorBasis<T> {
    rows: Vec<Mask<T>>,
}

impl<T> XorBasis<T> {
    pub fn new() -> Self {
        Self { rows: vec![] }
    }

    pub fn insert(&mut self, mask: &Mask<T>) -> bool {
        let mut leftover = mask.clone();
        for row in self.rows.iter() {
            let pivot = row.indices().next().unwrap();
            if leftover.contains(pivot) {
                leftover.xor_assign(row);
            }
        }
        let is_new = leftover.indices().next().is_some();
        if is_new {
            self.rows.push(leftover);
        }
        is_new
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }
}
