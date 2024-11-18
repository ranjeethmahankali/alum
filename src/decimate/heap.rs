use std::cmp::Ordering;

pub struct Heap<T>
where
    T: Copy + Clone + PartialOrd,
{
    items: Vec<T>,
}

fn parent(index: usize) -> Option<usize> {
    if index > 0 {
        Some((index - 1) >> 1)
    } else {
        None
    }
}

fn left(index: usize) -> usize {
    (index << 1) + 1
}

fn right(index: usize) -> usize {
    (index << 1) + 2
}

impl<T> Heap<T>
where
    T: Copy + Clone + PartialOrd,
{
    pub fn new() -> Self {
        Heap { items: Vec::new() }
    }

    fn sift_up(&mut self, index: usize) {
        if index >= self.items.len() {
            return;
        }
        let item = self.items[index];
        let mut index = index;
        while let Some(pi) = parent(index) {
            let mut stop = true;
            if let Some(ord) = self.items[pi].partial_cmp(&item) {
                if let Ordering::Less = ord {
                    stop = false;
                    self.items[index] = self.items[pi];
                    index = pi;
                }
            }
            if stop {
                break;
            }
        }
        self.items[index] = item;
    }

    pub fn push(&mut self, val: T) {
        let idx = self.items.len();
        self.items.push(val);
        self.sift_up(idx);
    }

    pub fn pop(&mut self) -> Option<T> {
        todo!()
    }

    pub fn clear(&mut self) {
        self.items.clear();
    }
}
