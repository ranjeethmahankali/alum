use std::{cmp::Ordering, ops::Range};

pub trait HeapItem: PartialOrd {
    fn index(&self) -> usize;
}

pub struct Heap<T>
where
    T: HeapItem,
{
    items: Vec<T>,
    track: Box<[Option<usize>]>,
}

const fn heap_parent(index: usize) -> Option<usize> {
    if index > 0 {
        Some((index - 1) >> 1)
    } else {
        None
    }
}

const fn heap_children(index: usize) -> Range<usize> {
    let off = index << 1;
    (off + 1)..(off + 3)
}

impl<T> Heap<T>
where
    T: HeapItem,
{
    pub fn new(num_items: usize) -> Self {
        Heap {
            items: Vec::with_capacity(num_items),
            track: vec![None; num_items].into_boxed_slice(),
        }
    }

    fn sift_up(&mut self, index: usize) {
        let mut index = index;
        while let Some(pi) = heap_parent(index) {
            if let Some(Ordering::Less) = self.items[index].partial_cmp(&self.items[pi]) {
                self.track[self.items[index].index()] = Some(pi);
                self.track[self.items[pi].index()] = Some(index);
                self.items.swap(index, pi);
                index = pi;
            } else {
                break;
            }
        }
    }

    fn sift_down(&mut self, index: usize) {
        let mut index = index;
        while index < self.len() {
            match heap_children(index).fold(None, |prev, ci| {
                if ci < self.len() {
                    match prev {
                        Some(prev) => match self.items[ci].partial_cmp(&self.items[prev]) {
                            Some(Ordering::Less) => Some(ci),
                            _ => Some(prev),
                        },
                        None => Some(ci),
                    }
                } else {
                    prev
                }
            }) {
                Some(child) => match self.items[index].partial_cmp(&self.items[child]) {
                    Some(Ordering::Less) => break,
                    _ => {
                        self.track[self.items[index].index()] = Some(child);
                        self.track[self.items[child].index()] = Some(index);
                        self.items.swap(index, child);
                        index = child;
                    }
                },
                None => break,
            }
        }
    }

    pub fn insert(&mut self, val: T) {
        match self.track[val.index()] {
            Some(index) => {
                // Update existing item.
                self.track[val.index()] = Some(index);
                self.items[index] = val;
                self.sift_down(index);
                self.sift_up(index);
            }
            None => {
                // Push new item.
                self.track[val.index()] = Some(self.items.len());
                self.items.push(val);
                self.sift_up(self.len() - 1)
            }
        }
    }

    pub fn remove(&mut self, index: usize) {
        let last = self.len() - 1;
        if index == last {
            if let Some(val) = self.items.pop() {
                self.track[val.index()] = None;
            }
        } else {
            self.track[self.items[index].index()] = Some(last);
            self.track[self.items[last].index()] = Some(index);
            self.items.swap(index, last);
            if let Some(val) = self.items.pop() {
                self.track[val.index()] = None;
            }
            self.sift_down(index);
            self.sift_up(index);
        }
    }

    pub fn pop(&mut self) -> Option<T> {
        match self.items.pop() {
            Some(last) => {
                self.track[last.index()] = None;
                if self.items.is_empty() {
                    Some(last)
                } else {
                    self.track[self.items[0].index()] = None;
                    self.track[last.index()] = Some(0);
                    let out = Some(std::mem::replace(&mut self.items[0], last));
                    self.sift_down(0);
                    out
                }
            }
            None => None,
        }
    }

    pub fn clear(&mut self) {
        self.items.clear();
        self.track.fill(None);
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }
}

#[cfg(test)]
mod test {
    use super::{Heap, HeapItem};
    use std::fmt::Debug;

    fn drain_heap<T>(mut heap: Heap<T>) -> Vec<T>
    where
        T: Clone + Copy + HeapItem,
    {
        let mut out = Vec::with_capacity(heap.len());
        while let Some(val) = heap.pop() {
            out.push(val);
        }
        out
    }

    fn heap_from_slice<T>(vals: &[T], size: usize) -> Heap<T>
    where
        T: Clone + HeapItem + Debug,
    {
        let mut heap = Heap::new(usize::max(size, vals.len()));
        for val in vals.iter() {
            heap.insert(val.clone());
        }
        heap
    }

    impl HeapItem for i32 {
        fn index(&self) -> usize {
            *self as usize
        }
    }

    #[test]
    fn t_heap_push_two() {
        let mut heap = Heap::new(2);
        heap.insert(5);
        heap.insert(2);
        assert_eq!(vec![2, 5], heap.items);
    }

    #[test]
    fn t_heap_push_many() {
        // Push integers in a weird order, and expect them to come out sorted.
        let heap = heap_from_slice(&[8, 1, 5, 3, 9, 2, 6, 4, 0, 7], 10);
        assert_eq!(10, heap.len());
        assert_eq!(&(0..10).collect::<Vec<_>>(), &drain_heap(heap));
    }

    #[test]
    fn t_heap_remove_one() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 0, 6], 10);
        assert_eq!(10, heap.len());
        heap.remove(3);
        assert_eq!(&vec![0, 1, 2, 4, 5, 6, 7, 8, 9], &drain_heap(heap));
    }

    #[test]
    fn t_heap_remove_two() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 0, 6], 10);
        assert_eq!(10, heap.len());
        heap.remove(3);
        heap.remove(6);
        assert_eq!(&vec![0, 1, 2, 4, 5, 7, 8, 9], &drain_heap(heap));
    }

    #[test]
    fn t_heap_update_one() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 6], 10);
        assert_eq!(10, heap.len());
        heap.insert(0);
        assert_eq!(&vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], &drain_heap(heap));
    }

    struct Item {
        idx: usize,
        cost: f64,
    }

    impl PartialEq for Item {
        fn eq(&self, other: &Self) -> bool {
            self.idx == other.idx
        }
    }

    impl PartialOrd for Item {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            self.cost.partial_cmp(&other.cost)
        }
    }

    impl HeapItem for Item {
        fn index(&self) -> usize {
            self.idx
        }
    }

    #[test]
    fn t_heap_update_two() {
        let mut heap = Heap::new(10);
        for item in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6].iter().map(|i| Item {
            idx: *i as usize,
            cost: *i as f64,
        }) {
            heap.insert(item);
        }
        assert_eq!(10, heap.len());
        heap.insert(Item { idx: 4, cost: -1.0 });
        heap.insert(Item { idx: 2, cost: 13.0 });
        assert_eq!(
            &vec![4, 0, 1, 3, 5, 6, 7, 8, 9, 2],
            &heap.items.iter().map(|item| item.idx).collect::<Vec<_>>()
        );
    }
}
