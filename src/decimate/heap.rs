use std::{cmp::Ordering, ops::Range};

pub struct Heap<T>
where
    T: PartialOrd,
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

fn children(index: usize) -> Range<usize> {
    let off = index << 1;
    (off + 1)..(off + 3)
}

impl<T> Heap<T>
where
    T: PartialOrd,
{
    pub fn new() -> Self {
        Heap { items: Vec::new() }
    }

    fn sift_up(&mut self, index: usize) -> usize {
        let mut index = index;
        while let Some(pi) = parent(index) {
            if let Some(Ordering::Less) = self.items[index].partial_cmp(&self.items[pi]) {
                self.items.swap(index, pi);
                index = pi;
            } else {
                break;
            }
        }
        index
    }

    fn sift_down(&mut self, index: usize) -> usize {
        let mut index = index;
        while index < self.len() {
            match children(index).fold(None, |prev, ci| {
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
                        self.items.swap(index, child);
                        index = child;
                    }
                },
                None => break,
            }
        }
        index
    }

    pub fn push(&mut self, val: T) -> usize {
        self.items.push(val);
        self.sift_up(self.len() - 1)
    }

    pub fn update(&mut self, index: usize, val: T) -> usize {
        self.items[index] = val;
        let down = self.sift_down(index);
        let up = self.sift_up(index);
        if down == index {
            up
        } else {
            down
        }
    }

    pub fn remove(&mut self, index: usize) {
        let last = self.len() - 1;
        if index == last {
            self.items.pop();
        } else {
            self.items.swap(index, last);
            self.items.pop();
            self.sift_down(index);
            self.sift_up(index);
        }
    }

    pub fn pop(&mut self) -> Option<T> {
        match self.items.pop() {
            Some(last) => {
                if self.items.is_empty() {
                    Some(last)
                } else {
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
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }
}

#[cfg(test)]
mod test {
    use super::Heap;
    use std::fmt::Debug;

    fn drain_heap<T>(mut heap: Heap<T>) -> Vec<T>
    where
        T: Clone + Copy + PartialOrd,
    {
        let mut out = Vec::with_capacity(heap.len());
        while let Some(val) = heap.pop() {
            out.push(val);
        }
        out
    }

    fn heap_from_slice<T>(vals: &[T]) -> Heap<T>
    where
        T: Clone + PartialOrd + Debug,
    {
        let mut heap = Heap::new();
        for val in vals.iter() {
            let index = heap.push(val.clone());
            assert_eq!(&heap.items[index], val);
        }
        heap
    }

    #[test]
    fn t_heap_push_two() {
        let mut heap = Heap::new();
        heap.push(5);
        heap.push(2);
        assert_eq!(vec![2, 5], heap.items);
    }

    #[test]
    fn t_heap_push_many() {
        // Push integers in a weird order, and expect them to come out sorted.
        let heap = heap_from_slice(&[8, 1, 5, 3, 9, 2, 6, 4, 0, 7]);
        assert_eq!(10, heap.len());
        assert_eq!(&(0..10).collect::<Vec<_>>(), &drain_heap(heap));
    }

    #[test]
    fn t_heap_remove_one() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 0, 6]);
        assert_eq!(10, heap.len());
        let index = heap
            .items
            .iter()
            .position(|i| *i == 3)
            .expect("Cannot find the number");
        heap.remove(index);
        assert_eq!(&vec![0, 1, 2, 4, 5, 6, 7, 8, 9], &drain_heap(heap));
    }

    #[test]
    fn t_heap_remove_two() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 0, 6]);
        assert_eq!(10, heap.len());
        heap.remove(
            heap.items
                .iter()
                .position(|i| *i == 3)
                .expect("Cannot find the number"),
        );
        heap.remove(
            heap.items
                .iter()
                .position(|i| *i == 6)
                .expect("Cannot find the number"),
        );
        assert_eq!(&vec![0, 1, 2, 4, 5, 7, 8, 9], &drain_heap(heap));
    }

    #[test]
    fn t_heap_update_one() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 0, 6]);
        assert_eq!(10, heap.len());
        heap.update(
            heap.items
                .iter()
                .position(|i| *i == 4)
                .expect("Cannot find the number"),
            -1,
        );
        assert_eq!(&vec![-1, 0, 1, 2, 3, 5, 6, 7, 8, 9], &drain_heap(heap));
    }

    #[test]
    fn t_heap_update_two() {
        let mut heap = heap_from_slice(&[4, 3, 5, 8, 2, 9, 1, 7, 0, 6]);
        assert_eq!(10, heap.len());
        heap.update(
            heap.items
                .iter()
                .position(|i| *i == 4)
                .expect("Cannot find the number"),
            -1,
        );
        heap.update(
            heap.items
                .iter()
                .position(|i| *i == 2)
                .expect("Cannot find the number"),
            13,
        );
        assert_eq!(&vec![-1, 0, 1, 3, 5, 6, 7, 8, 9, 13], &drain_heap(heap));
    }
}
