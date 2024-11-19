use crate::Handle;
use std::{cmp::Ordering, ops::Range};

pub struct Heap<H, Cost>
where
    H: Handle,
    Cost: PartialOrd,
{
    items: Vec<(H, Cost)>,
    index_map: Box<[Option<usize>]>,
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

impl<H, Cost> Heap<H, Cost>
where
    H: Handle,
    Cost: PartialOrd,
{
    pub fn new(num_items: usize) -> Self {
        Heap {
            items: Vec::with_capacity(num_items),
            index_map: vec![None; num_items].into_boxed_slice(),
        }
    }

    fn compare(&self, i: usize, j: usize) -> Option<Ordering> {
        self.items[i].1.partial_cmp(&self.items[j].1)
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.index_map[self.items[i].0.index() as usize] = Some(j);
        self.index_map[self.items[j].0.index() as usize] = Some(i);
        self.items.swap(i, j);
    }

    fn remove_last(&mut self) -> Option<(H, Cost)> {
        let last = self.items.pop();
        if let Some((val, _cost)) = &last {
            self.index_map[val.index() as usize] = None;
        }
        last
    }

    fn sift_up(&mut self, index: usize) {
        let mut index = index;
        while let Some(pi) = heap_parent(index) {
            if let Some(Ordering::Less) = self.compare(index, pi) {
                self.swap(index, pi);
                index = pi;
            } else {
                break;
            }
        }
    }

    fn sift_down(&mut self, index: usize) {
        let mut index = index;
        while index < self.len() {
            match heap_children(index).fold(None, |prev: Option<usize>, ci| {
                if ci < self.len() {
                    match prev {
                        Some(prev) => match self.compare(ci, prev) {
                            Some(Ordering::Less) => Some(ci),
                            _ => Some(prev),
                        },
                        None => Some(ci),
                    }
                } else {
                    prev
                }
            }) {
                Some(child) => match self.compare(index, child) {
                    Some(Ordering::Less) => break,
                    _ => {
                        self.swap(index, child);
                        index = child;
                    }
                },
                None => break,
            }
        }
    }

    pub fn insert(&mut self, val: H, cost: Cost) {
        match self.index_map[val.index() as usize] {
            Some(index) => {
                // Update existing item.
                self.index_map[val.index() as usize] = Some(index);
                self.items[index] = (val, cost);
                self.sift_down(index);
                self.sift_up(index);
            }
            None => {
                // Push new item.
                self.index_map[val.index() as usize] = Some(self.items.len());
                self.items.push((val, cost));
                self.sift_up(self.len() - 1)
            }
        }
    }

    pub fn remove(&mut self, val: H) {
        let index = match self.index_map[val.index() as usize] {
            Some(i) => i,
            None => return, // The item isn't present in the heap.
        };
        let last = self.len() - 1;
        if index == last {
            self.remove_last();
        } else {
            self.swap(index, last);
            self.remove_last();
            self.sift_down(index);
            self.sift_up(index);
        }
    }

    pub fn pop(&mut self) -> Option<(H, Cost)> {
        if self.len() > 1 {
            self.swap(0, self.len() - 1);
            let out = self.remove_last();
            self.sift_down(0);
            out
        } else {
            self.remove_last()
        }
    }

    pub fn clear(&mut self) {
        self.items.clear();
        self.index_map.fill(None);
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }
}

#[cfg(test)]
mod test {
    use super::Heap;
    use crate::{Handle, VH};

    fn drain_heap<H, Cost>(mut heap: Heap<H, Cost>) -> Vec<H>
    where
        H: Handle,
        Cost: PartialOrd,
    {
        let mut out = Vec::with_capacity(heap.len());
        while let Some(val) = heap.pop() {
            out.push(val.0);
        }
        out
    }

    #[test]
    fn t_heap_push_two() {
        let mut heap: Heap<VH, f64> = Heap::new(2);
        heap.insert(0.into(), 5.0);
        heap.insert(1.into(), 2.0);
        assert_eq!(
            vec![1, 0],
            heap.items
                .iter()
                .map(|item| item.0.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_heap_push_many() {
        let mut heap: Heap<VH, f64> = Heap::new(10);
        for i in [8u32, 1, 5, 3, 9, 2, 6, 4, 0, 7] {
            heap.insert(i.into(), i as f64);
        }
        // Push integers in a weird order, and expect them to come out sorted.
        assert_eq!(10, heap.len());
        assert_eq!(
            &(0..10).collect::<Vec<_>>(),
            &drain_heap(heap)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_heap_remove_one() {
        let mut heap: Heap<VH, f64> = Heap::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6] {
            heap.insert(i.into(), i as f64);
        }
        assert_eq!(10, heap.len());
        heap.remove(3.into());
        assert_eq!(
            &vec![0, 1, 2, 4, 5, 6, 7, 8, 9],
            &drain_heap(heap)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_heap_remove_two() {
        let mut heap: Heap<VH, f64> = Heap::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6] {
            heap.insert(i.into(), i as f64);
        }
        assert_eq!(10, heap.len());
        heap.remove(3.into());
        heap.remove(6.into());
        assert_eq!(
            &vec![0, 1, 2, 4, 5, 7, 8, 9],
            &drain_heap(heap)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_heap_update_one() {
        let mut heap: Heap<VH, f64> = Heap::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 6] {
            heap.insert(i.into(), i as f64);
        }
        assert_eq!(9, heap.len());
        heap.insert(0.into(), 0.0);
        assert_eq!(
            &vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &drain_heap(heap)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_heap_update_two() {
        let mut heap: Heap<VH, f64> = Heap::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6] {
            heap.insert(i.into(), i as f64);
        }
        assert_eq!(10, heap.len());
        heap.insert(4.into(), -1.0);
        heap.insert(2.into(), 13.0);
        assert_eq!(
            vec![4, 0, 1, 3, 5, 6, 7, 8, 9, 2],
            drain_heap(heap)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }
}
