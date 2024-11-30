/*!
This priority queue uses an implementation of heap that efficiently allows
updating the priority of an element. It accomplishes this by keeping track of
the positions of all elements in the heap, and then using sift-up and sift-down
operations to update the heap.
*/

use crate::Handle;
use std::{cmp::Ordering, ops::Range};

pub struct Queue<H, Cost>
where
    H: Handle,
    Cost: PartialOrd,
{
    items: Vec<(H, Cost)>,
    map: Vec<Option<usize>>,
}

/// Get the index of the parent in the binary heap.
const fn heap_parent(index: usize) -> Option<usize> {
    if index > 0 {
        Some((index - 1) >> 1)
    } else {
        None
    }
}

/// Get the indices of the children in the binary heap.
const fn heap_children(index: usize) -> Range<usize> {
    let off = index << 1;
    (off + 1)..(off + 3)
}

impl<H, Cost> Queue<H, Cost>
where
    H: Handle,
    Cost: PartialOrd,
{
    /// Create a new heap with the number of elements.
    ///
    /// For example, if the heap is meant to be used to queue up vertices of a
    /// mesh, provide the number of vertices in that mesh.
    pub fn new(num_items: usize) -> Self {
        Queue {
            items: Vec::with_capacity(num_items),
            map: vec![None; num_items],
        }
    }

    /// Clear the queue.
    pub fn clear(&mut self) {
        self.items.clear();
        self.map.fill(None);
    }

    /// Number of items currently in the queue.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Check if this queue is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    // Low level functions for moving things around inside the heap.

    fn compare(&self, i: usize, j: usize) -> Option<Ordering> {
        self.items[i].1.partial_cmp(&self.items[j].1)
    }

    fn position(&self, val: H) -> Option<usize> {
        self.map.get(val.index() as usize).copied().flatten()
    }

    fn set_position(&mut self, val: H, index: usize) {
        let vi = val.index() as usize;
        if vi >= self.map.len() {
            self.map.resize(vi + 1, None);
        }
        unsafe {
            *self.map.get_unchecked_mut(vi) = Some(index);
        }
    }

    fn unset_position(&mut self, val: H) {
        if let Some(i) = self.map.get_mut(val.index() as usize) {
            *i = None;
        }
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.set_position(self.items[i].0, j);
        self.set_position(self.items[j].0, i);
        self.items.swap(i, j);
    }

    fn remove_last(&mut self) -> Option<(H, Cost)> {
        self.items.pop().inspect(|last| {
            self.unset_position(last.0);
        })
    }

    fn write(&mut self, index: usize, val: H, cost: Cost) {
        self.set_position(val, index);
        self.items[index] = (val, cost);
    }

    fn push(&mut self, val: H, cost: Cost) {
        self.set_position(val, self.items.len());
        self.items.push((val, cost));
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

    // High level functions to work with the queue.

    /// Insert the given handle `val` and `cost` into the heap.
    ///
    /// If the handle is already present in the heap, it will be updated to have
    /// the new cost.
    pub fn insert(&mut self, val: H, cost: Cost) {
        match self.position(val) {
            Some(index) => {
                // Update existing item.
                self.write(index, val, cost);
                self.sift_down(index);
                self.sift_up(index);
            }
            None => {
                // Push new item.
                self.push(val, cost);
                self.sift_up(self.len() - 1)
            }
        }
    }

    /// Remove an element from the queue, if it is present.
    pub fn remove(&mut self, val: H) {
        let index = match self.position(val) {
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

    /// Remove the item from the front of the queue and return it.
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
}

#[cfg(test)]
mod test {
    use super::Queue;
    use crate::{Handle, VH};

    fn drain_queue<H, Cost>(mut heap: Queue<H, Cost>) -> Vec<H>
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
    fn t_queue_push_two() {
        let mut queue: Queue<VH, f64> = Queue::new(2);
        queue.insert(0.into(), 5.0);
        queue.insert(1.into(), 2.0);
        assert_eq!(
            vec![1, 0],
            queue
                .items
                .iter()
                .map(|item| item.0.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_queue_push_many() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        for i in [8u32, 1, 5, 3, 9, 2, 6, 4, 0, 7] {
            queue.insert(i.into(), i as f64);
        }
        // Push integers in a weird order, and expect them to come out sorted.
        assert_eq!(10, queue.len());
        assert_eq!(
            &(0..10).collect::<Vec<_>>(),
            &drain_queue(queue)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_queue_remove_one() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6] {
            queue.insert(i.into(), i as f64);
        }
        assert_eq!(10, queue.len());
        queue.remove(3.into());
        assert_eq!(
            &vec![0, 1, 2, 4, 5, 6, 7, 8, 9],
            &drain_queue(queue)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_queue_remove_two() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6] {
            queue.insert(i.into(), i as f64);
        }
        assert_eq!(10, queue.len());
        queue.remove(3.into());
        queue.remove(6.into());
        assert_eq!(
            &vec![0, 1, 2, 4, 5, 7, 8, 9],
            &drain_queue(queue)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_queue_update_one() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 6] {
            queue.insert(i.into(), i as f64);
        }
        assert_eq!(9, queue.len());
        queue.insert(0.into(), 0.0);
        assert_eq!(
            &vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &drain_queue(queue)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_queue_update_two() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        for i in [4, 3, 5, 8, 2, 9, 1, 7, 0, 6] {
            queue.insert(i.into(), i as f64);
        }
        assert_eq!(10, queue.len());
        queue.insert(4.into(), -1.0);
        queue.insert(2.into(), 13.0);
        assert_eq!(
            vec![4, 0, 1, 3, 5, 6, 7, 8, 9, 2],
            drain_queue(queue)
                .iter()
                .map(|v| v.index())
                .collect::<Vec<_>>()
        );
    }
}
