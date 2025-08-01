/*!
This priority queue uses an implementation of heap that efficiently allows
updating the priority of an element. It accomplishes this by keeping track of
the positions of all elements in the heap, and then using sift-up and sift-down
operations to update the heap.
*/

use crate::Handle;
use std::{cmp::Ordering, ops::Range};

/// Priorty queue implementation for mesh elements.
///
/// When you pop an element, the one with the lowest cost is removed from the
/// queue and returned.
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
    use super::{Queue, heap_parent};
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

    // Empty queue operations
    #[test]
    fn t_empty_queue_pop() {
        let mut queue: Queue<VH, f64> = Queue::new(0);
        assert_eq!(None, queue.pop());
        assert!(queue.is_empty());
    }

    #[test]
    fn t_empty_queue_remove() {
        let mut queue: Queue<VH, f64> = Queue::new(0);
        queue.remove(0.into()); // Should not panic
        assert!(queue.is_empty());
    }

    #[test]
    fn t_empty_queue_insert_remove() {
        let mut queue: Queue<VH, f64> = Queue::new(0);
        queue.insert(5.into(), 10.0);
        assert_eq!(1, queue.len());
        queue.remove(5.into());
        assert!(queue.is_empty());
    }

    // Single element operations
    #[test]
    fn t_single_element_update_cost() {
        let mut queue: Queue<VH, f64> = Queue::new(1);
        queue.insert(0.into(), 5.0);
        queue.insert(0.into(), 10.0); // Update cost
        assert_eq!(Some((0.into(), 10.0)), queue.pop());
        assert!(queue.is_empty());
    }

    #[test]
    fn t_single_element_remove_reinsert() {
        let mut queue: Queue<VH, f64> = Queue::new(1);
        queue.insert(0.into(), 5.0);
        queue.remove(0.into());
        assert!(queue.is_empty());
        queue.insert(0.into(), 10.0);
        assert_eq!(Some((0.into(), 10.0)), queue.pop());
    }

    // Handle index edge cases
    #[test]
    fn t_large_handle_indices() {
        // This ensures the actual indices of the handle are not used in memory related arithmetic.
        let mut queue: Queue<VH, f64> = Queue::new(2);
        queue.insert(1000.into(), 1.0);
        queue.insert(2000.into(), 2.0);
        let result = drain_queue(queue);
        assert_eq!(
            vec![1000, 2000],
            result.iter().map(|v| v.index()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_handle_index_zero_vs_others() {
        let mut queue: Queue<VH, f64> = Queue::new(3);
        queue.insert(0.into(), 3.0);
        queue.insert(1.into(), 1.0);
        queue.insert(2.into(), 2.0);
        let result = drain_queue(queue);
        assert_eq!(
            vec![1, 2, 0],
            result.iter().map(|v| v.index()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_sparse_handle_indices() {
        let mut queue: Queue<VH, f64> = Queue::new(3);
        queue.insert(100.into(), 3.0);
        queue.insert(50.into(), 1.0);
        queue.insert(200.into(), 2.0);
        let result = drain_queue(queue);
        assert_eq!(
            vec![50, 200, 100],
            result.iter().map(|v| v.index()).collect::<Vec<_>>()
        );
    }

    // Cost comparison edge cases
    #[test]
    fn t_identical_costs() {
        let mut queue: Queue<VH, f64> = Queue::new(3);
        queue.insert(0.into(), 5.0);
        queue.insert(1.into(), 5.0);
        queue.insert(2.into(), 5.0);
        let result = drain_queue(queue);
        assert_eq!(3, result.len());
        // All should have same cost, order may vary but all should be present
        let mut indices: Vec<_> = result.iter().map(|v| v.index()).collect();
        indices.sort();
        assert_eq!(vec![0, 1, 2], indices);
    }

    #[test]
    fn t_special_float_values() {
        let mut queue: Queue<VH, f64> = Queue::new(4);
        queue.insert(0.into(), f64::INFINITY);
        queue.insert(1.into(), f64::NEG_INFINITY);
        queue.insert(2.into(), 0.0);
        queue.insert(3.into(), -0.0);

        // NEG_INFINITY should come first, then zeros, then INFINITY
        let result: Vec<_> = (0..4).map(|_| queue.pop().unwrap().0.index()).collect();
        assert_eq!(1, result[0]); // NEG_INFINITY first
        assert_eq!(0, result[3]); // INFINITY last
    }

    #[test]
    fn t_nan_costs() {
        let mut queue: Queue<VH, f64> = Queue::new(2);
        queue.insert(0.into(), f64::NAN);
        queue.insert(1.into(), 5.0);

        // With NaN, behavior is implementation-defined but shouldn't panic
        let result1 = queue.pop();
        let result2 = queue.pop();
        assert!(result1.is_some());
        assert!(result2.is_some());
    }

    #[test]
    fn t_negative_and_zero_costs() {
        let mut queue: Queue<VH, f64> = Queue::new(4);
        queue.insert(0.into(), -5.0);
        queue.insert(1.into(), 0.0);
        queue.insert(2.into(), 5.0);
        queue.insert(3.into(), -10.0);

        let result = drain_queue(queue);
        assert_eq!(
            vec![3, 0, 1, 2],
            result.iter().map(|v| v.index()).collect::<Vec<_>>()
        );
    }

    // Complex operation sequences
    #[test]
    fn t_insert_remove_reinsert_different_cost() {
        let mut queue: Queue<VH, f64> = Queue::new(5);
        queue.insert(0.into(), 10.0);
        queue.insert(1.into(), 5.0);
        queue.remove(0.into());
        queue.insert(0.into(), 1.0); // Lower cost

        let result = drain_queue(queue);
        assert_eq!(
            vec![0, 1],
            result.iter().map(|v| v.index()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_multiple_updates_same_element() {
        let mut queue: Queue<VH, f64> = Queue::new(3);
        queue.insert(0.into(), 10.0);
        queue.insert(1.into(), 5.0);
        queue.insert(0.into(), 15.0); // Higher cost
        queue.insert(0.into(), 2.0); // Lower cost
        queue.insert(0.into(), 7.0); // Medium cost

        let result = drain_queue(queue);
        assert_eq!(
            vec![1, 0],
            result.iter().map(|v| v.index()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn t_remove_nonexistent_twice() {
        let mut queue: Queue<VH, f64> = Queue::new(5);
        queue.insert(0.into(), 5.0);
        queue.remove(1.into()); // Remove non-existent
        queue.remove(1.into()); // Remove non-existent again
        assert_eq!(1, queue.len());
        assert_eq!(Some((0.into(), 5.0)), queue.pop());
    }

    #[test]
    fn t_update_vs_insert_nonexistent() {
        let mut queue: Queue<VH, f64> = Queue::new(5);
        queue.insert(0.into(), 5.0);
        queue.insert(1.into(), 10.0); // This should insert, not update
        assert_eq!(2, queue.len());
    }

    #[test]
    fn t_alternating_insert_remove() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        // Insert 0,1,2
        for i in 0..3 {
            queue.insert(i.into(), i as f64);
        }
        // Remove 1, insert 3
        queue.remove(1.into());
        queue.insert(3.into(), 3.0);
        // Remove 0, insert 4
        queue.remove(0.into());
        queue.insert(4.into(), 1.5);
        let result = drain_queue(queue);
        let mut indices: Vec<_> = result.iter().map(|v| v.index()).collect();
        indices.sort();
        assert_eq!(vec![2, 3, 4], indices);
    }

    // Heap property verification helper
    fn verify_heap_property<H: Handle, Cost: PartialOrd>(queue: &Queue<H, Cost>) -> bool {
        for i in 0..queue.len() {
            if let Some(parent_idx) = heap_parent(i) {
                if let Some(std::cmp::Ordering::Greater) = queue.compare(parent_idx, i) {
                    return false; // Parent is greater than child - heap property violated
                }
            }
        }
        true
    }

    fn verify_position_map<H: Handle, Cost: PartialOrd>(queue: &Queue<H, Cost>) -> bool {
        for (actual_pos, (handle, _)) in queue.items.iter().enumerate() {
            if let Some(mapped_pos) = queue.position(*handle) {
                if mapped_pos != actual_pos {
                    return false; // Position map inconsistency
                }
            } else {
                return false; // Handle not found in position map
            }
        }
        true
    }

    #[test]
    fn t_heap_property_after_complex_operations() {
        let mut queue: Queue<VH, f64> = Queue::new(20);
        // Insert many elements
        for i in [
            15, 3, 8, 12, 1, 18, 7, 4, 11, 9, 2, 14, 6, 10, 5, 13, 16, 17, 19,
        ] {
            queue.insert(i.into(), i as f64);
            assert!(
                verify_heap_property(&queue),
                "Heap property violated after inserting {}",
                i
            );
            assert!(
                verify_position_map(&queue),
                "Position map inconsistent after inserting {}",
                i
            );
        }
        // Remove some elements
        for i in [3, 12, 7, 11, 14, 16] {
            queue.remove(i.into());
            assert!(
                verify_heap_property(&queue),
                "Heap property violated after removing {}",
                i
            );
            assert!(
                verify_position_map(&queue),
                "Position map inconsistent after removing {}",
                i
            );
        }
        // Update some elements
        queue.insert(15.into(), 0.5); // Make it smallest
        assert!(verify_heap_property(&queue));
        assert!(verify_position_map(&queue));
        queue.insert(1.into(), 25.0); // Make it largest
        assert!(verify_heap_property(&queue));
        assert!(verify_position_map(&queue));
    }

    #[test]
    fn t_position_map_consistency() {
        let mut queue: Queue<VH, f64> = Queue::new(10);
        // Insert elements
        for i in [5, 2, 8, 1, 9, 3, 7] {
            queue.insert(i.into(), i as f64);
        }
        // Verify every element in items has correct position in map
        assert!(verify_position_map(&queue));
        // Remove some and verify again
        queue.remove(5.into());
        queue.remove(8.into());
        assert!(verify_position_map(&queue));
        // Update some and verify again
        queue.insert(2.into(), 10.0);
        assert!(verify_position_map(&queue));
    }

    // Stress testing
    #[test]
    fn t_stress_random_operations() {
        let mut queue: Queue<VH, f64> = Queue::new(100);
        let mut expected_elements = std::collections::HashSet::new();
        // Insert 50 elements
        for i in 0..50 {
            queue.insert(i.into(), (i * 7 % 23) as f64); // Some pseudo-random costs
            expected_elements.insert(i);
        }
        // Remove every 3rd element
        for i in (0..50).step_by(3) {
            queue.remove(i.into());
            expected_elements.remove(&i);
        }
        // Update every 5th remaining element
        for i in (0..50).step_by(5) {
            if expected_elements.contains(&i) {
                queue.insert(i.into(), (i as f64) * 0.5);
            }
        }
        // Verify length matches expected
        assert_eq!(expected_elements.len(), queue.len());
        // Verify all expected elements are present and in sorted order
        let mut last_cost = f64::NEG_INFINITY;
        while let Some((handle, cost)) = queue.pop() {
            assert!(expected_elements.remove(&handle.index()));
            assert!(cost >= last_cost, "Queue not properly ordered");
            last_cost = cost;
        }
        assert!(
            expected_elements.is_empty(),
            "Some elements missing from queue"
        );
    }

    #[test]
    fn t_stress_many_identical_costs() {
        let mut queue: Queue<VH, f64> = Queue::new(100);
        // Insert many elements with same cost
        for i in 0..50 {
            queue.insert(i.into(), 5.0);
        }
        // Remove some
        for i in (10..30).step_by(2) {
            queue.remove(i.into());
        }
        // All remaining should have cost 5.0
        let mut count = 0;
        while let Some((_, cost)) = queue.pop() {
            assert_eq!(5.0, cost);
            count += 1;
        }
        // Should have 50 - 10 = 40 elements remaining
        assert_eq!(40, count);
    }
}
