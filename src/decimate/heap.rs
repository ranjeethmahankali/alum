use std::{cmp::Ordering, ops::Range};

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

fn children(index: usize) -> Range<usize> {
    let off = index << 1;
    (off + 1)..(off + 3)
}

impl<T> Heap<T>
where
    T: Copy + Clone + PartialOrd,
{
    pub fn new() -> Self {
        Heap { items: Vec::new() }
    }

    fn sift_up(&mut self, index: usize) -> usize {
        let item = self.items[index];
        let mut index = index;
        while let Some(pi) = parent(index) {
            if let Some(Ordering::Less) = self.items[pi].partial_cmp(&item) {
                self.items[index] = self.items[pi];
                index = pi;
            } else {
                break;
            }
        }
        self.items[index] = item;
        index
    }

    fn sift_down(&mut self, index: usize) -> usize {
        let item = self.items[index];
        let mut index = index;
        while index < self.items.len() {
            match children(index).fold(None, |prev, ci| {
                if ci < self.items.len() {
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
                Some(child) => match item.partial_cmp(&self.items[child]) {
                    Some(Ordering::Less) => break,
                    _ => {
                        self.items[index] = self.items[child];
                        index = child;
                    }
                },
                None => break,
            }
        }
        self.items[index] = item;
        index
    }

    pub fn push(&mut self, val: T) -> usize {
        self.items.push(val);
        self.sift_up(self.items.len() - 1)
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
        let last = self.items.len() - 1;
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
                if self.items.len() > 1 {
                    let out = Some(std::mem::replace(&mut self.items[0], last));
                    self.sift_down(0);
                    out
                } else {
                    self.items.first().copied()
                }
            }
            None => None,
        }
    }

    pub fn clear(&mut self) {
        self.items.clear();
    }
}

#[cfg(test)]
mod test {
    use super::Heap;

    #[test]
    fn t_heap_push() {
        // Push integers in a weird order, and expect them to come out sorted.
        let mut heap = Heap::new();
        for i in [8, 1, 5, 3, 9, 2, 6, 4, 7] {
            eprintln!("Pushing {}", i);
            heap.push(i);
        }
        eprintln!("Done pushing");
        let mut sorted = Vec::new();
        while let Some(i) = heap.pop() {
            eprintln!("Popped {}", i);
            sorted.push(i);
        }
        assert_eq!(&sorted, &(0..10).collect::<Vec<_>>());
    }
}
