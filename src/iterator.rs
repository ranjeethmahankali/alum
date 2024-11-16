use crate::{
    element::{EH, FH, HH, VH},
    topol::Topology,
    HasTopology,
};
use std::{marker::PhantomData, ptr::NonNull};

/*
Halfedges are at the core of all iterators, because halfedges have the most
amount of connectivity information than any other type of mesh element. There
are two base types of iterators: Radial and Loop. Radial iterators go over the
halfedges rotated around their tail vertex. Loop iterators traverse the
next/previous halfedge links, hence traverse loops. Both of these iterator types
can go either counter-clockwise or clockwise. Further more, they can either
borrow the mesh immutably or mutably. All other iterators such as face-vertex,
and vertex-vertex etc. can be expressed in terms of these base iterators using
`.map()`.
 */

/// Iterator over halfedges radiating outward from a vertex.
struct RadialHalfedgeIter<'a, const CCW: bool> {
    topol: &'a Topology,
    hstart: Option<HH>,
    hcurrent: Option<HH>,
}

impl<'a, const CCW: bool> RadialHalfedgeIter<'a, CCW> {
    fn new(topol: &'a Topology, h: Option<HH>) -> Self {
        RadialHalfedgeIter {
            topol,
            hstart: h,
            hcurrent: h,
        }
    }
}

impl<'a> Iterator for RadialHalfedgeIter<'a, true> {
    type Item = HH;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.prev(self.topol).opposite();
                self.hcurrent = match self.hstart {
                    Some(start) if start != next => Some(next),
                    _ => None,
                };
                Some(current)
            }
            None => None,
        }
    }
}

impl<'a> Iterator for RadialHalfedgeIter<'a, false> {
    type Item = HH;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.opposite().next(self.topol);
                self.hcurrent = match self.hstart {
                    Some(start) if start != next => Some(next),
                    _ => None,
                };
                Some(current)
            }
            None => None,
        }
    }
}

/// Iterator for halfedges in a loop.
struct LoopHalfedgeIter<'a, const CCW: bool> {
    topol: &'a Topology,
    hstart: HH,
    hcurrent: Option<HH>,
}

impl<'a, const CCW: bool> LoopHalfedgeIter<'a, CCW> {
    fn new(topol: &'a Topology, h: HH) -> Self {
        LoopHalfedgeIter {
            topol,
            hstart: h,
            hcurrent: Some(h),
        }
    }
}

impl<'a> Iterator for LoopHalfedgeIter<'a, true> {
    type Item = HH;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.next(self.topol);
                self.hcurrent = if next == self.hstart {
                    None
                } else {
                    Some(next)
                };
                Some(current)
            }
            None => None,
        }
    }
}

impl<'a> Iterator for LoopHalfedgeIter<'a, false> {
    type Item = HH;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.prev(self.topol);
                self.hcurrent = if next == self.hstart {
                    None
                } else {
                    Some(next)
                };
                Some(current)
            }
            None => None,
        }
    }
}

struct RadialHalfedgeIterMut<'a, const CCW: bool, T> {
    reference: NonNull<T>,
    topol: &'a mut Topology,
    hstart: Option<HH>,
    hcurrent: Option<HH>,
    _phantom: PhantomData<&'a mut T>,
}

impl<'a, T> Iterator for RadialHalfedgeIterMut<'a, true, T> {
    type Item = (&'a mut T, HH);

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.prev(self.topol).opposite();
                self.hcurrent = match self.hstart {
                    Some(start) if start != next => Some(next),
                    _ => None,
                };
                Some((unsafe { self.reference.as_mut() }, current))
            }
            None => None,
        }
    }
}

impl<'a, T> Iterator for RadialHalfedgeIterMut<'a, false, T> {
    type Item = (&'a mut T, HH);

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.opposite().next(self.topol);
                self.hcurrent = match self.hstart {
                    Some(start) if start != next => Some(next),
                    _ => None,
                };
                Some((unsafe { self.reference.as_mut() }, current))
            }
            None => None,
        }
    }
}

struct LoopHalfedgeIterMut<'a, const CCW: bool, T> {
    reference: NonNull<T>,
    topol: &'a mut Topology,
    hstart: HH,
    hcurrent: Option<HH>,
    _phantom: PhantomData<&'a mut T>,
}

impl<'a, T> Iterator for LoopHalfedgeIterMut<'a, true, T> {
    type Item = (&'a mut T, HH);

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.next(self.topol);
                self.hcurrent = if next == self.hstart {
                    None
                } else {
                    Some(next)
                };
                Some((unsafe { self.reference.as_mut() }, current))
            }
            None => None,
        }
    }
}

impl<'a, T> Iterator for LoopHalfedgeIterMut<'a, false, T> {
    type Item = (&'a mut T, HH);

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = current.prev(self.topol);
                self.hcurrent = if next == self.hstart {
                    None
                } else {
                    Some(next)
                };
                Some((unsafe { self.reference.as_mut() }, current))
            }
            None => None,
        }
    }
}

/// This trait defines iterators to traverse the topology of a mesh.
pub trait HasIterators: HasTopology {
    /*!
    # Working with mutable iterators

    Iterators such as [`Self::loop_ccw_iter`] and other `*_iter` iterators
    borrow the mesh immutably. These are useful when iterating over the elements
    of the mesh without needing to modify the mesh. While the mesh is borrowed
    by the iterators, the borrow checker will not let you borrow the mesh again
    mutably, for as long as the iterator is alive. This is a problem if you're
    trying to modify the mesh while iterating over it's elements. In such
    scenarios, mutable iterators are useful. These functions end with
    `*_iter_mut`.

    The `*_iter_mut` functions borrow the mesh mutably. Similar to immutable
    iterators, while the mutable iterator is alive, the borrow checker won't let
    you mutably borrow the mesh again. But, unlike the immutable iterators, the
    mutable iterators don't yield just the elements of the mesh. Instead the
    mutable iterators yield a tuple containing a mutable reference to the mesh,
    and the mesh element. You can modify the mesh using the mutable reference
    yielded by the mutable iterator. So essentially, the borrow checker is happy
    because the iterator borrows the mesh mutably and to modify the mesh, you
    borrow the mesh from the iterator. The borrows are chained instead of being
    simultaenous. This keeps the borrow checker happy, and ensures safety to
    some extent. Below is some example code that iterates over the vertices of
    face with index 2, and modifies their positions.

    ```rust
    use alum::{alum_glam::PolyMeshF32, FH, HasIterators, HasTopology};

    let mut boxmesh = PolyMeshF32::unit_box().expect("Cannot create a box");
    assert_eq!(1.0, boxmesh.try_calc_volume().expect("Cannot compute volume"));
    // Modify the mesh while using a mutable iterator - pull points closer to origin.
    let f: FH = 2.into();
    for (mesh, v) in boxmesh.fv_ccw_iter_mut(f) {
        // Inside this loop, while the iterator is alive, we cannot borrow `boxmesh`
        // because the iterator already borrowed `boxmesh` mutably. Instead we will use
        // the `mesh` mutable reference yielded by the iterator along with the halfedge.
        let mut pos = mesh.point(v).expect("Cannot read position");
        pos *= 0.75;
        mesh.set_point(v, pos);
    }
    // The iterator is not alive anymore, so we can borrow `boxmesh` again.
    boxmesh.check_topology().expect("Topological errors found");
    // Volume is smaller than one because we pulled the vertices closer to the origin.
    assert!(1.0 > boxmesh.try_calc_volume().expect("Cannot compute volume"));
    ```

    Despite enforcing the borrow checker rules, the mutable iterators can lead
    to problems when used incorrectly. Modifying the topology of the mesh while
    using mutable iterators is NOT advised, as this can lead to topological
    errors. Only do this if you're know what you're doing. This is akin to
    mutably iterating over a linked list while modifying the links between the
    elements of the linked list. You can easily create cycles, infinite loops
    and other problems if you're not careful.
     */

    /// Iterator over the outgoing halfedges around a vertex, going counter-clockwise.
    fn voh_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> {
        RadialHalfedgeIter::<true>::new(self.topology(), v.halfedge(self.topology()))
    }

    /// Iterator over the outgoing halfedges around a vertex, going clockwise
    fn voh_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> {
        RadialHalfedgeIter::<false>::new(self.topology(), v.halfedge(self.topology()))
    }

    /// Iterator over the incoming halfedges around a vertex, going
    /// counter-clockwise
    fn vih_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> {
        self.voh_ccw_iter(v).map(|h| h.opposite())
    }

    /// Iterator over the incoming halfedges around a vertex, going clockwise
    fn vih_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> {
        self.voh_cw_iter(v).map(|h| h.opposite())
    }

    /// Iterator over the faces incident on a vertex, going counter-clockwise.
    fn vf_ccw_iter(&self, v: VH) -> impl Iterator<Item = FH> {
        self.voh_ccw_iter(v).filter_map(|h| h.face(self.topology()))
    }

    /// Iterator over the faces incident on a vertex, going clockwise.
    fn vf_cw_iter(&self, v: VH) -> impl Iterator<Item = FH> {
        self.voh_cw_iter(v).filter_map(|h| h.face(self.topology()))
    }

    /// Iterator over the neighboring vertices around the given vertex, going
    /// counter-clockwise.
    fn vv_ccw_iter(&self, v: VH) -> impl Iterator<Item = VH> {
        self.voh_ccw_iter(v).map(|h| h.head(self.topology()))
    }

    /// Iterator over the neighboring vertices around the given vertex, going
    /// clockwise.
    fn vv_cw_iter(&self, v: VH) -> impl Iterator<Item = VH> {
        self.voh_cw_iter(v).map(|h| h.head(self.topology()))
    }

    /// Iterator over the incident edges around an vertex, going counter-clockwise.
    fn ve_ccw_iter(&self, v: VH) -> impl Iterator<Item = EH> {
        self.voh_ccw_iter(v).map(|h| h.edge())
    }

    /// Iterator over the incident edges around a vertex, going clockwise.
    fn ve_cw_iter(&self, v: VH) -> impl Iterator<Item = EH> {
        self.voh_cw_iter(v).map(|h| h.edge())
    }

    /// Iterator over the two vertices incident on the given edge.
    fn ev_iter(&self, e: EH) -> impl Iterator<Item = VH> {
        self.eh_iter(e).map(|h| h.head(self.topology()))
    }

    /// Iterator over the two halfedges corresponding to an edge.
    fn eh_iter(&self, e: EH) -> impl Iterator<Item = HH> {
        [false, true].iter().map(move |flag| e.halfedge(*flag))
    }

    /// Iterator over the faces incident on an edge.
    fn ef_iter(&self, e: EH) -> impl Iterator<Item = FH> {
        self.eh_iter(e).filter_map(|h| h.face(self.topology()))
    }

    /// Iterator over the halfedges of a face loop, going counter-clockwise.
    fn fh_ccw_iter(&self, f: FH) -> impl Iterator<Item = HH> {
        LoopHalfedgeIter::<true>::new(self.topology(), f.halfedge(self.topology()))
    }

    /// Iterator over the halfedges of a face loop, going clockwise.
    fn fh_cw_iter(&self, f: FH) -> impl Iterator<Item = HH> {
        LoopHalfedgeIter::<false>::new(self.topology(), f.halfedge(self.topology()))
    }

    /// Iterator over the vertices incident on a face, going counter-clockwise.
    fn fv_ccw_iter(&self, f: FH) -> impl Iterator<Item = VH> {
        self.fh_ccw_iter(f).map(|h| h.head(self.topology()))
    }

    /// Iterator over the vertices incident on a face, going clockwise.
    fn fv_cw_iter(&self, f: FH) -> impl Iterator<Item = VH> {
        self.fh_cw_iter(f).map(|h| h.head(self.topology()))
    }

    /// Iterator over the edges incident on a face, going counter-clockwise.
    fn fe_ccw_iter(&self, f: FH) -> impl Iterator<Item = EH> {
        self.fh_ccw_iter(f).map(|h| h.edge())
    }

    /// Iterator over the edges incident on a face, going clockwise.
    fn fe_cw_iter(&self, f: FH) -> impl Iterator<Item = EH> {
        self.fh_cw_iter(f).map(|h| h.edge())
    }

    /// Iterator over the neighboring faces arouund the given face, going
    /// counter-clockwise.
    ///
    /// This includes the faces connected via a shared edge, but not those
    /// connected via a shared vertex.
    fn ff_ccw_iter(&self, f: FH) -> impl Iterator<Item = FH> {
        self.fh_ccw_iter(f)
            .filter_map(|h| h.opposite().face(self.topology()))
    }

    /// Iterator over the neighboring faces around the given face, going
    /// clockwise.
    ///
    /// This includes the faces connected via a shared edge, but not those
    /// connected via a shared vertex.
    fn ff_cw_iter(&self, f: FH) -> impl Iterator<Item = FH> {
        self.fh_cw_iter(f)
            .filter_map(|h| h.opposite().face(self.topology()))
    }

    /// This is similar to [`Self::voh_ccw_iter`] around the tail of the given
    /// halfedge, except this iterator starts at the provided halfedge.
    ///
    /// This is equivalent to a circular shifted [`Self::voh_ccw_iter`] of the
    /// vertex at the tail of this halfedge.
    fn ccw_rotate_iter(&self, h: HH) -> impl Iterator<Item = HH> {
        RadialHalfedgeIter::<true>::new(self.topology(), Some(h))
    }

    /// This is similar to [`Self::voh_cw_iter`] around the tail of the given
    /// halfedge, except this iterator starts at the provided halfedge.
    ///
    /// This is equivalent to a circular shifted [`Self::voh_cw_iter`] of the
    /// vertex at the tail of this halfedge.
    fn cw_rotate_iter(&self, h: HH) -> impl Iterator<Item = HH> {
        RadialHalfedgeIter::<false>::new(self.topology(), Some(h))
    }

    /// Counter-clockwise iterator over the halfedges in a loop.
    ///
    /// The iterator will start at the given halfedge. If the halfedge has an
    /// incident face, this iterator is equivalent to a circular shifted
    /// [`Self::fh_ccw_iter`] of the incident face. If the halfedge is on the
    /// boundary, this iterator goes over the boundary loop counter-clockwise.
    fn loop_ccw_iter(&self, h: HH) -> impl Iterator<Item = HH> {
        LoopHalfedgeIter::<true>::new(self.topology(), h)
    }

    /// Counter-clockwise iterator over the halfedges in a loop.
    ///
    /// The iterator will start at the given halfedge. If the halfedge has an
    /// incident face, this iterator is equivalent to a circular shifted
    /// [`Self::fh_cw_iter`] of the incident face. If the halfedge is on the
    /// boundary, this iterator goes over the boundary loop clockwise.
    fn loop_cw_iter(&self, h: HH) -> impl Iterator<Item = HH> {
        LoopHalfedgeIter::<false>::new(self.topology(), h)
    }

    /// Iterator over the outgoing halfedges around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn voh_ccw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, HH)> {
        let h = v.halfedge(self.topology());
        RadialHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: h,
            hcurrent: h,
            _phantom: PhantomData,
        }
    }

    /// Iterator over the outgoing halfedges around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn voh_cw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, HH)> {
        let h = v.halfedge(self.topology());
        RadialHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: h,
            hcurrent: h,
            _phantom: PhantomData,
        }
    }

    /// Iterator over the incoming halfedges around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn vih_ccw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, HH)> {
        self.voh_ccw_iter_mut(v)
            .map(|(mesh, h)| (mesh, h.opposite()))
    }

    /// Iterator over the incoming halfedges around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn vih_cw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, HH)> {
        self.voh_cw_iter_mut(v)
            .map(|(mesh, h)| (mesh, h.opposite()))
    }

    /// Mutable iterator over the neighboring vertices around the given vertex,
    /// going counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn vv_ccw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, VH)> {
        self.voh_ccw_iter_mut(v).map(|(mesh, h)| {
            let v = h.head(mesh);
            (mesh, v)
        })
    }

    /// Mutable iterator over the neighboring vertices around the given vertex,
    /// going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn vv_cw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, VH)> {
        self.voh_cw_iter_mut(v).map(|(mesh, h)| {
            let v = h.head(mesh);
            (mesh, v)
        })
    }

    /// Iterator over the incident edges around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn ve_ccw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, EH)> {
        self.voh_ccw_iter_mut(v).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the incident edges around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn ve_cw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, EH)> {
        self.voh_cw_iter_mut(v).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the incident faces around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn vf_ccw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, FH)> {
        self.voh_ccw_iter_mut(v)
            .filter_map(|(mesh, h)| h.face(mesh).map(|f| (mesh, f)))
    }

    /// Iterator over the incident faces around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn vf_cw_iter_mut(&mut self, v: VH) -> impl Iterator<Item = (&mut Self, FH)> {
        self.voh_cw_iter_mut(v)
            .filter_map(|(mesh, h)| h.face(mesh).map(|f| (mesh, f)))
    }

    /// Iterator over the vertices incident on a face, going counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn fv_ccw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, VH)> {
        self.fh_ccw_iter_mut(f).map(|(mesh, h)| {
            let v = h.head(mesh);
            (mesh, v)
        })
    }

    /// Iterator over the vertices incident on a face, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn fv_cw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, VH)> {
        self.fh_cw_iter_mut(f).map(|(mesh, h)| {
            let v = h.head(mesh);
            (mesh, v)
        })
    }

    /// Iterator over the halfedges of a face loop, going counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn fh_ccw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, HH)> {
        let h = f.halfedge(self);
        LoopHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: h,
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }

    /// Iterator over the halfedges in a face loop, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn fh_cw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, HH)> {
        let h = f.halfedge(self);
        LoopHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: h,
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }

    /// Iterator over the edges incident on a face, going counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn fe_ccw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, EH)> {
        self.fh_ccw_iter_mut(f).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the edges incident on a face, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn fe_cw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, EH)> {
        self.fh_cw_iter_mut(f).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the neighboring faces around a given face, going
    /// counter-clockwise.
    ///
    /// This includes the faces connected via a shared edge, but not those
    /// connected via shared edge.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn ff_ccw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, FH)> {
        self.fh_ccw_iter_mut(f)
            .filter_map(|(mesh, h)| h.opposite().face(mesh).map(|f| (mesh, f)))
    }

    /// Iterator over the neighboring faces around a given face, going
    /// clockwise.
    ///
    /// This includes the faces connected via a shared edge, but not those
    /// connected via shared edge.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    fn ff_cw_iter_mut(&mut self, f: FH) -> impl Iterator<Item = (&mut Self, FH)> {
        self.fh_cw_iter_mut(f)
            .filter_map(|(mesh, h)| h.opposite().face(mesh).map(|f| (mesh, f)))
    }

    /// This is similar to [`Self::voh_ccw_iter_mut`] around the tail of the
    /// given halfedge, except this iterator starts at the provided halfedge.
    ///
    /// This is equivalent to a circular shifted [`Self::voh_ccw_iter_mut`] of
    /// the vertex at the tail of the given halfedge. A mutable iterator allows
    /// modifying the mesh while iterating over its elements. See [this
    /// section](#working-with-mutable-iterators) for an overview on
    /// mutable iterators.
    fn ccw_rotate_iter_mut(&mut self, h: HH) -> impl Iterator<Item = (&mut Self, HH)> {
        RadialHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: Some(h),
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }

    /// This is similar to [`Self::voh_cw_iter_mut`] around the tail of the
    /// given halfedge, except this iterator starts at the provided halfedge.
    ///
    /// This is equivalent to a circular shifted [`Self::voh_cw_iter_mut`] of
    /// the vertex at the tail of the given halfedge. A mutable iterator allows
    /// modifying the mesh while iterating over its elements. See [this
    /// section](#working-with-mutable-iterators) for an overview on
    /// mutable iterators.
    fn cw_rotate_iter_mut(&mut self, h: HH) -> impl Iterator<Item = (&mut Self, HH)> {
        RadialHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: Some(h),
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }

    /// Counter-clockwise iterator over the halfedges in a loop.
    ///
    /// The iterator will start at the given halfedge. If the halfedge has an
    /// incident face, this iterator is equivalent to a circular shifted
    /// [`Self::fh_ccw_iter_mut`] of the incident face. If the halfedge is on
    /// the boundary, this iterator goes over the boundary loop
    /// counter-clockwise. A mutable iterator allows modifying the mesh while
    /// iterating over its elements. See [this
    /// section](#working-with-mutable-iterators) for an overview on
    /// mutable iterators.
    fn loop_ccw_iter_mut(&mut self, h: HH) -> impl Iterator<Item = (&mut Self, HH)> {
        LoopHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: h,
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }

    /// Counter-clockwise iterator over the halfedges in a loop.
    ///
    /// The iterator will start at the given halfedge. If the halfedge has an
    /// incident face, this iterator is equivalent to a circular shifted
    /// [`Self::fh_cw_iter_mut`] of the incident face. If the halfedge is on the
    /// boundary, this iterator goes over the boundary loop clockwise. A mutable
    /// iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for
    /// an overview on mutable iterators.
    fn loop_cw_iter_mut(&mut self, h: HH) -> impl Iterator<Item = (&mut Self, HH)> {
        LoopHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: self.topology_mut(),
            hstart: h,
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }
    /// Iterator over the vertex triplets that represent a triangulation of the
    /// given face.
    ///
    /// The triangulation does not take the shape of the face into account. It
    /// only accounts for the topology of the face.
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, Handle, HasTopology, HasIterators};
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_quad_face(0.into(), 1.into(), 2.into(), 3.into());
    /// assert_eq!(mesh.triangulated_face_vertices(0.into())
    ///                .flatten()
    ///                .map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [3, 0, 1, 3, 1, 2]);
    /// ```
    fn triangulated_face_vertices(&self, f: FH) -> impl Iterator<Item = [VH; 3]> {
        let hstart = f.halfedge(self);
        let vstart = hstart.tail(self);
        self.loop_ccw_iter(hstart.next(self))
            .take_while(move |h| h.head(self) != vstart)
            .map(move |h| [vstart, h.tail(self), h.head(self)])
    }

    /// Iterator over the vertex triplets that represent a triangulation of this
    /// mesh. The triangulation of a face does not take it's shape into
    /// account. It only accounts for the topology.
    ///
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, Handle, HasTopology, HasIterators};
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_quad_face(0.into(), 1.into(), 2.into(), 3.into());
    /// assert_eq!(mesh.triangulated_vertices()
    ///                .flatten()
    ///                .map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [3, 0, 1, 3, 1, 2]);
    /// ```
    fn triangulated_vertices(&self) -> impl Iterator<Item = [VH; 3]> {
        self.faces()
            .flat_map(move |f| self.triangulated_face_vertices(f))
    }

    /// Find a halfedge spanning the vertices `from` and `to`, if one exists
    fn find_halfedge(&self, from: VH, to: VH) -> Option<HH> {
        self.voh_ccw_iter(from)
            .find(|h| h.head(self.topology()) == to)
    }

    /// Adjust the outgoing halfedge of a vertex, if necessary.
    ///
    /// If a vertex is on the boundary, it's outgoing halfedge MUST be a
    /// boundary halfedge.
    fn adjust_outgoing_halfedge(&mut self, v: VH) {
        let h = self.voh_ccw_iter(v).find(|h| h.is_boundary(self));
        if let Some(h) = h {
            self.topology_mut().vertex_mut(v).halfedge = Some(h);
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{
        alum_glam::PolyMeshF32,
        element::{Handle, VH},
        iterator::HasIterators,
        topol::{HasTopology, TopolCache, Topology},
    };
    use arrayvec::ArrayVec;

    /**
     * Makes a box with the following topology.
     * ```text
     *
     *      7-----------6
     *     /|          /|
     *    / |         / |
     *   4-----------5  |
     *   |  |        |  |
     *   |  3--------|--2
     *   | /         | /
     *   |/          |/
     *   0-----------1
     * ```
     */
    fn quad_box() -> Topology {
        let mut topol = Topology::with_capacity(8, 12, 6);
        let verts: Vec<_> = (0..8)
            .map(|_| topol.add_vertex().expect("Unable to add a vertex").index())
            .collect();
        assert_eq!(verts, (0u32..8).collect::<Vec<_>>());
        let mut cache = TopolCache::default();
        let faces: Vec<_> = [
            [0u32, 3, 2, 1],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
            [4, 5, 6, 7],
        ]
        .iter()
        .map(|indices| {
            topol
                .add_face(
                    &indices.iter().map(|idx| (*idx).into()).collect::<Vec<_>>(),
                    &mut cache,
                )
                .expect("Unable to add a face")
        })
        .collect();
        assert_eq!(faces, (0u32..6).map(|i| i.into()).collect::<Vec<_>>());
        assert_eq!(topol.num_vertices(), 8);
        assert_eq!(topol.num_halfedges(), 24);
        assert_eq!(topol.num_edges(), 12);
        assert_eq!(topol.num_faces(), 6);
        topol
    }

    fn loop_mesh() -> Topology {
        /*

                            12---------13---------14---------15
                           /          /          /          /
                          /   f5     /   f6     /    f7    /
                         /          /          /          /
                        /          /          /          /
                       8----------9----------10---------11
                      /          /          /          /
                     /    f3    /          /    f4    /
                    /          /          /          /
                   /          /          /          /
                  4----------5----------6----------7
                 /          /          /          /
                /   f0     /    f1    /    f2    /
               /          /          /          /
              /          /          /          /
             0----------1----------2----------3
        */
        let mut topol = Topology::with_capacity(16, 24, 8);
        let mut cache = TopolCache::default();
        for _ in 0u32..16 {
            topol.add_vertex().expect("Unable to add vertex");
        }
        for fvi in [
            [0u32, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [4, 5, 9, 8],
            [6, 7, 11, 10],
            [8, 9, 13, 12],
            [9, 10, 14, 13],
            [10, 11, 15, 14],
        ] {
            let vs = fvi.iter().map(|i| (*i).into()).collect::<ArrayVec<VH, 4>>();
            topol.add_face(&vs, &mut cache).expect("Unable to add face");
        }
        topol
    }

    #[test]
    fn t_box_vv_ccw_iter() {
        let qbox = quad_box();
        for (vi, vis) in [
            (0u32, [4u32, 3, 1]),
            (1u32, [2u32, 5, 0]),
            (2u32, [3u32, 6, 1]),
            (3u32, [0u32, 7, 2]),
            (4u32, [5u32, 7, 0]),
            (5u32, [6u32, 4, 1]),
            (6u32, [7u32, 5, 2]),
            (7u32, [4u32, 6, 3]),
        ] {
            assert_eq!(
                qbox.vv_ccw_iter(vi.into())
                    .map(|v| v.index())
                    .collect::<Vec<_>>(),
                vis
            );
        }
    }

    #[test]
    fn t_box_vv_cw_iter() {
        let qbox = quad_box();
        for (vi, fis) in [
            (0u32, [4, 1, 3]),
            (1u32, [2, 0, 5]),
            (2u32, [3, 1, 6]),
            (3u32, [0, 2, 7]),
            (4u32, [5, 0, 7]),
            (5u32, [6, 1, 4]),
            (6u32, [7, 2, 5]),
            (7u32, [4, 3, 6]),
        ] {
            assert_eq!(
                qbox.vv_cw_iter(vi.into())
                    .map(|x| x.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_box_vih_iter() {
        let qbox = quad_box();
        for v in qbox.vertices() {
            assert!(qbox
                .vih_ccw_iter(v)
                .all(|h| h.head(&qbox) == v && h.tail(&qbox) != v));
            assert!(qbox
                .vih_cw_iter(v)
                .all(|h| h.head(&qbox) == v && h.tail(&qbox) != v));
        }
    }

    #[test]
    fn t_box_voh_iter() {
        let qbox = quad_box();
        for v in qbox.vertices() {
            assert!(qbox
                .voh_ccw_iter(v)
                .all(|h| h.tail(&qbox) == v && h.head(&qbox) != v));
            assert!(qbox
                .voh_cw_iter(v)
                .all(|h| h.tail(&qbox) == v && h.head(&qbox) != v));
        }
    }

    #[test]
    fn t_box_vf_ccw_iter() {
        let qbox = quad_box();
        for (vi, fis) in [
            (0u32, [4u32, 0, 1]),
            (1u32, [2u32, 1, 0]),
            (2u32, [3u32, 2, 0]),
            (3u32, [4u32, 3, 0]),
            (4u32, [5u32, 4, 1]),
            (5u32, [5u32, 1, 2]),
            (6u32, [5u32, 2, 3]),
            (7u32, [5u32, 3, 4]),
        ] {
            assert_eq!(
                qbox.vf_ccw_iter(vi.into())
                    .map(|x| x.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_box_vf_cw_iter() {
        let qbox = quad_box();
        for (vi, fis) in [
            (0u32, [4u32, 1, 0]),
            (1, [2, 0, 1]),
            (2, [3, 0, 2]),
            (3, [4, 0, 3]),
            (4, [5, 1, 4]),
            (5, [5, 2, 1]),
            (6, [5, 3, 2]),
            (7, [5, 4, 3]),
        ] {
            assert_eq!(
                qbox.vf_cw_iter(vi.into())
                    .map(|x| x.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_box_fv_ccw_iter() {
        let qbox = quad_box();
        for (fi, vis) in [
            (0u32, [0, 3, 2, 1]),
            (1u32, [0, 1, 5, 4]),
            (2u32, [1, 2, 6, 5]),
            (3u32, [2, 3, 7, 6]),
            (4u32, [3, 0, 4, 7]),
            (5u32, [4, 5, 6, 7]),
        ] {
            assert_eq!(
                qbox.fv_ccw_iter(fi.into())
                    .map(|x| x.index())
                    .collect::<Vec<_>>(),
                vis
            );
        }
    }

    #[test]
    fn t_box_fv_cw_iter() {
        let qbox = quad_box();
        for (fi, vis) in [
            (0u32, [0, 1, 2, 3]),
            (1u32, [0, 4, 5, 1]),
            (2u32, [1, 5, 6, 2]),
            (3u32, [2, 6, 7, 3]),
            (4u32, [3, 7, 4, 0]),
            (5u32, [4, 7, 6, 5]),
        ] {
            assert_eq!(
                qbox.fv_cw_iter(fi.into())
                    .map(|x| x.index())
                    .collect::<Vec<_>>(),
                vis
            );
        }
    }

    #[test]
    fn t_box_ff_ccw_iter() {
        let qbox = quad_box();
        for (fi, fis) in [
            (0u32, [1, 4, 3, 2]),
            (1u32, [4, 0, 2, 5]),
            (2u32, [1, 0, 3, 5]),
            (3u32, [2, 0, 4, 5]),
            (4u32, [3, 0, 1, 5]),
            (5u32, [4, 1, 2, 3]),
        ] {
            assert_eq!(
                qbox.ff_ccw_iter(fi.into())
                    .map(|f| f.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_box_ff_cw_iter() {
        let qbox = quad_box();
        for (fi, fis) in [
            (0u32, [1, 2, 3, 4]),
            (1u32, [4, 5, 2, 0]),
            (2u32, [1, 5, 3, 0]),
            (3u32, [2, 5, 4, 0]),
            (4u32, [3, 5, 1, 0]),
            (5u32, [4, 3, 2, 1]),
        ] {
            assert_eq!(
                qbox.ff_cw_iter(fi.into())
                    .map(|f| f.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_loop_mesh_vf_ccw_iter() {
        let topol = loop_mesh();
        for (v, fis) in [
            (0u32, vec![0u32]),
            (1, vec![1, 0]),
            (2, vec![2, 1]),
            (3, vec![2]),
            (4, vec![0, 3]),
            (5, vec![3, 0, 1]),
            (6, vec![1, 2, 4]),
            (7, vec![4, 2]),
            (8, vec![3, 5]),
            (9, vec![6, 5, 3]),
            (10, vec![4, 7, 6]),
            (11, vec![7, 4]),
            (12, vec![5]),
            (13, vec![5, 6]),
            (14, vec![6, 7]),
            (15, vec![7]),
        ] {
            assert_eq!(
                topol
                    .vf_ccw_iter(v.into())
                    .map(|i| i.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_loop_mesh_vf_cw_iter() {
        let topol = loop_mesh();
        for (v, fis) in [
            (0u32, vec![0u32]),
            (1, vec![0, 1]),
            (2, vec![1, 2]),
            (3, vec![2]),
            (4, vec![3, 0]),
            (5, vec![1, 0, 3]),
            (6, vec![4, 2, 1]),
            (7, vec![2, 4]),
            (8, vec![5, 3]),
            (9, vec![3, 5, 6]),
            (10, vec![6, 7, 4]),
            (11, vec![4, 7]),
            (12, vec![5]),
            (13, vec![6, 5]),
            (14, vec![7, 6]),
            (15, vec![7]),
        ] {
            assert_eq!(
                topol
                    .vf_cw_iter(v.into())
                    .map(|i| i.index())
                    .collect::<Vec<_>>(),
                fis
            );
        }
    }

    #[test]
    fn t_box_mesh_voh_ccw_iter_mut() {
        // I have other tests checking the mutability of the mesh etc. This is
        // just to make sure the mutable iterator walks the same topology as the
        // immutable iterator.
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create a box");
        let expected = mesh
            .vertices()
            .flat_map(|v| mesh.voh_ccw_iter(v))
            .collect::<Vec<_>>();
        let mut actual = Vec::new();
        for v in mesh.vertices() {
            actual.extend(mesh.voh_ccw_iter_mut(v).map(|(_mesh, h)| h));
        }
        assert_eq!(expected, actual);
    }

    #[test]
    fn t_box_mesh_voh_cw_iter_mut() {
        // I have other tests checking the mutability of the mesh etc. This is
        // just to make sure the mutable iterator walks the same topology as the
        // immutable iterator.
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create a box");
        let expected = mesh
            .vertices()
            .flat_map(|v| mesh.voh_cw_iter(v))
            .collect::<Vec<_>>();
        let mut actual = Vec::new();
        for v in mesh.vertices() {
            actual.extend(mesh.voh_cw_iter_mut(v).map(|(_mesh, h)| h));
        }
        assert_eq!(expected, actual);
    }

    #[test]
    fn t_box_mesh_fh_ccw_iter_mut() {
        // I have other tests checking the mutability of the mesh etc. This is
        // just to make sure the mutable iterator walks the same topology as the
        // immutable iterator.
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create a box");
        let expected = mesh
            .faces()
            .flat_map(|v| mesh.fh_ccw_iter(v))
            .collect::<Vec<_>>();
        let mut actual = Vec::new();
        for v in mesh.faces() {
            actual.extend(mesh.fh_ccw_iter_mut(v).map(|(_mesh, h)| h));
        }
        assert_eq!(expected, actual);
    }

    #[test]
    fn t_box_mesh_fh_cw_iter_mut() {
        // I have other tests checking the mutability of the mesh etc. This is
        // just to make sure the mutable iterator walks the same topology as the
        // immutable iterator.
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create a box");
        let expected = mesh
            .faces()
            .flat_map(|v| mesh.fh_cw_iter(v))
            .collect::<Vec<_>>();
        let mut actual = Vec::new();
        for v in mesh.faces() {
            actual.extend(mesh.fh_cw_iter_mut(v).map(|(_mesh, h)| h));
        }
        assert_eq!(expected, actual);
    }

    #[test]
    fn t_box_mesh_voh_ccw_iter_mut_delete_faces() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create a box");
        // Checking to make sure I can modify the mesh while iterating over it's
        // elements.  The borrow checking and safety should still be enforced
        // because I can only modify the mesh via the mutable references that
        // the iteator provides.
        for v in mesh.vertices() {
            // using a random condition to delete some faces, to make sure I can
            // modify the mesh. If you try changing any of the `m` inside this
            // loop to `mesh`, the borrow checker should complain and the code
            // should not compile.
            for (m, h) in mesh.voh_ccw_iter_mut(v) {
                if let Some(f) = h.face(m) {
                    if (f.index() + h.index()) % 2 != 0 {
                        m.delete_face(f, true).expect("Cannot delete face");
                    }
                }
            }
        }
        mesh.garbage_collection()
            .expect("Garbage collection failed");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(2, mesh.num_faces());
        assert_eq!(8, mesh.num_edges());
        assert_eq!(16, mesh.num_halfedges());
        assert_eq!(8, mesh.num_vertices());
    }
}
