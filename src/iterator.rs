use crate::{
    element::{EH, FH, HH, VH},
    topol::Topology,
    Adaptor, PolyMeshT,
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

/* Immutable iterators. */

pub(crate) fn vv_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = VH> + use<'_> {
    voh_ccw_iter(topol, v).map(|h| h.head(topol))
}

pub(crate) fn vv_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = VH> + use<'_> {
    voh_cw_iter(topol, v).map(|h| h.head(topol))
}

pub(crate) fn vih_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    voh_ccw_iter(topol, v).map(|h| h.opposite())
}

pub(crate) fn vih_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    voh_cw_iter(topol, v).map(|h| h.opposite())
}

pub(crate) fn voh_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    RadialHalfedgeIter::<true>::new(topol, v.halfedge(topol))
}

pub(crate) fn voh_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    RadialHalfedgeIter::<false>::new(topol, v.halfedge(topol))
}

pub(crate) fn ve_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = EH> + use<'_> {
    voh_ccw_iter(topol, v).map(|h| h.edge())
}

pub(crate) fn ve_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = EH> + use<'_> {
    voh_cw_iter(topol, v).map(|h| h.edge())
}

pub(crate) fn vf_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = FH> + use<'_> {
    voh_ccw_iter(topol, v).filter_map(|h| h.face(topol))
}

pub(crate) fn vf_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = FH> + use<'_> {
    voh_cw_iter(topol, v).filter_map(|h| h.face(topol))
}

pub(crate) fn ev_iter(topol: &Topology, e: EH) -> impl Iterator<Item = VH> + use<'_> {
    eh_iter(e).map(|h| h.head(topol))
}

pub(crate) fn eh_iter(e: EH) -> impl Iterator<Item = HH> {
    [false, true].iter().map(move |flag| e.halfedge(*flag))
}

pub(crate) fn ef_iter(topol: &Topology, e: EH) -> impl Iterator<Item = FH> + use<'_> {
    eh_iter(e).filter_map(|h| h.face(topol))
}

pub(crate) fn fv_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = VH> + use<'_> {
    fh_ccw_iter(topol, f).map(|h| h.head(topol))
}

pub(crate) fn fv_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = VH> + use<'_> {
    fh_cw_iter(topol, f).map(|h| h.head(topol))
}

pub(crate) fn fh_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = HH> + use<'_> {
    LoopHalfedgeIter::<true>::new(topol, f.halfedge(topol))
}

pub(crate) fn fh_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = HH> + use<'_> {
    LoopHalfedgeIter::<false>::new(topol, f.halfedge(topol))
}

pub(crate) fn fe_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = EH> + use<'_> {
    fh_ccw_iter(topol, f).map(|h| h.edge())
}

pub(crate) fn fe_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = EH> + use<'_> {
    fh_cw_iter(topol, f).map(|h| h.edge())
}

pub(crate) fn ff_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = FH> + use<'_> {
    fh_ccw_iter(topol, f).filter_map(|h| h.opposite().face(topol))
}

pub(crate) fn ff_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = FH> + use<'_> {
    fh_cw_iter(topol, f).filter_map(|h| h.opposite().face(topol))
}

pub(crate) fn ccw_rotate_iter(topol: &Topology, h: HH) -> impl Iterator<Item = HH> + use<'_> {
    RadialHalfedgeIter::<true>::new(topol, Some(h))
}

pub(crate) fn cw_rotate_iter(topol: &Topology, h: HH) -> impl Iterator<Item = HH> + use<'_> {
    RadialHalfedgeIter::<false>::new(topol, Some(h))
}

pub(crate) fn loop_ccw_iter(topol: &Topology, h: HH) -> impl Iterator<Item = HH> + use<'_> {
    LoopHalfedgeIter::<true>::new(topol, h)
}

pub(crate) fn loop_cw_iter(topol: &Topology, h: HH) -> impl Iterator<Item = HH> + use<'_> {
    LoopHalfedgeIter::<false>::new(topol, h)
}

/* Mutable iterators.*/

pub(crate) fn loop_ccw_iter_mut(
    topol: &mut Topology,
    h: HH,
) -> impl Iterator<Item = (&mut Topology, HH)> + use<'_> {
    LoopHalfedgeIterMut::<true, Topology> {
        reference: topol.into(),
        topol,
        hstart: h,
        hcurrent: Some(h),
        _phantom: PhantomData,
    }
}

pub(crate) fn fh_ccw_iter_mut(
    topol: &mut Topology,
    f: FH,
) -> impl Iterator<Item = (&mut Topology, HH)> + use<'_> {
    let h = f.halfedge(topol);
    LoopHalfedgeIterMut::<true, Topology> {
        reference: topol.into(),
        topol,
        hstart: h,
        hcurrent: Some(h),
        _phantom: PhantomData,
    }
}

// Expose immutable iterators as member functions of the mesh.
impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    /// Iterator over the outgoing halfedges around a vertex, going counter-clockwise.
    pub fn voh_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        voh_ccw_iter(&self.topol, v)
    }

    /// Iterator over the outgoing halfedges around a vertex, going clockwise
    pub fn voh_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        voh_cw_iter(&self.topol, v)
    }

    /// Iterator over the incoming halfedges around a vertex, going
    /// counter-clockwise
    pub fn vih_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        vih_ccw_iter(&self.topol, v)
    }

    /// Iterator over the incoming halfedges around a vertex, going clockwise
    pub fn vih_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        vih_cw_iter(&self.topol, v)
    }

    /// Iterator over the faces incident on a vertex, going counter-clockwise.
    pub fn vf_ccw_iter(&self, v: VH) -> impl Iterator<Item = FH> + use<'_, A, DIM> {
        vf_ccw_iter(&self.topol, v)
    }

    /// Iterator over the faces incident on a vertex, going clockwise.
    pub fn vf_cw_iter(&self, v: VH) -> impl Iterator<Item = FH> + use<'_, A, DIM> {
        vf_cw_iter(&self.topol, v)
    }

    /// Iterator over the neighboring vertices around the given vertex, going
    /// counter-clockwise.
    pub fn vv_ccw_iter(&self, v: VH) -> impl Iterator<Item = VH> + use<'_, A, DIM> {
        vv_ccw_iter(&self.topol, v)
    }

    /// Iterator over the neighboring vertices around the given vertex, going
    /// clockwise.
    pub fn vv_cw_iter(&self, v: VH) -> impl Iterator<Item = VH> + use<'_, A, DIM> {
        vv_cw_iter(&self.topol, v)
    }

    /// Iterator over the incident edges around an vertex, going counter-clockwise.
    pub fn ve_ccw_iter(&self, v: VH) -> impl Iterator<Item = EH> + use<'_, A, DIM> {
        ve_ccw_iter(&self.topol, v)
    }

    /// Iterator over the incident edges around a vertex, going clockwise.
    pub fn ve_cw_iter(&self, v: VH) -> impl Iterator<Item = EH> + use<'_, A, DIM> {
        ve_cw_iter(&self.topol, v)
    }

    /// Iterator over the two vertices incident on the given edge.
    pub fn ev_iter(&self, e: EH) -> impl Iterator<Item = VH> + use<'_, A, DIM> {
        ev_iter(&self.topol, e)
    }

    /// Iterator over the two halfedges corresponding to an edge.
    pub fn eh_iter(&self, e: EH) -> impl Iterator<Item = HH> + use<A, DIM> {
        eh_iter(e)
    }

    /// Iterator over the faces incident on an edge.
    pub fn ef_iter(&self, e: EH) -> impl Iterator<Item = FH> + use<'_, A, DIM> {
        ef_iter(&self.topol, e)
    }

    /// Iterator over the halfedges of a face loop, going counter-clockwise.
    pub fn fh_ccw_iter(&self, f: FH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        fh_ccw_iter(&self.topol, f)
    }

    /// Iterator over the halfedges of a face loop, going clockwise.
    pub fn fh_cw_iter(&self, f: FH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        fh_cw_iter(&self.topol, f)
    }

    /// Iterator over the vertices incident on a face, going counter-clockwise.
    pub fn fv_ccw_iter(&self, f: FH) -> impl Iterator<Item = VH> + use<'_, A, DIM> {
        fv_ccw_iter(&self.topol, f)
    }

    /// Iterator over the vertices incident on a face, going clockwise.
    pub fn fv_cw_iter(&self, f: FH) -> impl Iterator<Item = VH> + use<'_, A, DIM> {
        fv_cw_iter(&self.topol, f)
    }

    /// Iterator over the edges incident on a face, going counter-clockwise.
    pub fn fe_ccw_iter(&self, f: FH) -> impl Iterator<Item = EH> + use<'_, A, DIM> {
        fe_ccw_iter(&self.topol, f)
    }

    /// Iterator over the edges incident on a face, going clockwise.
    pub fn fe_cw_iter(&self, f: FH) -> impl Iterator<Item = EH> + use<'_, A, DIM> {
        fe_cw_iter(&self.topol, f)
    }

    /// Iterator over the neighboring faces arouund the given face, going
    /// counter-clockwise.
    ///
    /// This includes the faces connected via a shared edge, but not those
    /// connected via a shared vertex.
    pub fn ff_ccw_iter(&self, f: FH) -> impl Iterator<Item = FH> + use<'_, A, DIM> {
        ff_ccw_iter(&self.topol, f)
    }

    /// Iterator over the neighboring faces around the given face, going
    /// clockwise.
    ///
    /// This includes the faces connected via a shared edge, but not those
    /// connected via a shared vertex.
    pub fn ff_cw_iter(&self, f: FH) -> impl Iterator<Item = FH> + use<'_, A, DIM> {
        ff_cw_iter(&self.topol, f)
    }

    /// This is similar to [`Self::voh_ccw_iter`] around the tail of the given
    /// halfedge, except this iterator starts at the provided halfedge.
    ///
    /// This is equivalent to a circular shifted [`Self::voh_ccw_iter`] of the
    /// vertex at the tail of this halfedge.
    pub fn ccw_rotate_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        ccw_rotate_iter(&self.topol, h)
    }

    /// This is similar to [`Self::voh_cw_iter`] around the tail of the given
    /// halfedge, except this iterator starts at the provided halfedge.
    ///
    /// This is equivalent to a circular shifted [`Self::voh_cw_iter`] of the
    /// vertex at the tail of this halfedge.
    pub fn cw_rotate_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        cw_rotate_iter(&self.topol, h)
    }

    /// Counter-clockwise iterator over the halfedges in a loop.
    ///
    /// The iterator will start at the given halfedge. If the halfedge has an
    /// incident face, this iterator is equivalent to a circular shifted
    /// [`Self::fh_ccw_iter`] of the incident face. If the halfedge is on the
    /// boundary, this iterator goes over the boundary loop counter-clockwise.
    pub fn loop_ccw_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        loop_ccw_iter(&self.topol, h)
    }

    /// Counter-clockwise iterator over the halfedges in a loop.
    ///
    /// The iterator will start at the given halfedge. If the halfedge has an
    /// incident face, this iterator is equivalent to a circular shifted
    /// [`Self::fh_cw_iter`] of the incident face. If the halfedge is on the
    /// boundary, this iterator goes over the boundary loop clockwise.
    pub fn loop_cw_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, A, DIM> {
        loop_cw_iter(&self.topol, h)
    }
}

// Expose mutable iterators as member functions.
impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
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
    use alum::{alum_glam::PolyMeshF32, FH};

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

    /// Mutable iterator over the neighboring vertices around the given vertex,
    /// going counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn vv_ccw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, VH)> + use<'_, A, DIM> {
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
    pub fn vv_cw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, VH)> + use<'_, A, DIM> {
        self.voh_cw_iter_mut(v).map(|(mesh, h)| {
            let v = h.head(mesh);
            (mesh, v)
        })
    }

    /// Iterator over the incoming halfedges around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn vih_ccw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        self.voh_ccw_iter_mut(v)
            .map(|(mesh, h)| (mesh, h.opposite()))
    }

    /// Iterator over the incoming halfedges around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn vih_cw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        self.voh_cw_iter_mut(v)
            .map(|(mesh, h)| (mesh, h.opposite()))
    }

    /// Iterator over the outgoing halfedges around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn voh_ccw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        let h = v.halfedge(&self.topol);
        RadialHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: &mut self.topol,
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
    pub fn voh_cw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        let h = v.halfedge(&self.topol);
        RadialHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: &mut self.topol,
            hstart: h,
            hcurrent: h,
            _phantom: PhantomData,
        }
    }

    /// Iterator over the incident edges around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn ve_ccw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, EH)> + use<'_, A, DIM> {
        self.voh_ccw_iter_mut(v).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the incident edges around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn ve_cw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, EH)> + use<'_, A, DIM> {
        self.voh_cw_iter_mut(v).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the incident faces around a vertex, going
    /// counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn vf_ccw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, FH)> + use<'_, A, DIM> {
        self.voh_ccw_iter_mut(v)
            .filter_map(|(mesh, h)| h.face(mesh).map(|f| (mesh, f)))
    }

    /// Iterator over the incident faces around a vertex, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn vf_cw_iter_mut(
        &mut self,
        v: VH,
    ) -> impl Iterator<Item = (&mut Self, FH)> + use<'_, A, DIM> {
        self.voh_cw_iter_mut(v)
            .filter_map(|(mesh, h)| h.face(mesh).map(|f| (mesh, f)))
    }

    /// Iterator over the vertices incident on a face, going counter-clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn fv_ccw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, VH)> + use<'_, A, DIM> {
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
    pub fn fv_cw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, VH)> + use<'_, A, DIM> {
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
    pub fn fh_ccw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        let h = f.halfedge(self);
        LoopHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: &mut self.topol,
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
    pub fn fh_cw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        let h = f.halfedge(self);
        LoopHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: &mut self.topol,
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
    pub fn fe_ccw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, EH)> + use<'_, A, DIM> {
        self.fh_ccw_iter_mut(f).map(|(mesh, h)| (mesh, h.edge()))
    }

    /// Iterator over the edges incident on a face, going clockwise.
    ///
    /// A mutable iterator allows modifying the mesh while iterating over its
    /// elements. See [this section](#working-with-mutable-iterators) for an
    /// overview on mutable iterators.
    pub fn fe_cw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, EH)> + use<'_, A, DIM> {
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
    pub fn ff_ccw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, FH)> + use<'_, A, DIM> {
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
    pub fn ff_cw_iter_mut(
        &mut self,
        f: FH,
    ) -> impl Iterator<Item = (&mut Self, FH)> + use<'_, A, DIM> {
        self.fh_cw_iter_mut(f)
            .filter_map(|(mesh, h)| h.opposite().face(mesh).map(|f| (mesh, f)))
    }

    /// This is similar to [`Self::voh_ccw_iter_mut`] around the tail of the
    /// given halfedge, except this iterator starts at the provided halfedge.
    ///
    /// Thisis equivalent to a circular shifted [`Self::voh_ccw_iter_mut`] of
    /// the vertex at the tail of the given halfedge.
    pub fn ccw_rotate_iter_mut(
        &mut self,
        h: HH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        RadialHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: &mut self.topol,
            hstart: Some(h),
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }

    /// This is similar to [`Self::voh_cw_iter_mut`] around the tail of the
    /// given halfedge, except this iterator starts at the provided halfedge.
    ///
    /// Thisis equivalent to a circular shifted [`Self::voh_cw_iter_mut`] of
    /// the vertex at the tail of the given halfedge.
    pub fn cw_rotate_iter_mut(
        &mut self,
        h: HH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        RadialHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: &mut self.topol,
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
    /// counter-clockwise.
    pub fn loop_ccw_iter_mut(
        &mut self,
        h: HH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        LoopHalfedgeIterMut::<true, Self> {
            reference: self.into(),
            topol: &mut self.topol,
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
    /// boundary, this iterator goes over the boundary loop clockwise.
    pub fn loop_cw_iter_mut(
        &mut self,
        h: HH,
    ) -> impl Iterator<Item = (&mut Self, HH)> + use<'_, A, DIM> {
        LoopHalfedgeIterMut::<false, Self> {
            reference: self.into(),
            topol: &mut self.topol,
            hstart: h,
            hcurrent: Some(h),
            _phantom: PhantomData,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{
        alum_glam::PolyMeshF32,
        element::{Handle, VH},
        iterator::{
            ff_ccw_iter, ff_cw_iter, fv_ccw_iter, fv_cw_iter, vf_ccw_iter, vf_cw_iter,
            vih_ccw_iter, vih_cw_iter, voh_ccw_iter, voh_cw_iter, vv_ccw_iter, vv_cw_iter,
        },
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
                vv_ccw_iter(&qbox, vi.into())
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
                vv_cw_iter(&qbox, vi.into())
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
            assert!(vih_ccw_iter(&qbox, v).all(|h| h.head(&qbox) == v && h.tail(&qbox) != v));
            assert!(vih_cw_iter(&qbox, v).all(|h| h.head(&qbox) == v && h.tail(&qbox) != v));
        }
    }

    #[test]
    fn t_box_voh_iter() {
        let qbox = quad_box();
        for v in qbox.vertices() {
            assert!(voh_ccw_iter(&qbox, v).all(|h| h.tail(&qbox) == v && h.head(&qbox) != v));
            assert!(voh_cw_iter(&qbox, v).all(|h| h.tail(&qbox) == v && h.head(&qbox) != v));
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
                vf_ccw_iter(&qbox, vi.into())
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
                vf_cw_iter(&qbox, vi.into())
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
                fv_ccw_iter(&qbox, fi.into())
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
                fv_cw_iter(&qbox, fi.into())
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
                ff_ccw_iter(&qbox, fi.into())
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
                ff_cw_iter(&qbox, fi.into())
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
                vf_ccw_iter(&topol, v.into())
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
                vf_cw_iter(&topol, v.into())
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
