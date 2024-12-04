/*!
# Mesh Decimation

This module provides a framework for mesh decimation, and functions to decimate
meshes such as [`decimate`](HasDecimation::decimate),
[`decimate_to_vertex_count`](HasDecimation::decimate_to_vertex_count),
[`decimate_to_face_count`](HasDecimation::decimate_to_face_count).

The decimation functions take a `decimater` as an argument. This decimater must
implement the [`Decimater<MeshT>`] trait. It is responsible for determining the
validity and priority of edge collapses. It can also keep track of the changes
being made to the mesh during decimation.

Out of the box, this module also provides the following implementations of
[`Decimater<MeshT>`]:

- [`EdgeLengthDecimater`](edge_length::EdgeLengthDecimater) prioritzes collapsing short edges.

- [`QuadricDecimater`](quadric::QuadricDecimater) prioritizes collapsing edges
  that minimize the quadric error. This decimater can be initialized as one of
  two variants, decided by the enum [`QuadricType`](quadric::QuadricType). The
  probabilistic quadrics use [this
  implementation](https://github.com/Philip-Trettner/probabilistic-quadrics) by Philip Trettner.

## Writing You Own Decimater Implementation

- [`Decimater<MeshT>::collapse_cost`] must return the cost of a collapse. The
  collapses with low cost will be prioritized.

- [`Decimater<MeshT>::before_collapse`], [`Decimater<MeshT>::after_collapse`]
  functions can be used to keep track of the mesh, update the locations of
  vertices, normals etc.

*/

pub mod edge_length;
pub mod quadric;

use crate::Queue;
use crate::{
    topol::Topology, Adaptor, EPropRef, EditableTopology, Error, HasIterators, PolyMeshT, Status,
    VPropRefMut, HH, VH,
};
use std::cmp::Ordering;

/// An implementation of this trait is required to control the mesh decimation.
///
/// You can either write your own implementation of this trait, or use one of
/// the implementations that come with this crate out of the box:
///
/// - [`Decimater<MeshT>::collapse_cost`] must return the cost of a collapse. The
///   collapses with low cost will be prioritized.
///
/// - [`Decimater<MeshT>::before_collapse`], [`Decimater<MeshT>::after_collapse`]
///   functions can be used to keep track of the mesh, update the locations of
///   vertices, normals etc.
pub trait Decimater<MeshT>
where
    MeshT: HasIterators,
    Self::Cost: PartialOrd,
{
    /// This is the type of cost returned by this decimater implementation. This
    /// must implement [`PartialOrd`] to allow for comparison, to prioritze edge
    /// collapses with the lowest cost.
    type Cost;

    /// The cost of a halfedge collapse.
    ///
    /// The cost of collapsing the tail of `h` toward the head of `h` must be
    /// returned. If `None` is returned, it means this collapse is not allowed.
    fn collapse_cost(&self, mesh: &MeshT, h: HH) -> Option<Self::Cost>;

    /// This can be used to keep track of the changes being done to the mesh.
    ///
    /// `h` is the edge that is about to be collapsed.
    fn before_collapse(&mut self, mesh: &MeshT, h: HH) -> Result<(), Error>;

    /// This can be used to keep track of the changes being done to the mesh.
    ///
    /// `v` is the vertex resulting from the edge that was just collapsed.
    fn after_collapse(&mut self, mesh: &MeshT, v: VH) -> Result<(), Error>;
}

/// Find the best collapse target for this vertex, queue up the vertex with the
/// associated cost, and return the target halfedge.
fn queue_vertex_collapse<MeshT, DecT>(
    mesh: &MeshT,
    v: VH,
    decimater: &DecT,
    heap: &mut Queue<VH, DecT::Cost>,
) -> Option<HH>
where
    MeshT: HasIterators,
    DecT: Decimater<MeshT>,
{
    // Find the best collapsible edge around this vertex, it's collapse cost and update the queue.
    match mesh.voh_ccw_iter(v).fold(None, |best, h| {
        match (best, decimater.collapse_cost(mesh, h)) {
            (None, None) => None,                                   // Found nothing.
            (None, Some(cost)) => Some((h, cost)),                  // Current wins by default.
            (Some((hbest, lowest)), None) => Some((hbest, lowest)), // Previous wins by default.
            (Some((hbest, lowest)), Some(cost)) => match cost.partial_cmp(&lowest) {
                Some(Ordering::Less) => Some((h, cost)), // Current wins by comparison
                _ => Some((hbest, lowest)),              // Previous wins.
            },
        }
    }) {
        Some((h, cost)) => {
            heap.insert(v, cost);
            Some(h)
        }
        None => {
            heap.remove(v);
            None
        }
    }
}

fn is_collapse_legal(
    mesh: &Topology,
    h: HH,
    estatus: &EPropRef<Status>,
    vstatus: &mut VPropRefMut<Status>,
) -> bool {
    let v = h.tail(mesh);
    !vstatus[v].feature() // Vertex not a feature.
        && !estatus[h.edge()].feature() // Edge not a feature.
        && mesh.vf_ccw_iter(v).take(2).count() == 2 // Has at leaset two incident faces.
        && mesh.check_edge_collapse(h, estatus, vstatus) // No topological errors.
}

/// Mesh types can support decimation by implementing this trait.
///
/// The built in types implement this trait.
pub trait HasDecimation: EditableTopology {
    /// Decimate the mesh by collapsing the halfedge, until either we run out of
    /// halfedges to collapse, or if the given predicate returns `false`.
    ///
    /// The predicate will be called with three arguments: number of collapses
    /// done so far, number of vertices left in the mesh, and number of faces
    /// left in the mesh.
    fn decimate_while<F, DecT>(&mut self, decimater: &mut DecT, pred: F) -> Result<usize, Error>
    where
        F: Fn(usize, usize, usize) -> bool,
        DecT: Decimater<Self>,
    {
        let mut collapse_targets = self.create_vertex_prop::<Option<HH>>(None);
        let mut collapse_targets = collapse_targets.try_borrow_mut()?;
        let mut cache = Vec::new();
        let mut vstatus = self.vertex_status_prop();
        let mut hstatus = self.halfedge_status_prop();
        let mut estatus = self.edge_status_prop();
        let mut fstatus = self.face_status_prop();
        // Initialize heap with vertex collapses.
        let mut heap = Queue::<VH, DecT::Cost>::new(self.num_vertices());
        {
            let vstatus = vstatus.try_borrow()?;
            for (v, target) in self
                .vertices()
                .zip(collapse_targets.iter_mut())
                .filter(|(v, _target)| !vstatus[*v].deleted())
            {
                *target = queue_vertex_collapse(self, v, decimater, &mut heap);
            }
        }
        let (mut nc, mut nv, mut nf) = (0usize, self.num_vertices(), self.num_faces());
        let mut one_ring = Vec::new();
        while let Some((v0, _cost)) = heap.pop() {
            if !pred(nc, nv, nf) {
                break;
            }
            // Get the collapse target and check if it is legal.
            let h = collapse_targets[v0].ok_or(Error::UndefinedCollapseTarget(v0))?;
            {
                let estatus = estatus.try_borrow()?;
                let mut vstatus = vstatus.try_borrow_mut()?;
                if !is_collapse_legal(self.topology(), h, &estatus, &mut vstatus) {
                    continue;
                }
            }
            one_ring.clear(); // Collect the one ring for later.
            one_ring.extend(self.vv_ccw_iter(v0));
            // Update stats.
            (nc, nv, nf) = (
                nc + 1,
                nv - 1,
                nf - if h.edge().is_boundary(self) { 1 } else { 2 },
            );
            let v1 = h.head(self);
            decimater.before_collapse(self, h)?;
            {
                let mut vstatus = vstatus.try_borrow_mut()?;
                let mut hstatus = hstatus.try_borrow_mut()?;
                let mut estatus = estatus.try_borrow_mut()?;
                let mut fstatus = fstatus.try_borrow_mut()?;
                self.collapse_edge(
                    h,
                    &mut vstatus,
                    &mut hstatus,
                    &mut estatus,
                    &mut fstatus,
                    &mut cache,
                );
            }
            decimater.after_collapse(self, v1)?;
            // Update heap.
            for &v in one_ring.iter() {
                collapse_targets[v] = queue_vertex_collapse(self, v, decimater, &mut heap);
            }
        }
        heap.clear();
        Ok(nc)
    }

    /// Decimate the mesh by collapsing the given number of edges.
    ///
    /// The `decimater` controls the priority of the edges to be collapsed, and
    /// can also keep track of the changes being made to the mesh. Read [this
    /// explanation](crate::decimate#mesh-decimation) for more on how to use the
    /// decimater.
    fn decimate(
        &mut self,
        decimater: &mut impl Decimater<Self>,
        num_collapses: usize,
    ) -> Result<usize, Error> {
        self.decimate_while(decimater, |n, _v, _f| n < num_collapses)
    }

    /// Try to decimate the mesh till the target vertex count is achieved.
    ///
    /// The decimation may stop sooner due to other criteria. The `decimater`
    /// controls the priority of the edges to be collapsed, and can also keep
    /// track of the changes being made to the mesh. Read [this
    /// explanation](crate::decimate#mesh-decimation) for more on how to use the
    /// decimater.
    fn decimate_to_vertex_count(
        &mut self,
        decimater: &mut impl Decimater<Self>,
        vert_target: usize,
    ) -> Result<usize, Error> {
        let vtarget = usize::min(self.num_vertices(), vert_target);
        self.decimate_while(decimater, |_n, v, _f| v > vtarget)
    }

    /// Try to decimate the mesh till the target mesh count is achieved.
    ///
    /// The decimation may stop sooner due to other criteria. The `decimater`
    /// controls the priority of the edges to be collapsed, and can also keep
    /// track of the changes being made to the mesh. Read [this
    /// explanation](crate::decimate#mesh-decimation) for more on how to use the
    /// decimater.
    fn decimate_to_face_count(
        &mut self,
        decimater: &mut impl Decimater<Self>,
        face_target: usize,
    ) -> Result<usize, Error> {
        let ftarget = usize::min(face_target, self.num_faces());
        self.decimate_while(decimater, |_n, _v, f| f > ftarget)
    }
}

impl<const DIM: usize, A> HasDecimation for PolyMeshT<DIM, A> where A: Adaptor<DIM> {}

impl HasDecimation for Topology {}
