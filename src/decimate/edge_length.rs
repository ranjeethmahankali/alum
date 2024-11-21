use crate::{DotProductAdaptor, Error, PolyMeshT, HH, VH};
use std::ops::{Mul, Sub};

use super::Decimater;

/// Decimater implementation that prioritizes the collapse of short edges.
pub struct EdgeLengthDecimater<const DIM: usize, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: PartialOrd,
{
    max_length: A::Scalar,
}

impl<const DIM: usize, A> EdgeLengthDecimater<DIM, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: PartialOrd + Mul<Output = A::Scalar>,
{
    /// Create a new EdgeLengthDecimater.
    ///
    /// Edges longer than `max_length` will not be collapsed.
    pub fn new(max_length: A::Scalar) -> Self {
        EdgeLengthDecimater {
            max_length: max_length * max_length,
        }
    }
}

impl<const DIM: usize, A> Decimater<PolyMeshT<DIM, A>> for EdgeLengthDecimater<DIM, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: PartialOrd + Mul<Output = A::Scalar>,
{
    type Cost = A::Scalar;

    /// Cost of collapsing an edge, i.e. the squared length of the edge.
    fn collapse_cost(&self, mesh: &PolyMeshT<DIM, A>, h: HH) -> Option<Self::Cost> {
        let points = mesh.points();
        let points = points.try_borrow().ok()?;
        let elen = mesh.calc_edge_length_sqr(h.edge(), &points);
        if elen < self.max_length {
            Some(elen)
        } else {
            None
        }
    }

    fn before_collapse(&mut self, _mesh: &PolyMeshT<DIM, A>, _h: HH) -> Result<(), Error> {
        Ok(()) // Do nothing.
    }

    fn after_collapse(&mut self, _mesh: &PolyMeshT<DIM, A>, _v: VH) -> Result<(), Error> {
        Ok(()) // Do nothing.
    }
}

#[cfg(test)]
mod test {
    use super::EdgeLengthDecimater;
    use crate::{obj::test::bunny_mesh, HasDecimation, HasIterators, HasTopology};

    #[test]
    fn t_bunny_edge_length_decimate() {
        let mut mesh = bunny_mesh();
        let (nv, ne, nf) = (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces());
        let nloops_before = {
            let mut visited = mesh.create_halfedge_prop(false);
            let mut visited = visited.try_borrow_mut().unwrap();
            mesh.halfedges().fold(0usize, |count, h| {
                if !h.is_boundary(&mesh) || visited[h] {
                    count
                } else {
                    for h2 in mesh.loop_ccw_iter(h) {
                        visited[h2] = true;
                    }
                    count + 1
                }
            })
        };
        let mut decimater = EdgeLengthDecimater::new(0.1);
        mesh.decimate(&mut decimater, 2000)
            .expect("Cannot decimate");
        mesh.garbage_collection()
            .expect("Failed to garbage collect");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(nv - 2000, mesh.num_vertices());
        assert!(mesh.num_edges() < ne);
        assert!(mesh.num_faces() < nf);
        let nloops_after = {
            let mut visited = mesh.create_halfedge_prop(false);
            let mut visited = visited.try_borrow_mut().unwrap();
            mesh.halfedges().fold(0usize, |count, h| {
                if !h.is_boundary(&mesh) || visited[h] {
                    count
                } else {
                    for h2 in mesh.loop_ccw_iter(h) {
                        visited[h2] = true;
                    }
                    count + 1
                }
            })
        };
        assert_eq!(nloops_before, nloops_after);
    }
}
