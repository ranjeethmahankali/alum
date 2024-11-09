use crate::{
    element::{FH, VH},
    error::Error,
    iterator,
    mesh::PolyMeshT,
    vector::TVec,
    Topology,
};

impl Topology {
    pub fn triangulated_face_vertices(&self, f: FH) -> impl Iterator<Item = [VH; 3]> + use<'_> {
        let hstart = self.face_halfedge(f);
        let vstart = self.from_vertex(hstart);
        iterator::loop_ccw_iter(self, self.next_halfedge(hstart))
            .take_while(move |h| self.to_vertex(*h) != vstart)
            .map(move |h| [vstart, self.from_vertex(h), self.to_vertex(h)])
    }

    pub fn triangulated_vertices(&self) -> impl Iterator<Item = [VH; 3]> + use<'_> {
        self.faces()
            .flat_map(move |f| self.triangulated_face_vertices(f))
    }

    pub fn triangulate_face(&mut self, f: FH) -> Result<(), Error> {
        let mut base = self.face_halfedge(f);
        let vstart = self.from_vertex(base);
        let prev = self.prev_halfedge(base);
        let mut next = self.next_halfedge(base);
        while self.to_vertex(self.next_halfedge(next)) != vstart {
            let next2 = self.next_halfedge(next);
            let fnew = self.new_face(base)?;
            let enew = self.new_edge(vstart, self.to_vertex(next), prev, next2, next, base)?;
            let hnew = self.edge_halfedge(enew, false);
            let ohnew = self.edge_halfedge(enew, true);
            // Link the triangle created.
            self.link_halfedges(base, next);
            self.link_halfedges(next, ohnew);
            self.link_halfedges(ohnew, base);
            // Set face handles.
            self.halfedge_mut(base).face = Some(fnew);
            self.halfedge_mut(next).face = Some(fnew);
            self.halfedge_mut(ohnew).face = Some(fnew);
            // Copy properties.
            self.hprops.copy(prev, ohnew)?;
            self.hprops.copy(prev, hnew)?;
            self.fprops.copy(f, fnew)?;
            // For next iteration.
            base = hnew;
            next = next2;
        }
        // Last face takes the original face handle.
        self.face_mut(f).halfedge = base;
        self.link_halfedges(base, next);
        self.link_halfedges(self.next_halfedge(next), base);
        self.halfedge_mut(base).face = Some(f);
        Ok(())
    }

    pub fn triangulate(&mut self) -> Result<(), Error> {
        for f in self.faces() {
            self.triangulate_face(f)?;
        }
        Ok(())
    }
}

impl<VecT, const DIM: usize> PolyMeshT<VecT, DIM>
where
    VecT: TVec<DIM>,
{
    pub fn triangulate_face(&mut self, f: FH) -> Result<(), Error> {
        self.topology_mut().triangulate_face(f)
    }

    pub fn triangulate(&mut self) -> Result<(), Error> {
        self.topology_mut().triangulate()
    }
}

#[cfg(test)]
mod test {
    use crate::{element::Handle, iterator, topol::test::quad_box};

    #[test]
    fn t_box_triangulated_indices() {
        let qbox = quad_box();
        assert_eq!(
            qbox.triangulated_vertices()
                .flatten()
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[
                1, 0, 3, 1, 3, 2, 4, 0, 1, 4, 1, 5, 5, 1, 2, 5, 2, 6, 6, 2, 3, 6, 3, 7, 7, 3, 0, 7,
                0, 4, 7, 4, 5, 7, 5, 6
            ]
        );
    }

    #[test]
    fn t_box_triangulate_face() {
        let qbox = {
            let mut mesh = quad_box();
            mesh.triangulate_face(5.into())
                .expect("Failed to triangulate face");
            mesh
        };
        assert_eq!(7, qbox.num_faces());
        assert_eq!(
            (2, 5),
            qbox.faces().fold((0usize, 0usize), |(t, q), f| {
                match qbox.face_valence(f) {
                    3 => (t + 1, q),
                    4 => (t, 1 + q),
                    _ => (t, q),
                }
            })
        );
        assert_eq!(13, qbox.num_edges());
        assert_eq!(26, qbox.num_halfedges());
        assert_eq!(
            iterator::fv_ccw_iter(&qbox, 5.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[5, 6, 7]
        );
        assert_eq!(
            iterator::fv_ccw_iter(&qbox, 6.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[4, 5, 7]
        )
    }

    #[test]
    fn t_box_triangulate() {
        let qbox = {
            let mut mesh = quad_box();
            mesh.triangulate().expect("Cannot triangulate mesh");
            mesh
        };
        assert_eq!(12, qbox.num_faces());
        assert_eq!(18, qbox.num_edges());
        assert_eq!(36, qbox.num_halfedges());
        assert_eq!(8, qbox.num_vertices());
        assert_eq!(
            qbox.faces()
                .flat_map(|f| iterator::fv_ccw_iter(&qbox, f))
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[
                3, 2, 1, 1, 5, 4, 2, 6, 5, 3, 7, 6, 0, 4, 7, 5, 6, 7, 0, 3, 1, 0, 1, 4, 1, 2, 5, 2,
                3, 6, 3, 0, 7, 4, 5, 7
            ]
        );
    }
}
