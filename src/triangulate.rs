use crate::{
    element::{FH, VH},
    iterator, Topology,
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
}

#[cfg(test)]
mod test {
    use crate::{element::Handle, topol::test::quad_box};

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
}
