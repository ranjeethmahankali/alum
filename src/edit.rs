use crate::{
    element::{Handle, FH, HH, VH},
    error::Error,
    iterator,
    mesh::PolyMeshT,
    status::Status,
    topol::TopolCache,
    vector::TVec,
    Topology,
};

impl Topology {
    pub fn check_edge_collapse(
        &self,
        h: HH,
        edge_status: &[Status],
        vertex_status: &mut [Status],
    ) -> bool {
        // Check if already deleted.
        if edge_status[self.halfedge_edge(h).index() as usize].deleted() {
            return false;
        }
        let oh = self.opposite_halfedge(h);
        let v0 = self.to_vertex(oh);
        let v1 = self.to_vertex(h);
        if vertex_status[v0.index() as usize].deleted()
            || vertex_status[v1.index() as usize].deleted()
        {
            return false;
        }
        let htriangle = match self.halfedge_face(h) {
            Some(f) => self.face_valence(f) == 3,
            None => false,
        };
        let ohtriangle = match self.halfedge_face(oh) {
            Some(f) => self.face_valence(f) == 3,
            None => false,
        };
        // Check if the faces are triangles and the vertices opposite the edge
        // on those triangles are actually the same vertex.
        let vl = if htriangle {
            let h1 = self.next_halfedge(h);
            let h2 = self.next_halfedge(h1);
            if self.is_boundary_halfedge(self.opposite_halfedge(h1))
                && self.is_boundary_halfedge(self.opposite_halfedge(h2))
            {
                return false;
            }
            Some(self.to_vertex(h1))
        } else {
            None
        };
        let vr = if ohtriangle {
            let h1 = self.next_halfedge(oh);
            let h2 = self.next_halfedge(h1);
            if self.is_boundary_halfedge(self.opposite_halfedge(h1))
                && self.is_boundary_halfedge(self.opposite_halfedge(h2))
            {
                return false;
            }
            Some(self.to_vertex(h1))
        } else {
            None
        };
        if let (Some(vl), Some(vr)) = (vl, vr) {
            if vl == vr {
                return false;
            }
        }
        // Check if we're collapsing across two different boundaries.
        if self.is_boundary_vertex(v0)
            && self.is_boundary_vertex(v1)
            && !self.is_boundary_halfedge(h)
            && !self.is_boundary_halfedge(oh)
        {
            return false;
        }
        // Check the 'Link condition' Edelsbrunner [2006]. The intersection of
        // the one rings of from and to vertices must be the left and right
        // vertices, and only if the corresponding faces are triangles.
        let vl = self.to_vertex(self.next_halfedge(h));
        let vr = self.to_vertex(self.next_halfedge(oh));
        for v in iterator::vv_ccw_iter(self, v0) {
            vertex_status[v.index() as usize].set_tagged(false);
        }
        for v in iterator::vv_ccw_iter(self, v1) {
            vertex_status[v.index() as usize].set_tagged(true);
        }
        for v in iterator::vv_ccw_iter(self, v0) {
            if vertex_status[v.index() as usize].tagged()
                && !(v == vl && htriangle)
                && !(v == vr && ohtriangle)
            {
                return false;
            }
        }
        // Check for folded faces that might degenerate.
        if htriangle {
            let h1 = self.opposite_halfedge(self.next_halfedge(h));
            let h2 = self.opposite_halfedge(self.prev_halfedge(h));
            match (self.halfedge_face(h1), self.halfedge_face(h2)) {
                (None, None) => return false, // This is redundant but just in case.
                (Some(fa), Some(fb)) if fa == fb && self.face_valence(fa) != 3 => return false,
                _ => {} // Do nothing.
            }
        }
        if ohtriangle {
            let h1 = self.opposite_halfedge(self.next_halfedge(oh));
            let h2 = self.opposite_halfedge(self.prev_halfedge(oh));
            match (self.halfedge_face(h1), self.halfedge_face(h2)) {
                (None, None) => return false, // This is redundant but just in case.
                (Some(fa), Some(fb)) if fa == fb && self.face_valence(fa) != 3 => return false,
                _ => {} // Do nothing.
            }
        }
        // Check again if left and right are the same vertex.
        if let Some(h) = self.vertex_halfedge(v0) {
            if vertex_status[self.to_vertex(h).index() as usize].tagged()
                && vl == vr
                && htriangle
                && ohtriangle
            {
                return false;
            }
        }
        true
    }

    pub fn try_check_edge_collapse(&self, h: HH) -> Result<bool, Error> {
        let estatus = self.estatus.try_borrow()?;
        // Clone the property to side step compile time borrow checker. The
        // runtime borrow checker is still in use, so not a problem.
        let mut vstatus = self.vstatus.clone();
        let mut vstatus = vstatus.try_borrow_mut()?;
        Ok(self.check_edge_collapse(h, &estatus, &mut vstatus))
    }

    fn collapse_degenerate_triangle(
        &mut self,
        h: HH,
        hstatus: &mut [Status],
        estatus: &mut [Status],
        fstatus: &mut [Status],
    ) {
        let h1 = self.next_halfedge(h);
        let o = self.opposite_halfedge(h);
        let o1 = self.opposite_halfedge(h1);
        let v0 = self.to_vertex(h);
        let v1 = self.to_vertex(h1);
        let fh = self.halfedge_face(h);
        let fo = self.halfedge_face(o);
        // Ensure the loop represents a collapsed triangle. Because this is a
        // private function, all callers inside the implementation, so we can be
        // confident and assert to catch and weed out any bugs in debug builds.
        debug_assert_eq!(self.next_halfedge(h1), h);
        debug_assert_ne!(h1, o);
        // Rewire halfedge -> halfedge.
        self.link_halfedges(h1, self.next_halfedge(o));
        self.link_halfedges(self.prev_halfedge(o), h1);
        // Rewire halfedge -> face.
        self.halfedge_mut(h1).face = fo;
        // Rewire vertex -> halfedge.
        self.vertex_mut(v0).halfedge = Some(h1);
        self.adjust_outgoing_halfedge(v0);
        self.vertex_mut(v1).halfedge = Some(o1);
        self.adjust_outgoing_halfedge(v1);
        // Rewire face -> halfedge.
        if let Some(fo) = fo {
            if self.face_halfedge(fo) == o {
                self.face_mut(fo).halfedge = h1;
            }
        }
        // Delete stuff.
        if let Some(fh) = fh {
            fstatus[fh.index() as usize].set_deleted(true);
        }
        estatus[self.halfedge_edge(h).index() as usize].set_deleted(true);
        hstatus[h.index() as usize].set_deleted(true);
        hstatus[o.index() as usize].set_deleted(true);
    }

    pub fn collapse_edge(
        &mut self,
        h: HH,
        vstatus: &mut [Status],
        hstatus: &mut [Status],
        estatus: &mut [Status],
        fstatus: &mut [Status],
        cache: &mut TopolCache,
    ) {
        // Collect neighboring topology.
        let hn = self.next_halfedge(h);
        let hp = self.prev_halfedge(h);
        let o = self.opposite_halfedge(h);
        let on = self.next_halfedge(o);
        let op = self.prev_halfedge(o);
        let fh = self.halfedge_face(h);
        let fo = self.halfedge_face(o);
        let vh = self.to_vertex(h);
        let vo = self.to_vertex(o);
        // Setup cache.
        let hcache: &mut Vec<HH> = &mut cache.halfedges;
        // Rewire halfedge -> vertex
        hcache.clear();
        hcache.extend(iterator::vih_ccw_iter(self, vo));
        for ih in hcache.drain(..) {
            self.halfedge_mut(ih).vertex = vh;
        }
        // Rewire halfedge -> halfedge
        self.link_halfedges(hp, hn);
        self.link_halfedges(op, on);
        // Rewire face -> halfedge
        if let Some(fh) = fh {
            self.face_mut(fh).halfedge = hn;
        }
        if let Some(fo) = fo {
            self.face_mut(fo).halfedge = on;
        }
        // Rewire vertex -> halfedge
        if self.vertex_halfedge(vh) == Some(o) {
            self.vertex_mut(vh).halfedge = Some(hn);
        }
        self.adjust_outgoing_halfedge(vh);
        self.vertex_mut(vo).halfedge = None;
        // Delete stuff
        estatus[self.halfedge_edge(h).index() as usize].set_deleted(true);
        vstatus[vo.index() as usize].set_deleted(true);
        hstatus[h.index() as usize].set_deleted(true);
        hstatus[o.index() as usize].set_deleted(true);
        // If the loops that used to contain the halfedges that were collapsed
        // and deleted had a valance of 3, they are now degenerate. So we need
        // to collapse those loops.
        if self.next_halfedge(hn) == hp {
            self.collapse_degenerate_triangle(hn, hstatus, estatus, fstatus);
        }
        if self.next_halfedge(on) == op {
            self.collapse_degenerate_triangle(on, hstatus, estatus, fstatus);
        }
    }

    pub fn try_collapse_edge(&mut self, h: HH, cache: &mut TopolCache) -> Result<(), Error> {
        let mut vstatus = self.vstatus.clone();
        let mut vstatus = vstatus.try_borrow_mut()?;
        let mut hstatus = self.hstatus.clone();
        let mut hstatus = hstatus.try_borrow_mut()?;
        let mut estatus = self.estatus.clone();
        let mut estatus = estatus.try_borrow_mut()?;
        let mut fstatus = self.fstatus.clone();
        let mut fstatus = fstatus.try_borrow_mut()?;
        Ok(self.collapse_edge(
            h,
            &mut vstatus,
            &mut hstatus,
            &mut estatus,
            &mut fstatus,
            cache,
        ))
    }

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
    use crate::{
        element::Handle,
        iterator,
        topol::{
            test::{loop_mesh, quad_box},
            TopolCache,
        },
    };

    #[test]
    fn t_box_check_edge_collapse() {
        let qbox = quad_box();
        for h in qbox.halfedges() {
            assert!(qbox
                .try_check_edge_collapse(h)
                .expect("Cannot check halfedge collapse"));
        }
    }

    #[test]
    fn t_loop_mesh_check_edge_collapse() {
        let mesh = loop_mesh();
        assert_eq!(48, mesh.num_halfedges());
        // The edges spanning different boundary loops cannot be
        // collapsed. There are 8 of them, i.e. 16 halfedges. The remaining 32
        // can be collapsed.
        assert_eq!(
            (32, 16),
            mesh.halfedges().fold((0usize, 0usize), |(can, cannot), h| {
                if mesh
                    .try_check_edge_collapse(h)
                    .expect("Cannot check collapse")
                {
                    (can + 1, cannot)
                } else {
                    (can, 1 + cannot)
                }
            })
        );
    }

    #[test]
    fn t_box_edge_collapse() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        let h = qbox
            .find_halfedge(5.into(), 6.into())
            .expect("Cannot find halfedge");
        qbox.try_collapse_edge(h, &mut cache)
            .expect("Cannot collapse edges");
        assert_eq!(qbox.num_faces(), 6);
        assert_eq!(
            (2, 4),
            qbox.faces()
                .fold((0usize, 0usize), |(t, q), f| match qbox.face_valence(f) {
                    3 => (t + 1, q),
                    4 => (t, q + 1),
                    _ => (t, q),
                })
        );
        {
            let estatus = qbox
                .estatus
                .try_borrow()
                .expect("Cannot borrow the edge status property");
            assert_eq!(
                (1, 11),
                qbox.edges().fold((0usize, 0usize), |(del, ndel), e| {
                    if estatus[e.index() as usize].deleted() {
                        (del + 1, ndel)
                    } else {
                        (del, 1 + ndel)
                    }
                })
            );
        }
        {
            let hstatus = qbox
                .hstatus
                .try_borrow()
                .expect("Cannot borrow the halfedge status property");
            assert_eq!(
                (2, 22),
                qbox.halfedges().fold((0usize, 0usize), |(del, ndel), h| {
                    if hstatus[h.index() as usize].deleted() {
                        (del + 1, ndel)
                    } else {
                        (del, 1 + ndel)
                    }
                })
            );
        }
    }

    #[test]
    fn t_box_double_edge_collapse() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        // Collapse two opposite edges of a face, to produce a triangular prism.
        let h = qbox
            .find_halfedge(5.into(), 6.into())
            .expect("Cannot find halfedge");
        qbox.try_collapse_edge(h, &mut cache)
            .expect("Cannot collapse edges");
        let h = qbox
            .find_halfedge(4.into(), 7.into())
            .expect("Cannot find halfedge");
        qbox.try_collapse_edge(h, &mut cache)
            .expect("Cannot collapse edge");
        {
            let estatus = qbox
                .estatus
                .try_borrow()
                .expect("Cannot borrow the edge status property");
            assert_eq!(
                (3, 9),
                qbox.edges().fold((0usize, 0usize), |(del, ndel), e| {
                    if estatus[e.index() as usize].deleted() {
                        (del + 1, ndel)
                    } else {
                        (del, 1 + ndel)
                    }
                })
            );
        }
        {
            let hstatus = qbox
                .hstatus
                .try_borrow()
                .expect("Cannot borrow the halfedge status property");
            assert_eq!(
                (6, 18),
                qbox.halfedges().fold((0usize, 0usize), |(del, ndel), h| {
                    if hstatus[h.index() as usize].deleted() {
                        (del + 1, ndel)
                    } else {
                        (del, 1 + ndel)
                    }
                })
            );
        }
        {
            let fstatus = qbox
                .fstatus
                .try_borrow()
                .expect("Cannot borrow the halfedge status property");
            assert_eq!(
                (1, 5),
                qbox.faces().fold((0usize, 0usize), |(del, ndel), h| {
                    if fstatus[h.index() as usize].deleted() {
                        (del + 1, ndel)
                    } else {
                        (del, 1 + ndel)
                    }
                })
            );
            assert_eq!(
                (2, 3),
                qbox.faces().fold((0usize, 0usize), |(t, q), f| {
                    if fstatus[f.index() as usize].deleted() {
                        (t, q)
                    } else {
                        match qbox.face_valence(f) {
                            3 => (t + 1, q),
                            4 => (t, q + 1),
                            _ => (t, q),
                        }
                    }
                })
            );
        }
        qbox.garbage_collection(&mut cache)
            .expect("Garbage collection failed");
        assert_eq!(6, qbox.num_vertices());
        assert_eq!(18, qbox.num_halfedges());
        assert_eq!(9, qbox.num_edges());
        assert_eq!(5, qbox.num_faces());
    }

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
