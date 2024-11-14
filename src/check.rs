use crate::{element::Handle, iterator, topol::Topology, Error, Status};

impl Topology {
    fn check_vertices(
        &self,
        vstatus: &[Status],
        hstatus: &[Status],
        hvisited: &mut [bool],
    ) -> Result<(), Error> {
        hvisited.fill(false);
        for v in self
            .vertices()
            .filter(|v| !vstatus[v.index() as usize].deleted())
        {
            if let Some(h) = self.vertex_halfedge(v) {
                // Invalid index, or deleted halfedge.
                if !self.is_valid_halfedge(h) {
                    return Err(Error::InvalidHalfedge(h));
                }
                if hstatus[h.index() as usize].deleted() {
                    return Err(Error::DeletedHalfedge(h));
                }
                // The outgoing halfedge must be a boundary halfedge, or none of the
                // halfedges are boundary.
                if !self.is_boundary_halfedge(h)
                    && iterator::voh_ccw_iter(self, v).any(|h| self.is_boundary_halfedge(h))
                {
                    return Err(Error::OutgoingHalfedgeNotBoundary(v));
                }
                // Outgoing halfedge must point back to this vertex.
                if self.tail_vertex(h) != v {
                    return Err(Error::InvalidOutgoingHalfedges(v));
                }
            }
            // Check ccw iterator.
            for h in iterator::voh_ccw_iter(self, v) {
                if std::mem::replace(&mut hvisited[h.index() as usize], true) {
                    return Err(Error::InvalidOutgoingHalfedges(v));
                }
            }
            // Check cw iterator.
            for h in iterator::voh_cw_iter(self, v) {
                if !std::mem::replace(&mut hvisited[h.index() as usize], false) {
                    return Err(Error::InvalidOutgoingHalfedges(v));
                }
            }
        }
        Ok(())
    }

    fn check_edges(
        &self,
        vstatus: &[Status],
        hstatus: &[Status],
        estatus: &[Status],
        fstatus: &[Status],
        hflags: &mut [bool],
    ) -> Result<(), Error> {
        for h in self
            .halfedges()
            .filter(|h| !hstatus[h.index() as usize].deleted())
        {
            // Check if degenerate.
            if self.tail_vertex(h) == self.head_vertex(h) {
                return Err(Error::DegenerateHalfedge(h));
            }
            let hedge = self.halfedge(h);
            // Check for deleted.
            if hstatus[hedge.prev.index() as usize].deleted() {
                return Err(Error::DeletedHalfedge(hedge.prev));
            }
            if hstatus[hedge.next.index() as usize].deleted() {
                return Err(Error::DeletedHalfedge(hedge.next));
            }
            if vstatus[hedge.vertex.index() as usize].deleted() {
                return Err(Error::DeletedVertex(hedge.vertex));
            }
            if let Some(f) = hedge.face {
                if fstatus[f.index() as usize].deleted() {
                    return Err(Error::DeletedFace(f));
                }
            }
            let e = self.halfedge_edge(h);
            if estatus[e.index() as usize].deleted() {
                return Err(Error::DeletedEdge(e));
            }
            // Check connctivity.
            if self.next_halfedge(hedge.prev) != h
                || self.prev_halfedge(hedge.next) != h
                || self.head_vertex(h) != self.tail_vertex(hedge.next)
                || self.tail_vertex(h) != self.head_vertex(hedge.prev)
            {
                return Err(Error::InvalidHalfedgeLink(h));
            }
        }
        // Check all loops.
        for h in self
            .halfedges()
            .filter(|h| !hstatus[h.index() as usize].deleted())
        {
            if hflags[h.index() as usize] {
                continue;
            }
            let f = self.halfedge_face(h);
            for h in iterator::loop_ccw_iter(self, h) {
                if std::mem::replace(&mut hflags[h.index() as usize], true) {
                    return Err(Error::InvalidLoopTopology(h));
                }
                if self.halfedge_face(h) != f {
                    return Err(Error::InconsistentFaceInLoop(h));
                }
            }
        }
        for h in self
            .halfedges()
            .filter(|h| !hstatus[h.index() as usize].deleted())
        {
            if !hflags[h.index() as usize] {
                continue;
            }
            let f = self.halfedge_face(h);
            for h in iterator::loop_ccw_iter(self, h) {
                if !std::mem::replace(&mut hflags[h.index() as usize], false) {
                    return Err(Error::InvalidLoopTopology(h));
                }
                if self.halfedge_face(h) != f {
                    return Err(Error::InconsistentFaceInLoop(h));
                }
            }
        }
        assert!(hflags.iter().all(|f| !f));
        Ok(())
    }

    fn check_faces(&self, hstatus: &[Status], fstatus: &[Status]) -> Result<(), Error> {
        for f in self
            .faces()
            .filter(|f| !fstatus[f.index() as usize].deleted())
        {
            let h = self.face_halfedge(f);
            if hstatus[h.index() as usize].deleted() {
                return Err(Error::DeletedHalfedge(h));
            }
            if self.halfedge_face(h) != Some(f) {
                return Err(Error::InvalidFaceHalfedgeLink(f, h));
            }
        }
        Ok(())
    }

    pub fn check(&self) -> Result<(), Error> {
        let vstatus = self.vstatus.try_borrow()?;
        let vstatus: &[Status] = &vstatus;
        let hstatus = self.hstatus.try_borrow()?;
        let hstatus: &[Status] = &hstatus;
        let fstatus = self.fstatus.try_borrow()?;
        let fstatus: &[Status] = &fstatus;
        let estatus = self.estatus.try_borrow()?;
        let estatus: &[Status] = &estatus;
        // To keep track of visited halfedges.
        let mut hvisited = vec![false; self.num_halfedges()].into_boxed_slice();
        self.check_vertices(vstatus, hstatus, &mut hvisited)?;
        self.check_edges(vstatus, hstatus, estatus, fstatus, &mut hvisited)?;
        self.check_faces(hstatus, fstatus)?;
        Ok(())
    }
}
