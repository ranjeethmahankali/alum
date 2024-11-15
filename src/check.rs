use crate::{element::Handle, iterator, topol::Topology, Adaptor, Error, PolyMeshT, Status};

fn check_vertices(
    mesh: &Topology,
    vstatus: &[Status],
    hstatus: &[Status],
    hvisited: &mut [bool],
) -> Result<(), Error> {
    hvisited.fill(false);
    for v in mesh
        .vertices()
        .filter(|v| !vstatus[v.index() as usize].deleted())
    {
        if let Some(h) = mesh.vertex_halfedge(v) {
            // Invalid index, or deleted halfedge.
            if !h.is_valid(mesh) {
                return Err(Error::InvalidHalfedge(h));
            }
            if hstatus[h.index() as usize].deleted() {
                return Err(Error::DeletedHalfedge(h));
            }
            // The outgoing halfedge must be a boundary halfedge, or none of the
            // halfedges are boundary.
            if !mesh.is_boundary_halfedge(h)
                && iterator::voh_ccw_iter(mesh, v).any(|h| mesh.is_boundary_halfedge(h))
            {
                return Err(Error::OutgoingHalfedgeNotBoundary(v));
            }
            // Outgoing halfedge must point back to this vertex.
            if mesh.tail_vertex(h) != v {
                return Err(Error::InvalidOutgoingHalfedges(v));
            }
        }
        // Check ccw iterator.
        for h in iterator::voh_ccw_iter(mesh, v) {
            if std::mem::replace(&mut hvisited[h.index() as usize], true) {
                return Err(Error::InvalidOutgoingHalfedges(v));
            }
        }
        // Check cw iterator.
        for h in iterator::voh_cw_iter(mesh, v) {
            if !std::mem::replace(&mut hvisited[h.index() as usize], false) {
                return Err(Error::InvalidOutgoingHalfedges(v));
            }
        }
    }
    Ok(())
}

fn check_edges(
    mesh: &Topology,
    vstatus: &[Status],
    hstatus: &[Status],
    estatus: &[Status],
    fstatus: &[Status],
    hflags: &mut [bool],
) -> Result<(), Error> {
    for h in mesh
        .halfedges()
        .filter(|h| !hstatus[h.index() as usize].deleted())
    {
        // Check if degenerate.
        if mesh.tail_vertex(h) == mesh.head_vertex(h) {
            return Err(Error::DegenerateHalfedge(h));
        }
        let hedge = mesh.halfedge(h);
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
        let e = h.edge();
        if estatus[e.index() as usize].deleted() {
            return Err(Error::DeletedEdge(e));
        }
        // Check connctivity.
        let head = mesh.head_vertex(h);
        let tail = mesh.tail_vertex(h);
        if mesh.next_halfedge(hedge.prev) != h
            || mesh.prev_halfedge(hedge.next) != h
            || head != mesh.tail_vertex(hedge.next)
            || tail != mesh.head_vertex(hedge.prev)
        {
            return Err(Error::InvalidHalfedgeLink(h));
        }
        // Halfedge must be found in the circulators around head and tail.
        if !iterator::voh_ccw_iter(mesh, tail).any(|hh| hh == h)
            || !iterator::vih_ccw_iter(mesh, head).any(|hh| hh == h)
        {
            return Err(Error::InvalidHalfedgeVertexLink(h));
        }
    }
    // Check all loops.
    for h in mesh
        .halfedges()
        .filter(|h| !hstatus[h.index() as usize].deleted())
    {
        if hflags[h.index() as usize] {
            continue;
        }
        let f = mesh.halfedge_face(h);
        for h in iterator::loop_ccw_iter(mesh, h) {
            if std::mem::replace(&mut hflags[h.index() as usize], true) {
                return Err(Error::InvalidLoopTopology(h));
            }
            if mesh.halfedge_face(h) != f {
                return Err(Error::InconsistentFaceInLoop(h));
            }
        }
    }
    for h in mesh
        .halfedges()
        .filter(|h| !hstatus[h.index() as usize].deleted())
    {
        if !hflags[h.index() as usize] {
            continue;
        }
        let f = mesh.halfedge_face(h);
        for h in iterator::loop_ccw_iter(mesh, h) {
            if !std::mem::replace(&mut hflags[h.index() as usize], false) {
                return Err(Error::InvalidLoopTopology(h));
            }
            if mesh.halfedge_face(h) != f {
                return Err(Error::InconsistentFaceInLoop(h));
            }
        }
    }
    assert!(hflags.iter().all(|f| !f));
    Ok(())
}

fn check_faces(mesh: &Topology, hstatus: &[Status], fstatus: &[Status]) -> Result<(), Error> {
    for f in mesh
        .faces()
        .filter(|f| !fstatus[f.index() as usize].deleted())
    {
        let h = mesh.face_halfedge(f);
        if hstatus[h.index() as usize].deleted() {
            return Err(Error::DeletedHalfedge(h));
        }
        if mesh.halfedge_face(h) != Some(f) {
            return Err(Error::InvalidFaceHalfedgeLink(f, h));
        }
    }
    Ok(())
}

impl Topology {
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
        check_vertices(self, vstatus, hstatus, &mut hvisited)?;
        check_edges(self, vstatus, hstatus, estatus, fstatus, &mut hvisited)?;
        check_faces(self, hstatus, fstatus)?;
        Ok(())
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    /// Check the topology of the mesh.
    ///
    /// This function will return an error if any errors are found in the topolgy.
    pub fn check_topology(&self) -> Result<(), Error> {
        self.topol.check()
    }
}
