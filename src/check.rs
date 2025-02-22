use crate::{
    EPropBuf, Error, FPropBuf, HPropBuf, Status, VPropBuf,
    element::Handle,
    iterator::HasIterators,
    topol::{HasTopology, Topology},
};

fn check_vertices(
    mesh: &Topology,
    vstatus: &VPropBuf<Status>,
    hstatus: &HPropBuf<Status>,
    hvisited: &mut [bool],
) -> Result<(), Error> {
    hvisited.fill(false);
    for v in mesh.vertices().filter(|v| !vstatus[*v].deleted()) {
        if let Some(h) = v.halfedge(mesh) {
            // Invalid index, or deleted halfedge.
            if !h.is_valid(mesh) {
                return Err(Error::InvalidHalfedge(h));
            }
            if hstatus[h].deleted() {
                return Err(Error::DeletedHalfedge(h));
            }
            // Outgoing halfedge must point back to this vertex.
            if h.tail(mesh) != v {
                return Err(Error::InvalidOutgoingHalfedges(v));
            }
        }
        // Check ccw iterator.
        for h in mesh.voh_ccw_iter(v) {
            if std::mem::replace(&mut hvisited[h.index() as usize], true) {
                return Err(Error::InvalidOutgoingHalfedges(v));
            }
        }
        // Check cw iterator.
        for h in mesh.voh_cw_iter(v) {
            if !std::mem::replace(&mut hvisited[h.index() as usize], false) {
                return Err(Error::InvalidOutgoingHalfedges(v));
            }
        }
        if let Some(h) = v.halfedge(mesh) {
            // The outgoing halfedge must be a boundary halfedge, or none of the
            // halfedges are boundary.
            if !h.is_boundary(mesh) && mesh.voh_ccw_iter(v).any(|h| h.is_boundary(mesh)) {
                return Err(Error::OutgoingHalfedgeNotBoundary(v));
            }
        }
    }
    Ok(())
}

fn check_edges(
    mesh: &Topology,
    vstatus: &VPropBuf<Status>,
    hstatus: &HPropBuf<Status>,
    estatus: &EPropBuf<Status>,
    fstatus: &FPropBuf<Status>,
    hflags: &mut [bool],
) -> Result<(), Error> {
    for h in mesh.halfedges().filter(|h| !hstatus[*h].deleted()) {
        // Check if degenerate.
        if h.tail(mesh) == h.head(mesh) {
            return Err(Error::DegenerateHalfedge(h));
        }
        let hedge = mesh.halfedge(h);
        // Check for deleted.
        if hstatus[hedge.prev].deleted() {
            return Err(Error::DeletedHalfedge(hedge.prev));
        }
        if hstatus[hedge.next].deleted() {
            return Err(Error::DeletedHalfedge(hedge.next));
        }
        if vstatus[hedge.vertex].deleted() {
            return Err(Error::DeletedVertex(hedge.vertex));
        }
        if let Some(f) = hedge.face {
            if fstatus[f].deleted() {
                return Err(Error::DeletedFace(f));
            }
        }
        let e = h.edge();
        if estatus[e].deleted() {
            return Err(Error::DeletedEdge(e));
        }
        // Check connctivity.
        let head = h.head(mesh);
        let tail = h.tail(mesh);
        if hedge.prev.next(mesh) != h
            || hedge.next.prev(mesh) != h
            || head != hedge.next.tail(mesh)
            || tail != hedge.prev.head(mesh)
        {
            return Err(Error::InvalidHalfedgeLink(h));
        }
        // Halfedge must be found in the circulators around head and tail.
        if !mesh.voh_ccw_iter(tail).any(|hh| hh == h) || !mesh.vih_ccw_iter(head).any(|hh| hh == h)
        {
            return Err(Error::InvalidHalfedgeVertexLink(h));
        }
    }
    // Check all loops.
    for h in mesh.halfedges().filter(|h| !hstatus[*h].deleted()) {
        if hflags[h.index() as usize] {
            continue;
        }
        let f = h.face(mesh);
        for h in mesh.loop_ccw_iter(h) {
            if std::mem::replace(&mut hflags[h.index() as usize], true) {
                return Err(Error::InvalidLoopTopology(h));
            }
            if h.face(mesh) != f {
                return Err(Error::InconsistentFaceInLoop(h));
            }
        }
    }
    for h in mesh.halfedges().filter(|h| !hstatus[*h].deleted()) {
        if !hflags[h.index() as usize] {
            continue;
        }
        let f = h.face(mesh);
        for h in mesh.loop_ccw_iter(h) {
            if !std::mem::replace(&mut hflags[h.index() as usize], false) {
                return Err(Error::InvalidLoopTopology(h));
            }
            if h.face(mesh) != f {
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
        let h = f.halfedge(mesh);
        if hstatus[h.index() as usize].deleted() {
            return Err(Error::DeletedHalfedge(h));
        }
        if h.face(mesh) != Some(f) {
            return Err(Error::InvalidFaceHalfedgeLink(f, h));
        }
    }
    Ok(())
}

impl Topology {
    pub fn check(&self) -> Result<(), Error> {
        let vstatus = self.vstatus.try_borrow()?;
        let hstatus = self.hstatus.try_borrow()?;
        let fstatus = self.fstatus.try_borrow()?;
        let estatus = self.estatus.try_borrow()?;
        // To keep track of visited halfedges.
        let mut hvisited = vec![false; self.num_halfedges()].into_boxed_slice();
        check_vertices(self, &vstatus, &hstatus, &mut hvisited)?;
        check_edges(self, &vstatus, &hstatus, &estatus, &fstatus, &mut hvisited)?;
        check_faces(self, &hstatus, &fstatus)?;
        Ok(())
    }
}
