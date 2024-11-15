use crate::{
    element::{Handle, EH, FH, HH, VH},
    error::Error,
    iterator::HasIterators,
    mesh::{Adaptor, PolyMeshT},
    status::Status,
    topol::{TopolCache, Topology},
    HasTopology,
};

impl Topology {
    /// Check if it is safe to collapse an edge.
    ///
    /// This function uses the edge-status and vertex-status(mutable) properties
    /// to keep track of the neighborhood. Use this with borrowed properties as
    /// arguments if you're calling this function in a hot loop, to avoid
    /// repeated borrows.
    pub fn check_edge_collapse(
        &self,
        h: HH,
        edge_status: &[Status],
        vertex_status: &mut [Status],
    ) -> bool {
        // Check if already deleted.
        if edge_status[h.edge().index() as usize].deleted() {
            return false;
        }
        let oh = h.opposite();
        let v0 = oh.head(self);
        let v1 = h.head(self);
        if vertex_status[v0.index() as usize].deleted()
            || vertex_status[v1.index() as usize].deleted()
        {
            return false;
        }
        let htriangle = match h.face(self) {
            Some(f) => f.valence(self) == 3,
            None => false,
        };
        let ohtriangle = match oh.face(self) {
            Some(f) => f.valence(self) == 3,
            None => false,
        };
        // Check if the faces are triangles and the vertices opposite the edge
        // on those triangles are actually the same vertex.
        let vl = if htriangle {
            let h1 = h.next(self);
            let h2 = h1.next(self);
            if h1.opposite().is_boundary(self) && h2.opposite().is_boundary(self) {
                return false;
            }
            Some(h1.head(self))
        } else {
            None
        };
        let vr = if ohtriangle {
            let h1 = oh.next(self);
            let h2 = h1.next(self);
            if h1.opposite().is_boundary(self) && h2.opposite().is_boundary(self) {
                return false;
            }
            Some(h1.head(self))
        } else {
            None
        };
        if let (Some(vl), Some(vr)) = (vl, vr) {
            if vl == vr {
                return false;
            }
        }
        // Check if we're collapsing across two different boundaries.
        if v0.is_boundary(self)
            && v1.is_boundary(self)
            && !h.is_boundary(self)
            && !oh.is_boundary(self)
        {
            return false;
        }
        // Check the 'Link condition' Edelsbrunner [2006]. The intersection of
        // the one rings of from and to vertices must be the left and right
        // vertices, and only if the corresponding faces are triangles.
        let vl = h.next(self).head(self);
        let vr = oh.next(self).head(self);
        for v in self.vv_ccw_iter(v0) {
            vertex_status[v.index() as usize].set_tagged(false);
        }
        for v in self.vv_ccw_iter(v1) {
            vertex_status[v.index() as usize].set_tagged(true);
        }
        for v in self.vv_ccw_iter(v0) {
            if vertex_status[v.index() as usize].tagged()
                && !(v == vl && htriangle)
                && !(v == vr && ohtriangle)
            {
                return false;
            }
        }
        // Check for folded faces that might degenerate.
        if htriangle {
            let h1 = h.next(self).opposite();
            let h2 = h.prev(self).opposite();
            match (h1.face(self), h2.face(self)) {
                (None, None) => return false, // This is redundant but just in case.
                (Some(fa), Some(fb)) if fa == fb && fa.valence(self) != 3 => return false,
                _ => {} // Do nothing.
            }
        }
        if ohtriangle {
            let h1 = oh.next(self).opposite();
            let h2 = oh.prev(self).opposite();
            match (h1.face(self), h2.face(self)) {
                (None, None) => return false, // This is redundant but just in case.
                (Some(fa), Some(fb)) if fa == fb && fa.valence(self) != 3 => return false,
                _ => {} // Do nothing.
            }
        }
        // Check again if left and right are the same vertex.
        if let Some(h) = v0.halfedge(self) {
            if vertex_status[h.head(self).index() as usize].tagged()
                && vl == vr
                && htriangle
                && ohtriangle
            {
                return false;
            }
        }
        true
    }

    /// This is the same as [`Self::check_edge_collapse`], except it will attempt to
    /// borrow the necessary properties and may return an error if it ccannot
    /// borrow the required properties.
    pub fn try_check_edge_collapse(&self, h: HH) -> Result<bool, Error> {
        let estatus = self.estatus.try_borrow()?;
        // Clone the property to side step compile time borrow checker. The
        // runtime borrow checker is still in use, so not a problem.
        let mut vstatus = self.vstatus.clone();
        let mut vstatus = vstatus.try_borrow_mut()?;
        Ok(self.check_edge_collapse(h, &estatus, &mut vstatus))
    }

    /// Sometimes after collapsing an edge, if the neighboring faces are
    /// triangles, we end up with degenerate loops / faces. This cleans up such
    /// loops.
    fn collapse_degenerate_triangle(
        &mut self,
        h: HH,
        hstatus: &mut [Status],
        estatus: &mut [Status],
        fstatus: &mut [Status],
    ) {
        let h1 = h.next(self);
        let o = h.opposite();
        let o1 = h1.opposite();
        let v0 = h.head(self);
        let v1 = h1.head(self);
        let fh = h.face(self);
        let fo = o.face(self);
        // Ensure the loop represents a collapsed triangle. Because this is a
        // private function, all callers inside the implementation, so we can be
        // confident and assert to catch and weed out any bugs in debug builds.
        debug_assert_eq!(h1.next(self), h);
        debug_assert_ne!(h1, o);
        // Rewire halfedge -> halfedge.
        self.link_halfedges(h1, o.next(self));
        self.link_halfedges(o.prev(self), h1);
        // Rewire halfedge -> face.
        self.halfedge_mut(h1).face = fo;
        // Rewire vertex -> halfedge.
        self.vertex_mut(v0).halfedge = Some(h1);
        self.adjust_outgoing_halfedge(v0);
        self.vertex_mut(v1).halfedge = Some(o1);
        self.adjust_outgoing_halfedge(v1);
        // Rewire face -> halfedge.
        if let Some(fo) = fo {
            if fo.halfedge(self) == o {
                self.face_mut(fo).halfedge = h1;
            }
        }
        // Delete stuff.
        if let Some(fh) = fh {
            fstatus[fh.index() as usize].set_deleted(true);
        }
        estatus[h.edge().index() as usize].set_deleted(true);
        hstatus[h.index() as usize].set_deleted(true);
        hstatus[o.index() as usize].set_deleted(true);
    }

    /// Collapse an edge.
    ///
    /// The vertex at the start of the halfedge will be deleted, and all the
    /// topology associated with that vertex will be rewired to the vertex that
    /// the halfedge is pointing towards. This function requires borrowed
    /// status properties of vertices, halfedges, edges and faces. This is
    /// useful when collapsing edges in a hot loop, to avoid repeated borrows of
    /// properties.
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
        let hn = h.next(self);
        let hp = h.prev(self);
        let o = h.opposite();
        let on = o.next(self);
        let op = o.prev(self);
        let fh = h.face(self);
        let fo = o.face(self);
        let vh = h.head(self);
        let vo = o.head(self);
        // Setup cache.
        let hcache: &mut Vec<HH> = &mut cache.halfedges;
        // Rewire halfedge -> vertex
        hcache.clear();
        hcache.extend(self.vih_ccw_iter(vo));
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
        if vh.halfedge(self) == Some(o) {
            self.vertex_mut(vh).halfedge = Some(hn);
        }
        self.adjust_outgoing_halfedge(vh);
        self.vertex_mut(vo).halfedge = None;
        // Delete stuff
        estatus[h.edge().index() as usize].set_deleted(true);
        vstatus[vo.index() as usize].set_deleted(true);
        hstatus[h.index() as usize].set_deleted(true);
        hstatus[o.index() as usize].set_deleted(true);
        // If the loops that used to contain the halfedges that were collapsed
        // and deleted had a valance of 3, they are now degenerate. So we need
        // to collapse those loops.
        if hn.next(self) == hp {
            self.collapse_degenerate_triangle(hn, hstatus, estatus, fstatus);
        }
        if on.next(self) == op {
            self.collapse_degenerate_triangle(on, hstatus, estatus, fstatus);
        }
    }

    /// This is the same as [`Self::collapse_edge`] except it will attempt to
    /// borrow all the required properties, and returns an error if borrowing
    /// fails.
    pub fn try_collapse_edge(&mut self, h: HH, cache: &mut TopolCache) -> Result<(), Error> {
        let mut vstatus = self.vstatus.clone();
        let mut vstatus = vstatus.try_borrow_mut()?;
        let mut hstatus = self.hstatus.clone();
        let mut hstatus = hstatus.try_borrow_mut()?;
        let mut estatus = self.estatus.clone();
        let mut estatus = estatus.try_borrow_mut()?;
        let mut fstatus = self.fstatus.clone();
        let mut fstatus = fstatus.try_borrow_mut()?;
        self.collapse_edge(
            h,
            &mut vstatus,
            &mut hstatus,
            &mut estatus,
            &mut fstatus,
            cache,
        );
        Ok(())
    }

    pub fn triangulate_face(&mut self, f: FH) -> Result<(), Error> {
        let mut base = f.halfedge(self);
        let vstart = base.tail(self);
        let prev = base.prev(self);
        let mut next = base.next(self);
        while next.next(self).head(self) != vstart {
            let next2 = next.next(self);
            let fnew = self.new_face(base)?;
            let enew = self.new_edge(vstart, next.head(self), prev, next2, next, base)?;
            let hnew = enew.halfedge(false);
            let ohnew = enew.halfedge(true);
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
        self.link_halfedges(next.next(self), base);
        self.halfedge_mut(base).face = Some(f);
        Ok(())
    }

    pub fn triangulate(&mut self) -> Result<(), Error> {
        for f in self.faces() {
            self.triangulate_face(f)?;
        }
        Ok(())
    }

    pub fn split_edge(&mut self, e: EH, v: VH, copy_props: bool) -> Result<EH, Error> {
        let (h0, h1) = e.halfedges();
        let vfrom = h0.tail(self);
        let (ph0, nh1) = (h0.prev(self), h1.next(self));
        let (f0, f1) = (h0.face(self), h1.face(self));
        // Create a new edge and rewire topology.
        let enew = self.new_edge(vfrom, v, ph0, h0, h1, nh1)?;
        let hnew = enew.halfedge(false);
        let ohnew = enew.halfedge(true);
        // Rewire halfedge -> vertex.
        self.halfedge_mut(h1).vertex = v;
        // Rewire halfedge -> halfedge.
        self.link_halfedges(hnew, h0);
        self.link_halfedges(h1, ohnew);
        self.link_halfedges(ph0, hnew);
        self.link_halfedges(ohnew, nh1);
        // Rewire halfedge -> face.
        self.halfedge_mut(hnew).face = f0;
        self.halfedge_mut(ohnew).face = f1;
        // Rewire vertex -> halfedge.
        self.vertex_mut(v).halfedge = Some(h0);
        self.adjust_outgoing_halfedge(v);
        if vfrom.halfedge(self) == Some(h0) {
            self.vertex_mut(vfrom).halfedge = Some(hnew);
            self.adjust_outgoing_halfedge(vfrom);
        }
        if copy_props {
            self.eprops.copy(e, enew)?;
            self.hprops.copy_many(&[h0, h1], &[hnew, ohnew])?;
        }
        Ok(enew)
    }

    pub fn swap_edge_ccw(&mut self, e: EH) -> bool {
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let (f, of) = match (h.face(self), oh.face(self)) {
            (Some(f), Some(of)) => (f, of),
            _ => return false, // Cannot swap boundary edge.
        };
        let hn = h.next(self);
        let on = oh.next(self);
        let v0 = h.tail(self);
        let v1 = h.head(self);
        // Check for degeneracy.
        if f == of || hn == oh || hn.head(self) == v0 || on == h || on.head(self) == v1 {
            return false;
        }
        let hp = h.prev(self);
        let op = oh.prev(self);
        let hnn = hn.next(self);
        let onn = on.next(self);
        let hnv = hn.head(self);
        let onv = on.head(self);
        // Rewire vertex -> halfedge.
        if v0.halfedge(self) == Some(h) {
            self.vertex_mut(v0).halfedge = Some(on);
        }
        if v1.halfedge(self) == Some(oh) {
            self.vertex_mut(v1).halfedge = Some(hn);
        }
        // Rewire halfedge -> vertex.
        self.halfedge_mut(h).vertex = hnv;
        self.halfedge_mut(oh).vertex = onv;
        // Rewire halfedge -> halfedge.
        self.link_halfedges(oh, onn);
        self.link_halfedges(op, hn);
        self.link_halfedges(hn, oh);
        self.link_halfedges(h, hnn);
        self.link_halfedges(hp, on);
        self.link_halfedges(on, h);
        // Rewire halfedge -> face.
        self.halfedge_mut(hn).face = Some(of);
        self.halfedge_mut(on).face = Some(f);
        // Rewire face -> halfedge.
        if f.halfedge(self) == hn {
            self.face_mut(f).halfedge = h;
        }
        if of.halfedge(self) == on {
            self.face_mut(of).halfedge = oh;
        }
        true
    }

    pub fn swap_edge_cw(&mut self, e: EH) -> bool {
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let (f, of) = match (h.face(self), oh.face(self)) {
            (Some(f), Some(of)) => (f, of),
            _ => return false, // Cannot swap boundary edge.
        };
        let hn = h.next(self);
        let on = oh.next(self);
        let v0 = h.tail(self);
        let v1 = h.head(self);
        // Check for degeneracy.
        if f == of || hn == oh || hn.head(self) == v0 || on == h || on.head(self) == v1 {
            return false;
        }
        let hp = h.prev(self);
        let op = oh.prev(self);
        let hpp = hp.prev(self);
        let opp = op.prev(self);
        let hpv = hp.tail(self);
        let opv = op.tail(self);
        // Rewire vertex -> halfedge.
        if v0.halfedge(self) == Some(h) {
            self.vertex_mut(v0).halfedge = Some(on);
        }
        if v1.halfedge(self) == Some(oh) {
            self.vertex_mut(v1).halfedge = Some(hn);
        }
        // Rewire halfedge -> vertex.
        self.halfedge_mut(h).vertex = opv;
        self.halfedge_mut(oh).vertex = hpv;
        // Rewire halfedge -> halfedge.
        self.link_halfedges(h, op);
        self.link_halfedges(op, hn);
        self.link_halfedges(hpp, h);
        self.link_halfedges(oh, hp);
        self.link_halfedges(hp, on);
        self.link_halfedges(opp, oh);
        // Rewire halfedge -> face.
        self.halfedge_mut(op).face = Some(f);
        self.halfedge_mut(hp).face = Some(of);
        // Rewire face -> halfedge.
        if f.halfedge(self) == hp {
            self.face_mut(f).halfedge = h;
        }
        if of.halfedge(self) == op {
            self.face_mut(of).halfedge = oh;
        }
        true
    }

    /// An edge is a unique link if it is the only edge connecting the two faces
    /// incident on it.
    ///
    /// The boundary is treated as one face. So the boundary edges can only be
    /// simple links if
    fn edge_is_unique_link(&self, e: EH) -> bool {
        let h = e.halfedge(false);
        let fo = h.opposite().face(self);
        self.loop_ccw_iter(h)
            .skip(1)
            .all(|h| h.opposite().face(self) != fo)
    }

    pub fn remove_edge(&mut self, e: EH) -> Result<FH, Error> {
        //    <--------- <----------v0<---------- <----------
        //   |                n0    ^|     p1                ^
        //   |                      ||                       |
        //   |                      ||                       |
        //   |         f0         h0||h1         f1          |
        //   |                      ||                       |
        //   |                      ||                       |
        //   v                p0    |v     n1                |
        //    ---------> ---------->v1----------> ---------->
        let mut estatus = self.estatus.clone();
        let mut estatus = estatus.try_borrow_mut()?;
        let mut fstatus = self.fstatus.clone();
        let mut fstatus = fstatus.try_borrow_mut()?;
        let mut hstatus = self.hstatus.clone();
        let mut hstatus = hstatus.try_borrow_mut()?;
        if estatus[e.index() as usize].deleted() {
            return Err(Error::DeletedEdge(e));
        }
        if !self.edge_is_unique_link(e) {
            return Err(Error::EdgeIsNotAUniqueLink(e));
        }
        let (h0, h1) = e.halfedges();
        let (f0, f1) = match (h0.face(self), h1.face(self)) {
            (Some(f0), Some(f1)) => (f0, f1),
            _ => return Err(Error::CannotRemoveBoundaryEdge(e)),
        };
        let (p0, p1) = (h0.prev(self), h1.prev(self));
        let (n0, n1) = (h0.next(self), h1.next(self));
        let (v0, v1) = (h0.head(self), h1.head(self));
        // Rewire vertex -> halfedge.
        if v0.halfedge(self) == Some(h1) {
            self.vertex_mut(v0).halfedge = Some(n0);
        }
        if v1.halfedge(self) == Some(h0) {
            self.vertex_mut(v1).halfedge = Some(n1);
        }
        // Rewire halfedge -> halfedge.
        self.link_halfedges(p0, n1);
        self.link_halfedges(p1, n0);
        // Rewire face -> halfedge. Keep f0 and delete f1.
        if f0.halfedge(self) == h0 {
            self.face_mut(f0).halfedge = p0;
        }
        // Rewire halfedge -> face for the loop of f1.
        for (mesh, h) in self.loop_ccw_iter_mut(n1).take_while(|(_, h)| *h != n0) {
            mesh.halfedge_mut(h).face = Some(f0);
        }
        estatus[e.index() as usize].set_deleted(true);
        hstatus[h0.index() as usize].set_deleted(true);
        hstatus[h1.index() as usize].set_deleted(true);
        fstatus[f1.index() as usize].set_deleted(true);
        Ok(f0)
    }

    pub fn insert_edge(&mut self, prev: HH, next: HH) -> Result<EH, Error> {
        //    <--------- <----------v1<---------- <----------
        //   |                next  ^|     n0                ^
        //   |                      ||                       |
        //   |                      ||                       |
        //   |                     h||oh                     |
        //   |                      ||                       |
        //   |                      ||                       |
        //   v                prev  |v     p1                |
        //    ---------> ---------->v0----------> ---------->
        let (p1, n0) = (prev.next(self), next.prev(self));
        if p1 == next || prev == next {
            return Err(Error::CannotInsertEdge(prev, next));
        }
        let f = prev.face(self);
        if f != next.face(self) {
            return Err(Error::HalfedgesNotInTheSameLoop(prev, next));
        }
        if f.is_none() {
            // Check if the halfedges are part of the same boundary loop. March
            // simultaenously starting from both prev and next halfedges, and
            // see if you arrive at the other halfedge. Marching from both
            // ensures we'll detect the loop as soon as possible.
            if !self
                .loop_ccw_iter(prev)
                .zip(self.loop_cw_iter(next))
                .any(|(n, p)| n == prev || p == next)
            {
                return Err(Error::HalfedgesNotInTheSameLoop(prev, next));
            }
        }
        let v0 = prev.head(self);
        let v1 = next.tail(self);
        let enew = self.new_edge(v0, v1, prev, next, n0, p1)?;
        let (h, oh) = enew.halfedges();
        // Rewire halfedge -> halfedge.
        self.link_halfedges(prev, h);
        self.link_halfedges(h, next);
        self.link_halfedges(n0, oh);
        self.link_halfedges(oh, p1);
        // Rewire face -> halfedge and halfedge -> face.
        if let Some(f) = f {
            let fnew = self.new_face(oh)?;
            let hf = f.halfedge(self);
            self.halfedge_mut(h).face = Some(f);
            for (mesh, h) in self.loop_ccw_iter_mut(oh) {
                if hf == h {
                    mesh.face_mut(f).halfedge = h;
                }
                mesh.halfedge_mut(h).face = Some(fnew);
            }
        } else {
            let fnew = self.new_face(h)?;
            for (mesh, h) in self.loop_ccw_iter_mut(h) {
                mesh.halfedge_mut(h).face = Some(fnew);
            }
        };
        self.adjust_outgoing_halfedge(v0);
        self.adjust_outgoing_halfedge(v1);
        Ok(enew)
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    /// Check if it is safe to collapse an edge.
    ///
    /// This function only checks for topological errors, and doesn't account
    /// for shape and other geometric properties. This function uses the
    /// edge-status and vertex-status(mutable) properties to keep track of the
    /// neighborhood. Use this with borrowed properties as arguments if you're
    /// calling this function in a hot loop, to avoid repeated borrows.
    pub fn check_edge_collapse(&self, h: HH, estatus: &[Status], vstatus: &mut [Status]) -> bool {
        self.topol.check_edge_collapse(h, estatus, vstatus)
    }

    /// This is the same as [`Self::check_edge_collapse`], except it will
    /// attempt to borrow the necessary properties and may return an error if it
    /// ccannot borrow the required properties.
    pub fn try_check_edge_collapse(&self, h: HH) -> Result<bool, Error> {
        self.topol.try_check_edge_collapse(h)
    }

    /// Collapse the given edge.
    ///
    /// The vertex at the start of the halfedge will be deleted, and all the
    /// topology associated with that vertex will be rewired to the vertex that
    /// the halfedge is pointing towards. This function requires borrowed
    /// status properties of vertices, halfedges, edges and faces. This is
    /// useful when collapsing edges in a hot loop, to avoid repeated borrows of
    /// properties.
    pub fn collapse_edge(
        &mut self,
        h: HH,
        vstatus: &mut [Status],
        hstatus: &mut [Status],
        estatus: &mut [Status],
        fstatus: &mut [Status],
    ) {
        self.topol
            .collapse_edge(h, vstatus, hstatus, estatus, fstatus, &mut self.cache)
    }

    /// This is the same as [`Self::collapse_edge`]. Except this function will
    /// try to borrow all the necessary properties, and return an error if the
    /// borrowing fails.
    pub fn try_collapse_edge(&mut self, h: HH) -> Result<(), Error> {
        self.topol.try_collapse_edge(h, &mut self.cache)
    }

    /// Triangulate a face.
    ///
    /// This does not take the geometry / shape of the face into account. This
    /// only accounts for the topology of the face.
    pub fn triangulate_face(&mut self, f: FH) -> Result<(), Error> {
        self.topol.triangulate_face(f)
    }

    /// Triangulate all faces in this mesh.
    ///
    /// This does not take the geometry / shape of the face into account. This
    /// only accounts for the topology of the face.
    pub fn triangulate(&mut self) -> Result<(), Error> {
        self.topol.triangulate()
    }

    /// Split an edge with a new vertex at the given position.
    ///
    /// A new vertex is inserted at the given position and is used to split the
    /// given edge. A new edge is created during this split. If successful, a
    /// tuple containing the new vertex and the new edge is returned.
    pub fn split_edge(
        &mut self,
        e: EH,
        pos: A::Vector,
        copy_props: bool,
    ) -> Result<(VH, EH), Error> {
        let v = self.add_vertex(pos)?;
        let enew = self.topol.split_edge(e, v, copy_props)?;
        Ok((v, enew))
    }

    /// Swap an edge counter-clockwise.
    ///
    /// If the edge is a boundary edge, or some other topological error is
    /// encountered, then mesh is unmodified and a `false` is
    /// returned. Otherwise a `true` is returned.
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, HasTopology, Handle, HasIterators};
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_tri_face(0.into(), 1.into(), 2.into()).expect("Cannot add face");
    /// mesh.add_tri_face(0.into(), 2.into(), 3.into()).expect("Cannot add face");
    /// assert_eq!(mesh.triangulated_vertices().flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [2, 0, 1, 3, 0, 2]);
    /// let e = mesh.find_halfedge(0.into(), 2.into())
    ///             .expect("Cannot find halfedge").edge();
    /// mesh.swap_edge_ccw(e);
    /// assert_eq!(mesh.triangulated_vertices().flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [3, 1, 2, 3, 0, 1]);
    /// ```
    pub fn swap_edge_ccw(&mut self, e: EH) -> bool {
        self.topol.swap_edge_ccw(e)
    }

    /// Swap an edge clockwise.
    ///
    /// If the edge is a boundary edge, or some other topological error is
    /// encountered, then mesh is unmodified and a `false` is
    /// returned. Otherwise a `true` is returned.
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, HasTopology, Handle, HasIterators};
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_tri_face(0.into(), 1.into(), 2.into()).expect("Cannot add face");
    /// mesh.add_tri_face(0.into(), 2.into(), 3.into()).expect("Cannot add face");
    /// assert_eq!(mesh.triangulated_vertices().flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [2, 0, 1, 3, 0, 2]);
    /// let e = mesh.find_halfedge(0.into(), 2.into())
    ///             .expect("Cannot find halfedge").edge();
    /// mesh.swap_edge_cw(e);
    /// assert_eq!(mesh.triangulated_vertices().flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [1, 3, 0, 3, 1, 2]);
    /// ```
    pub fn swap_edge_cw(&mut self, e: EH) -> bool {
        self.topol.swap_edge_cw(e)
    }

    /// Remove an edge and unite the two incident faces into one face.
    ///
    /// ```text
    ///     +---------+---------+             +---------+---------+
    ///     |         |         |             |                   |
    ///     |         |         |             |                   |
    ///     |    f0   e    f1   |     ====>   |         f0        |  (f1 is deleted)
    ///     |         |         |             |                   |
    ///     |         |         |             |                   |
    ///     +---------+---------+             +---------+---------+
    /// ```
    pub fn remove_edge(&mut self, e: EH) -> Result<FH, Error> {
        self.topol.remove_edge(e)
    }

    /// Insert a new edge panning the end of `prev` and the start of `next` to
    /// split the loop into two.
    ///
    /// ```text
    ///     <--------- <----------  <---------- <----------
    ///    |              next    ^|                       ^
    ///    |                      ||                       |
    ///    |                      ||                       |
    ///    |              new edge||                       |
    ///    |                      ||                       |
    ///    |                      ||                       |
    ///    v              prev    |v                       |
    ///     ---------> ---------->  ----------> ---------->
    /// ```
    ///
    /// The loop is split into two as shown in this diagram by inserting the new
    /// edge to span `prev` and `next`. The two halfedges must be part of the
    /// same loop. If this loop is a boundary loop, a new face is inserted in
    /// the loop containing `prev` and `next`. If the original loop contains a
    /// face, the face is split into two new faces. The existing face will
    /// remain valid and correspond to the loop containing `prev` and `next`.
    pub fn insert_edge(&mut self, prev: HH, next: HH) -> Result<EH, Error> {
        self.topol.insert_edge(prev, next)
    }
}

#[cfg(test)]
mod test {
    use crate::{
        element::Handle,
        iterator::HasIterators,
        topol::{
            test::{loop_mesh, quad_box},
            HasTopology, TopolCache,
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
                .fold((0usize, 0usize), |(t, q), f| match f.valence(&qbox) {
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
                        match f.valence(&qbox) {
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
            qbox.faces()
                .flat_map(|f| qbox.triangulated_face_vertices(f))
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
                match f.valence(&qbox) {
                    3 => (t + 1, q),
                    4 => (t, 1 + q),
                    _ => (t, q),
                }
            })
        );
        assert_eq!(13, qbox.num_edges());
        assert_eq!(26, qbox.num_halfedges());
        assert_eq!(
            qbox.fv_ccw_iter(5.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[5, 6, 7]
        );
        assert_eq!(
            qbox.fv_ccw_iter(6.into())
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
                .flat_map(|f| qbox.fv_ccw_iter(f))
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[
                3, 2, 1, 1, 5, 4, 2, 6, 5, 3, 7, 6, 0, 4, 7, 5, 6, 7, 0, 3, 1, 0, 1, 4, 1, 2, 5, 2,
                3, 6, 3, 0, 7, 4, 5, 7
            ]
        );
    }

    #[test]
    fn t_box_split_edge() {
        let mut qbox = quad_box();
        let e = qbox
            .find_halfedge(4.into(), 5.into())
            .expect("Cannot find halfedge")
            .edge();
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let v = qbox.add_vertex().expect("Cannotr add vertex");
        let enew = qbox.split_edge(e, v, false).expect("Cannot split edge");
        let hnew = enew.halfedge(false);
        let ohnew = enew.halfedge(true);
        assert_eq!(oh.head(&qbox), ohnew.tail(&qbox));
        assert_eq!(oh.head(&qbox), v);
        assert_eq!(hnew.head(&qbox), h.tail(&qbox));
        assert_eq!(hnew.head(&qbox), v);
        assert_eq!(
            ohnew.next(&qbox),
            qbox.find_halfedge(5.into(), 6.into())
                .expect("Cannot find halfedge")
        );
        assert_eq!(
            hnew.prev(&qbox),
            qbox.find_halfedge(1.into(), 5.into())
                .expect("Cannot find halfedge")
        );
        assert_eq!(h.prev(&qbox), hnew);
        assert_eq!(oh.next(&qbox), ohnew);
        assert_eq!(h.face(&qbox), hnew.face(&qbox));
        assert_eq!(oh.face(&qbox), ohnew.face(&qbox));
        qbox.check().expect("Topological errors found");
    }

    #[test]
    fn t_box_split_edge_copy_props() {
        let mut qbox = quad_box();
        let e = qbox
            .find_halfedge(5.into(), 6.into())
            .expect("Cannot find halfedge")
            .edge();
        // Set properties.
        let mut eprop = qbox.create_edge_prop::<usize>(0);
        eprop.set(e, 123).expect("Cannot set property");
        let mut hprop = qbox.create_halfedge_prop::<usize>(0);
        hprop
            .set(e.halfedge(true), 234)
            .expect("Cannot set property");
        hprop
            .set(e.halfedge(false), 345)
            .expect("Cannot set property");
        assert_eq!(
            1,
            qbox.edges()
                .filter(|e| eprop.get(*e).expect("Cannot read property") != usize::default())
                .count()
        );
        assert_eq!(
            2,
            qbox.halfedges()
                .filter(|h| hprop.get(*h).expect("Cannot read property") != usize::default())
                .count()
        );
        // Do the split and check topology.
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let v = qbox.add_vertex().expect("Cannotr add vertex");
        let enew = qbox.split_edge(e, v, true).expect("Cannot split edge");
        let hnew = enew.halfedge(false);
        let ohnew = enew.halfedge(true);
        assert_eq!(oh.head(&qbox), ohnew.tail(&qbox));
        assert_eq!(oh.head(&qbox), v);
        assert_eq!(hnew.head(&qbox), h.tail(&qbox));
        assert_eq!(hnew.head(&qbox), v);
        assert_eq!(h.prev(&qbox), hnew);
        assert_eq!(oh.next(&qbox), ohnew);
        assert_eq!(h.face(&qbox), hnew.face(&qbox));
        assert_eq!(oh.face(&qbox), ohnew.face(&qbox));
        // Check properties.
        assert_eq!(
            2,
            qbox.edges()
                .filter(|e| eprop.get(*e).expect("Cannot read property") != usize::default())
                .count()
        );
        assert_eq!(
            4,
            qbox.halfedges()
                .filter(|h| hprop.get(*h).expect("Cannot read property") != usize::default())
                .count()
        );
        let eprop = eprop.try_borrow().expect("Cannot borrow edge");
        assert_eq!(eprop[e.index() as usize], eprop[enew.index() as usize]);
        let hprop = hprop.try_borrow().expect("Cannot borrow halfedge");
        assert_eq!(hprop[h.index() as usize], hprop[hnew.index() as usize]);
        assert_eq!(hprop[oh.index() as usize], hprop[ohnew.index() as usize]);
        qbox.check().expect("Topological errors found");
    }

    #[test]
    fn t_box_swap_edge_ccw() {
        let mut qbox = quad_box();
        qbox.triangulate().expect("Cannot triangulate the mesh");
        let h = qbox
            .find_halfedge(5.into(), 7.into())
            .expect("Cannot find halfedge");
        let e = h.edge();
        assert!(qbox.swap_edge_ccw(e), "Cannot swap edge");
        assert_eq!(
            qbox.faces()
                .flat_map(|f| qbox.fv_ccw_iter(f))
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[
                3, 2, 1, 1, 5, 4, 2, 6, 5, 3, 7, 6, 0, 4, 7, 6, 7, 4, 0, 3, 1, 0, 1, 4, 1, 2, 5, 2,
                3, 6, 3, 0, 7, 4, 5, 6
            ]
        );
        qbox.check().expect("Topological errors found");
    }

    #[test]
    fn t_box_swap_edge_cw() {
        let mut qbox = quad_box();
        qbox.triangulate().expect("Cannot triangulate the mesh");
        let h = qbox
            .find_halfedge(5.into(), 7.into())
            .expect("Cannot find halfedge");
        let e = h.edge();
        assert!(qbox.swap_edge_cw(e), "Cannot swap edge");
        assert_eq!(
            qbox.faces()
                .flat_map(|f| qbox.fv_ccw_iter(f))
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            &[
                3, 2, 1, 1, 5, 4, 2, 6, 5, 3, 7, 6, 0, 4, 7, 4, 5, 6, 0, 3, 1, 0, 1, 4, 1, 2, 5, 2,
                3, 6, 3, 0, 7, 4, 6, 7
            ]
        );
        qbox.check().expect("Topological errors found");
    }

    #[test]
    fn t_box_remove_edge() {
        let mut cache = TopolCache::default();
        let mut qbox = quad_box();
        qbox.triangulate().expect("Cannot triangulate the mesh");
        qbox.remove_edge(
            qbox.find_halfedge(5.into(), 7.into())
                .expect("Cannot find halfedge")
                .edge(),
        )
        .expect("Cannot remove edge");
        qbox.garbage_collection(&mut cache)
            .expect("Cannot garbage collect");
        qbox.check().expect("Topological errors found");
        assert_eq!(11, qbox.num_faces());
        assert_eq!(17, qbox.num_edges());
        assert_eq!(8, qbox.num_vertices());
        assert_eq!(
            (10, 1),
            qbox.faces().fold((0usize, 0usize), |(tris, quads), f| {
                match f.valence(&qbox) {
                    3 => (tris + 1, quads),
                    4 => (tris, quads + 1),
                    _ => (tris, quads),
                }
            })
        );
        qbox.remove_edge(
            qbox.find_halfedge(1.into(), 4.into())
                .expect("Cannot find halfedge")
                .edge(),
        )
        .expect("Cannot remove edge");
        qbox.garbage_collection(&mut cache)
            .expect("Cannot garbage collect");
        qbox.check().expect("Topological errors found");
        assert_eq!(10, qbox.num_faces());
        assert_eq!(16, qbox.num_edges());
        assert_eq!(8, qbox.num_vertices());
        assert_eq!(
            (8, 2),
            qbox.faces().fold((0usize, 0usize), |(tris, quads), f| {
                match f.valence(&qbox) {
                    3 => (tris + 1, quads),
                    4 => (tris, quads + 1),
                    _ => (tris, quads),
                }
            })
        );
    }

    #[test]
    fn t_loop_mesh_insert_edge_in_hole() {
        let mut mesh = loop_mesh();
        let e = mesh
            .insert_edge(
                mesh.find_halfedge(5.into(), 6.into())
                    .expect("Cannot find halfedge"),
                mesh.find_halfedge(9.into(), 5.into())
                    .expect("Cannot find halfedge"),
            )
            .expect("Cannot insert halfedge");
        let (h, oh) = e.halfedges();
        mesh.check().expect("Topological errors found");
        assert!(oh.is_boundary(&mesh));
        assert!(!h.is_boundary(&mesh));
        assert_eq!(3, mesh.loop_ccw_iter(h).count());
        assert_eq!(3, mesh.loop_ccw_iter(oh).count());
        assert_eq!(
            (1, 8),
            mesh.faces().fold((0usize, 0usize), |(tris, quads), f| {
                match f.valence(&mesh) {
                    3 => (tris + 1, quads),
                    4 => (tris, quads + 1),
                    _ => (tris, quads),
                }
            })
        );
    }

    #[test]
    fn t_box_insert_edge() {
        let mut mesh = quad_box();
        let e = mesh
            .insert_edge(
                mesh.find_halfedge(4.into(), 5.into())
                    .expect("Cannot find halfedge"),
                mesh.find_halfedge(7.into(), 4.into())
                    .expect("Cannot find halfedge"),
            )
            .expect("Cannot insert edge");
        let (h, oh) = e.halfedges();
        mesh.check().expect("Topological errors found");
        assert!(!h.is_boundary(&mesh));
        assert!(!oh.is_boundary(&mesh));
        assert_eq!(3, mesh.loop_ccw_iter(h).count());
        assert_eq!(3, mesh.loop_ccw_iter(oh).count());
    }
}
