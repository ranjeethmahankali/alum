use crate::{
    Adaptor, EPropBuf, FPropBuf, HPropBuf, HasTopology, PolyMeshT, VPropBuf,
    element::{EH, FH, HH, VH},
    error::Error,
    iterator::HasIterators,
    status::Status,
    topol::Topology,
};

/// This trait defines functions to edit the topology of a mesh, with typical
/// editing operations such as edge swaps, collapses etc.
pub trait EditableTopology: HasIterators {
    /// Check if it is safe to collapse an edge.
    ///
    /// This function only checks for topological errors, and doesn't account
    /// for shape and other geometric properties. This function uses the
    /// edge-status and vertex-status(mutable) properties to keep track of the
    /// neighborhood. Use this with borrowed properties as arguments if you're
    /// calling this function in a hot loop, to avoid repeated borrows.
    fn check_edge_collapse(
        &self,
        h: HH,
        edge_status: &EPropBuf<Status>,
        vertex_status: &mut VPropBuf<Status>,
    ) -> bool {
        // Check if already deleted.
        if edge_status[h.edge()].deleted() {
            return false;
        }
        let oh = h.opposite();
        let v0 = oh.head(self);
        let v1 = h.head(self);
        if vertex_status[v0].deleted() || vertex_status[v1].deleted() {
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
        let v0_boundary = v0.is_boundary(self);
        let v1_boundary = v1.is_boundary(self);
        // Check if we're collapsing across two different boundaries.
        if v0_boundary && v1_boundary && !h.is_boundary(self) && !oh.is_boundary(self) {
            return false;
        }
        // A boundary vertex should not be collapsed inward.
        if v0_boundary && !v1_boundary {
            return false;
        }
        // Check the 'Link condition' Edelsbrunner [2006]. The intersection of
        // the one rings of from and to vertices must be the left and right
        // vertices, and only if the corresponding faces are triangles.
        let vl = h.next(self).head(self);
        let vr = oh.next(self).head(self);
        for v in self.vv_ccw_iter(v0) {
            vertex_status[v].set_tagged(false);
        }
        for v in self.vv_ccw_iter(v1) {
            vertex_status[v].set_tagged(true);
        }
        for v in self.vv_ccw_iter(v0) {
            if vertex_status[v].tagged() && !(v == vl && htriangle) && !(v == vr && ohtriangle) {
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
            if vertex_status[h.head(self)].tagged() && vl == vr && htriangle && ohtriangle {
                return false;
            }
        }
        true
    }

    /// This is the same as [`Self::check_edge_collapse`], except it will
    /// attempt to borrow the necessary properties and may return an error if it
    /// ccannot borrow the required properties.
    fn try_check_edge_collapse(&self, h: HH) -> Result<bool, Error> {
        let topol = self.topology();
        let estatus = topol.estatus.try_borrow()?;
        // Clone the property to side step compile time borrow checker. The
        // runtime borrow checker is still in use, so not a problem.
        let mut vstatus = topol.vstatus.clone();
        let mut vstatus = vstatus.try_borrow_mut()?;
        Ok(topol.check_edge_collapse(h, &estatus, &mut vstatus))
    }

    /// Sometimes after collapsing an edge, if the neighboring faces are
    /// triangles, we end up with degenerate loops / faces. This cleans up such
    /// loops.
    fn collapse_degenerate_triangle(
        &mut self,
        h: HH,
        hstatus: &mut HPropBuf<Status>,
        estatus: &mut EPropBuf<Status>,
        fstatus: &mut FPropBuf<Status>,
    ) {
        let topol = self.topology_mut();
        let h1 = h.next(topol);
        let o = h.opposite();
        let o1 = h1.opposite();
        let on = o.next(topol);
        let op = o.prev(topol);
        let v0 = h.head(topol);
        let v1 = h1.head(topol);
        let fh = h.face(topol);
        let fo = o.face(topol);
        // Ensure the loop represents a collapsed triangle. Because this is a
        // private function, all callers inside the implementation, so we can be
        // confident and assert to catch and weed out any bugs in debug builds.
        debug_assert_eq!(h1.next(topol), h);
        debug_assert_ne!(h1, o);
        // Rewire halfedge -> halfedge.
        topol.link_halfedges(h1, on);
        topol.link_halfedges(op, h1);
        // Rewire halfedge -> face.
        topol.halfedge_mut(h1).face = fo;
        // Rewire vertex -> halfedge.
        topol.vertex_mut(v0).halfedge = Some(h1);
        topol.adjust_outgoing_halfedge(v0);
        topol.vertex_mut(v1).halfedge = Some(o1);
        topol.adjust_outgoing_halfedge(v1);
        // Rewire face -> halfedge.
        if let Some(fo) = fo {
            if fo.halfedge(topol) == o {
                topol.face_mut(fo).halfedge = h1;
            }
        }
        // Delete stuff.
        if let Some(fh) = fh {
            fstatus[fh].set_deleted(true);
        }
        estatus[h.edge()].set_deleted(true);
        hstatus[h].set_deleted(true);
        hstatus[o].set_deleted(true);
    }

    /// Collapse the given edge.
    ///
    /// The vertex at the start of the halfedge will be deleted, and all the
    /// topology associated with that vertex will be rewired to the vertex that
    /// the halfedge is pointing towards. This function requires borrowed status
    /// properties of vertices, halfedges, edges and faces. This is useful when
    /// collapsing edges in a hot loop, to avoid repeated borrows of
    /// properties. `hcache` is used for temporary storage. Reusing the same
    /// `hcache` for many collapses avoids repeated allocations.
    fn collapse_edge(
        &mut self,
        h: HH,
        vstatus: &mut VPropBuf<Status>,
        hstatus: &mut HPropBuf<Status>,
        estatus: &mut EPropBuf<Status>,
        fstatus: &mut FPropBuf<Status>,
        hcache: &mut Vec<HH>,
    ) {
        let topol = self.topology_mut();
        // Collect neighboring topology.
        let hn = h.next(topol);
        let hp = h.prev(topol);
        let o = h.opposite();
        let on = o.next(topol);
        let op = o.prev(topol);
        let fh = h.face(topol);
        let fo = o.face(topol);
        let vh = h.head(topol);
        let vo = o.head(topol);
        // Rewire halfedge -> vertex
        hcache.clear();
        hcache.extend(topol.vih_ccw_iter(vo));
        for ih in hcache.drain(..) {
            topol.halfedge_mut(ih).vertex = vh;
        }
        // Rewire halfedge -> halfedge
        topol.link_halfedges(hp, hn);
        topol.link_halfedges(op, on);
        // Rewire face -> halfedge
        if let Some(fh) = fh {
            topol.face_mut(fh).halfedge = hn;
        }
        if let Some(fo) = fo {
            topol.face_mut(fo).halfedge = on;
        }
        // Rewire vertex -> halfedge
        if vh.halfedge(topol) == Some(o) {
            topol.vertex_mut(vh).halfedge = Some(hn);
        }
        topol.adjust_outgoing_halfedge(vh);
        topol.vertex_mut(vo).halfedge = None;
        // Delete stuff
        estatus[h.edge()].set_deleted(true);
        vstatus[vo].set_deleted(true);
        hstatus[h].set_deleted(true);
        hstatus[o].set_deleted(true);
        // If the loops that used to contain the halfedges that were collapsed
        // and deleted had a valance of 3, they are now degenerate. So we need
        // to collapse those loops.
        if hn.next(topol) == hp {
            topol.collapse_degenerate_triangle(hn, hstatus, estatus, fstatus);
        }
        if on.next(topol) == op {
            topol.collapse_degenerate_triangle(on, hstatus, estatus, fstatus);
        }
    }

    /// This is the same as [`Self::collapse_edge`]. Except this function will
    /// try to borrow all the necessary properties, and return an error if the
    /// borrowing fails.
    fn try_collapse_edge(&mut self, h: HH, hcache: &mut Vec<HH>) -> Result<(), Error> {
        let topol = self.topology_mut();
        let mut vstatus = topol.vstatus.clone();
        let mut vstatus = vstatus.try_borrow_mut()?;
        let mut hstatus = topol.hstatus.clone();
        let mut hstatus = hstatus.try_borrow_mut()?;
        let mut estatus = topol.estatus.clone();
        let mut estatus = estatus.try_borrow_mut()?;
        let mut fstatus = topol.fstatus.clone();
        let mut fstatus = fstatus.try_borrow_mut()?;
        topol.collapse_edge(
            h,
            &mut vstatus,
            &mut hstatus,
            &mut estatus,
            &mut fstatus,
            hcache,
        );
        Ok(())
    }

    /// Triangulate a face.
    ///
    /// This does not take the geometry / shape of the face into account. This
    /// only accounts for the topology of the face.
    fn triangulate_face(&mut self, f: FH) -> Result<(), Error> {
        let topol = self.topology_mut();
        let mut base = f.halfedge(topol);
        let vstart = base.tail(topol);
        let prev = base.prev(topol);
        let mut next = base.next(topol);
        while next.next(topol).head(topol) != vstart {
            let next2 = next.next(topol);
            let fnew = topol.new_face(base)?;
            let enew = topol.new_edge(vstart, next.head(topol), prev, next2, next, base)?;
            let hnew = enew.halfedge(false);
            let ohnew = enew.halfedge(true);
            // Link the triangle created.
            topol.link_halfedges(base, next);
            topol.link_halfedges(next, ohnew);
            topol.link_halfedges(ohnew, base);
            // Set face handles.
            topol.halfedge_mut(base).face = Some(fnew);
            topol.halfedge_mut(next).face = Some(fnew);
            topol.halfedge_mut(ohnew).face = Some(fnew);
            // Copy properties.
            topol.hprops.copy(prev, ohnew)?;
            topol.hprops.copy(prev, hnew)?;
            topol.fprops.copy(f, fnew)?;
            // For next iteration.
            base = hnew;
            next = next2;
        }
        // Last face takes the original face handle.
        topol.face_mut(f).halfedge = base;
        topol.link_halfedges(base, next);
        topol.link_halfedges(next.next(topol), base);
        topol.halfedge_mut(base).face = Some(f);
        Ok(())
    }

    /// Triangulate all faces in this mesh.
    ///
    /// This does not take the geometry / shape of the face into account. This
    /// only accounts for the topology of the face.
    fn triangulate(&mut self) -> Result<(), Error> {
        for f in self.faces() {
            self.triangulate_face(f)?;
        }
        Ok(())
    }

    /// Split an edge with a new vertex at the given position.
    ///
    /// A new vertex is inserted at the given position and is used to split the
    /// edge. Say, the existing edge spans vertices `a` and `b`. After the
    /// split, this edge will span vertices `v` and `b`, and a new edge is
    /// created that spans vertices `a` and `v`. If successful, a halfedge
    /// pointing from `a` to `v` is returned.
    fn split_edge(&mut self, h: HH, v: VH, copy_props: bool) -> Result<HH, Error> {
        let topol = self.topology_mut();
        let oh = h.opposite();
        let e = h.edge();
        let vfrom = h.tail(topol);
        let (ph, on) = (h.prev(topol), oh.next(topol));
        let (f, of) = (h.face(topol), oh.face(topol));
        // Create a new edge and rewire topology.
        let enew = topol.new_edge(vfrom, v, ph, h, oh, on)?;
        let hnew = enew.halfedge(false);
        let ohnew = enew.halfedge(true);
        // Rewire halfedge -> vertex.
        topol.halfedge_mut(oh).vertex = v;
        // Rewire halfedge -> halfedge.
        topol.link_halfedges(hnew, h);
        topol.link_halfedges(oh, ohnew);
        topol.link_halfedges(ph, hnew);
        topol.link_halfedges(ohnew, on);
        // Rewire halfedge -> face.
        topol.halfedge_mut(hnew).face = f;
        topol.halfedge_mut(ohnew).face = of;
        // Rewire vertex -> halfedge.
        topol.vertex_mut(v).halfedge = Some(h);
        topol.adjust_outgoing_halfedge(v);
        if vfrom.halfedge(topol) == Some(h) {
            topol.vertex_mut(vfrom).halfedge = Some(hnew);
            topol.adjust_outgoing_halfedge(vfrom);
        }
        if copy_props {
            topol.eprops.copy(e, enew)?;
            topol.hprops.copy_many(&[h, oh], &[hnew, ohnew])?;
        }
        Ok(hnew)
    }

    /// Split the face by connecting all incident vertices to the given vertex.
    ///
    /// This will triangulate the faces by connecting all the incident vertices
    /// to the given vertex. The given vertex must be isolated to produce valid
    /// topology, otherwise [`Error::ComplexVertex`] is returned.
    fn split_face(&mut self, f: FH, v: VH, copy_props: bool) -> Result<(), Error> {
        let topol = self.topology_mut();
        // After we're done, this vertex will be in the interior of the
        // face. The vertex must be isolated. Which means it must be isolated
        // before we start.
        if v.halfedge(topol).is_some() {
            return Err(Error::ComplexVertex(v));
        }
        let valence = f.valence(topol) as u32;
        let hend = f.halfedge(topol);
        let num_edges_before = topol.num_edges() as u32;
        let mut ei = 0u32;
        let mut hh = hend.next(topol);
        let hold = {
            // Predetermining the indices of edges we haven't created yet. We will create them later.
            let enext: EH = (num_edges_before + ((ei + 1) % valence)).into();
            let eprev: EH = (num_edges_before + ((ei + valence - 1) % valence)).into();
            let enew = topol.new_edge(
                hend.head(topol),
                v,
                hend,
                eprev.halfedge(true),
                enext.halfedge(false),
                hh,
            )?;
            ei += 1;
            enew.halfedge(false)
        };
        topol.link_halfedges(hend, hold);
        topol.halfedge_mut(hold).face = Some(f);
        let mut hold = hold.opposite();
        while hh != hend {
            let hnext = hh.next(topol);
            let fnew = topol.new_face(hh)?;
            // Predetermining the indices of edges we may not have created yet.
            let enext: EH = (num_edges_before + ((ei + 1) % valence)).into();
            let enew = topol.new_edge(hh.head(topol), v, hh, hold, enext.halfedge(false), hnext)?;
            ei += 1;
            let hnew = enew.halfedge(false);
            topol.link_halfedges(hnew, hold);
            topol.link_halfedges(hold, hh);
            topol.link_halfedges(hh, hnew);
            topol.halfedge_mut(hnew).face = Some(fnew);
            topol.halfedge_mut(hold).face = Some(fnew);
            topol.halfedge_mut(hh).face = Some(fnew);
            hold = hnew.opposite();
            hh = hnext;
        }
        topol.link_halfedges(hold, hend);
        topol.link_halfedges(hend.next(topol), hold);
        topol.halfedge_mut(hold).face = Some(f);
        topol.vertex_mut(v).halfedge = Some(hold);
        if copy_props {
            for (mesh, fnew) in topol.vf_ccw_iter_mut(v).filter(|(_m, fnew)| *fnew != f) {
                mesh.fprops.copy(f, fnew)?;
            }
        }
        Ok(())
    }

    /// Swap an edge counter-clockwise.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology, Handle, HasIterators, EditableTopology};
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_tri_face(0.into(), 1.into(), 2.into()).expect("Cannot add face");
    /// mesh.add_tri_face(0.into(), 2.into(), 3.into()).expect("Cannot add face");
    /// let fstatus = mesh.face_status_prop();
    /// let fstatus = fstatus.try_borrow().unwrap();
    /// assert_eq!(mesh.triangulated_vertices(&fstatus).flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [2, 0, 1, 3, 0, 2]);
    /// let e = mesh.find_halfedge(0.into(), 2.into())
    ///             .expect("Cannot find halfedge").edge();
    /// mesh.swap_edge_ccw(e);
    /// assert_eq!(mesh.triangulated_vertices(&fstatus).flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [3, 1, 2, 3, 0, 1]);
    /// ```
    fn swap_edge_ccw(&mut self, e: EH) -> Result<(), Error> {
        let topol = self.topology_mut();
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let (f, of) = match (h.face(topol), oh.face(topol)) {
            (Some(f), Some(of)) => (f, of),
            _ => return Err(Error::CannotSwapBoundaryEdge(e)), // Cannot swap boundary edge.
        };
        let hn = h.next(topol);
        let on = oh.next(topol);
        let v0 = h.tail(topol);
        let v1 = h.head(topol);
        // Check for degeneracy.
        if f == of || hn == oh || hn.head(topol) == v0 || on == h || on.head(topol) == v1 {
            return Err(Error::DegenerateEdge(e));
        }
        let hp = h.prev(topol);
        let op = oh.prev(topol);
        let hnn = hn.next(topol);
        let onn = on.next(topol);
        let hnv = hn.head(topol);
        let onv = on.head(topol);
        // Rewire vertex -> halfedge.
        if v0.halfedge(topol) == Some(h) {
            topol.vertex_mut(v0).halfedge = Some(on);
        }
        if v1.halfedge(topol) == Some(oh) {
            topol.vertex_mut(v1).halfedge = Some(hn);
        }
        // Rewire halfedge -> vertex.
        topol.halfedge_mut(h).vertex = hnv;
        topol.halfedge_mut(oh).vertex = onv;
        // Rewire halfedge -> halfedge.
        topol.link_halfedges(oh, onn);
        topol.link_halfedges(op, hn);
        topol.link_halfedges(hn, oh);
        topol.link_halfedges(h, hnn);
        topol.link_halfedges(hp, on);
        topol.link_halfedges(on, h);
        // Rewire halfedge -> face.
        topol.halfedge_mut(hn).face = Some(of);
        topol.halfedge_mut(on).face = Some(f);
        // Rewire face -> halfedge.
        if f.halfedge(topol) == hn {
            topol.face_mut(f).halfedge = h;
        }
        if of.halfedge(topol) == on {
            topol.face_mut(of).halfedge = oh;
        }
        Ok(())
    }

    /// Swap an edge clockwise.
    ///
    /// If the edge is a boundary edge, or some other topological error is
    /// encountered, then mesh is unmodified and a `false` is
    /// returned. Otherwise a `true` is returned.
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology, Handle, HasIterators, EditableTopology};
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_tri_face(0.into(), 1.into(), 2.into()).expect("Cannot add face");
    /// mesh.add_tri_face(0.into(), 2.into(), 3.into()).expect("Cannot add face");
    /// let fstatus = mesh.face_status_prop();
    /// let fstatus = fstatus.try_borrow().unwrap();
    /// assert_eq!(mesh.triangulated_vertices(&fstatus).flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [2, 0, 1, 3, 0, 2]);
    /// let e = mesh.find_halfedge(0.into(), 2.into())
    ///             .expect("Cannot find halfedge").edge();
    /// mesh.swap_edge_cw(e);
    /// assert_eq!(mesh.triangulated_vertices(&fstatus).flatten().map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [1, 3, 0, 3, 1, 2]);
    /// ```
    fn swap_edge_cw(&mut self, e: EH) -> bool {
        let topol = self.topology_mut();
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let (f, of) = match (h.face(topol), oh.face(topol)) {
            (Some(f), Some(of)) => (f, of),
            _ => return false, // Cannot swap boundary edge.
        };
        let hn = h.next(topol);
        let on = oh.next(topol);
        let v0 = h.tail(topol);
        let v1 = h.head(topol);
        // Check for degeneracy.
        if f == of || hn == oh || hn.head(topol) == v0 || on == h || on.head(topol) == v1 {
            return false;
        }
        let hp = h.prev(topol);
        let op = oh.prev(topol);
        let hpp = hp.prev(topol);
        let opp = op.prev(topol);
        let hpv = hp.tail(topol);
        let opv = op.tail(topol);
        // Rewire vertex -> halfedge.
        if v0.halfedge(topol) == Some(h) {
            topol.vertex_mut(v0).halfedge = Some(on);
        }
        if v1.halfedge(topol) == Some(oh) {
            topol.vertex_mut(v1).halfedge = Some(hn);
        }
        // Rewire halfedge -> vertex.
        topol.halfedge_mut(h).vertex = opv;
        topol.halfedge_mut(oh).vertex = hpv;
        // Rewire halfedge -> halfedge.
        topol.link_halfedges(h, op);
        topol.link_halfedges(op, hn);
        topol.link_halfedges(hpp, h);
        topol.link_halfedges(oh, hp);
        topol.link_halfedges(hp, on);
        topol.link_halfedges(opp, oh);
        // Rewire halfedge -> face.
        topol.halfedge_mut(op).face = Some(f);
        topol.halfedge_mut(hp).face = Some(of);
        // Rewire face -> halfedge.
        if f.halfedge(topol) == hp {
            topol.face_mut(f).halfedge = h;
        }
        if of.halfedge(topol) == op {
            topol.face_mut(of).halfedge = oh;
        }
        true
    }

    /// An edge is a unique link if it is the only edge connecting the two faces
    /// incident on it.
    ///
    /// The boundary is treated as one face. So a boundary edge can be a simple
    /// links only if it is the lone boundary edge, incident on the opposite
    /// face.
    fn edge_is_unique_link(&self, e: EH) -> bool {
        let h = e.halfedge(false);
        let fo = h.opposite().face(self);
        self.loop_ccw_iter(h)
            .skip(1)
            .all(|h| h.opposite().face(self) != fo)
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
    fn remove_edge(
        &mut self,
        e: EH,
        hstatus: &mut HPropBuf<Status>,
        estatus: &mut EPropBuf<Status>,
        fstatus: &mut FPropBuf<Status>,
    ) -> Result<FH, Error> {
        //    <--------- <----------v0<---------- <----------
        //   |                n0    ^|     p1                ^
        //   |                      ||                       |
        //   |                      ||                       |
        //   |         f0         h0||h1         f1          |
        //   |                      ||                       |
        //   |                      ||                       |
        //   v                p0    |v     n1                |
        //    ---------> ---------->v1----------> ---------->
        let topol = self.topology_mut();
        if estatus[e].deleted() {
            return Err(Error::DeletedEdge(e));
        }
        if !topol.edge_is_unique_link(e) {
            return Err(Error::EdgeIsNotAUniqueLink(e));
        }
        let (h0, h1) = e.halfedges();
        let (f0, f1) = match (h0.face(topol), h1.face(topol)) {
            (Some(f0), Some(f1)) => (f0, f1),
            _ => return Err(Error::CannotRemoveBoundaryEdge(e)),
        };
        let (p0, p1) = (h0.prev(topol), h1.prev(topol));
        let (n0, n1) = (h0.next(topol), h1.next(topol));
        let (v0, v1) = (h0.head(topol), h1.head(topol));
        // Rewire vertex -> halfedge.
        if v0.halfedge(topol) == Some(h1) {
            topol.vertex_mut(v0).halfedge = Some(n0);
        }
        if v1.halfedge(topol) == Some(h0) {
            topol.vertex_mut(v1).halfedge = Some(n1);
        }
        // Rewire halfedge -> halfedge.
        topol.link_halfedges(p0, n1);
        topol.link_halfedges(p1, n0);
        // Rewire face -> halfedge. Keep f0 and delete f1.
        if f0.halfedge(topol) == h0 {
            topol.face_mut(f0).halfedge = p0;
        }
        // Rewire halfedge -> face for the loop of f1.
        for (mesh, h) in topol.loop_ccw_iter_mut(n1).take_while(|(_, h)| *h != n0) {
            mesh.halfedge_mut(h).face = Some(f0);
        }
        estatus[e].set_deleted(true);
        hstatus[h0].set_deleted(true);
        hstatus[h1].set_deleted(true);
        fstatus[f1].set_deleted(true);
        Ok(f0)
    }

    /// Remove an edge and unite the two incident faces into one face. This is
    /// the same as `Self::remove_edge` except it will internally try to borrow
    /// the required properties without the caller having to supply them.
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
    fn try_remove_edge(&mut self, e: EH) -> Result<FH, Error> {
        let topol = self.topology_mut();
        let mut estatus = topol.estatus.clone();
        let mut estatus = estatus.try_borrow_mut()?;
        let mut fstatus = topol.fstatus.clone();
        let mut fstatus = fstatus.try_borrow_mut()?;
        let mut hstatus = topol.hstatus.clone();
        let mut hstatus = hstatus.try_borrow_mut()?;
        self.remove_edge(e, &mut hstatus, &mut estatus, &mut fstatus)
    }

    /// Insert a new edge panning the end of `prev` and the start of `next` to
    /// split the loop into two.
    ///
    /// ```text
    ///     <--------- <----------  <---------- <----------
    ///    |              next    ^|                       ^
    ///    |                      ||                       |
    ///    |                      ||                       |
    ///    |              new-edge||                       |
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
    /// face, the face is split into two faces. The existing face will remain
    /// valid and correspond to the loop containing `prev` and `next`, and a new
    /// face is created for the right loop.
    ///
    /// So in both cases, i.e. whether or not the loop is a boundary loop, or
    /// has a face, 1 new face is created during this operation. This can be
    /// avoided by optionally supplying a `newface`. If `newface` has `Some`
    /// face in it, that will be reused instead of creating a new face. The face
    /// supplied for reuse MUST be a deleted face, otherwise an error is
    /// returned.
    fn insert_edge(&mut self, prev: HH, next: HH, newface: Option<FH>) -> Result<EH, Error> {
        //    <--------- <----------v1<---------- <----------
        //   |                next  ^|     n0                ^
        //   |                      ||                       |
        //   |                      ||                       |
        //   |          f         h0||h1       fnew          |
        //   |                      ||                       |
        //   |                      ||                       |
        //   v                prev  |v     p1                |
        //    ---------> ---------->v0----------> ---------->
        let topol = self.topology_mut();
        let (p1, n0) = (prev.next(topol), next.prev(topol));
        if p1 == next || prev == next {
            return Err(Error::CannotInsertEdge(prev, next));
        }
        let f = prev.face(topol);
        if f != next.face(topol) {
            return Err(Error::HalfedgesNotInTheSameLoop(prev, next));
        }
        // Check to make sure the new face, if provided is a deleted face.
        if let Some(fnew) = newface {
            let mut fstatus = topol.fstatus.try_borrow_mut()?;
            if !fstatus[fnew].deleted() {
                return Err(Error::NotDeletedFace(fnew));
            }
            fstatus[fnew].set_deleted(false);
        }
        if f.is_none() {
            // Check if the halfedges are part of the same boundary loop. March
            // simultaenously starting from both prev and next halfedges, and
            // see if you arrive at the other halfedge. Marching from both
            // ensures we'll detect the loop as soon as possible.
            if !topol
                .loop_ccw_iter(prev)
                .zip(topol.loop_cw_iter(next))
                .any(|(n, p)| n == prev || p == next)
            {
                return Err(Error::HalfedgesNotInTheSameLoop(prev, next));
            }
        }
        let v0 = prev.head(topol);
        let v1 = next.tail(topol);
        let enew = topol.new_edge(v0, v1, prev, next, n0, p1)?;
        let (h0, h1) = enew.halfedges();
        // Rewire halfedge -> halfedge.
        topol.link_halfedges(prev, h0);
        topol.link_halfedges(h0, next);
        topol.link_halfedges(n0, h1);
        topol.link_halfedges(h1, p1);
        // Rewire face -> halfedge and halfedge -> face.
        match f {
            Some(f) => {
                let fnew = match newface {
                    Some(fnew) => {
                        topol.face_mut(fnew).halfedge = h1;
                        fnew
                    }
                    None => topol.new_face(h1)?,
                };
                let hf = f.halfedge(topol);
                topol.halfedge_mut(h0).face = Some(f);
                for (mesh, h) in topol.loop_ccw_iter_mut(h1) {
                    if hf == h {
                        mesh.face_mut(f).halfedge = h0;
                    }
                    mesh.halfedge_mut(h).face = Some(fnew);
                }
            }
            None => {
                let fnew = match newface {
                    Some(fnew) => {
                        topol.face_mut(fnew).halfedge = h0;
                        fnew
                    }
                    None => topol.new_face(h0)?,
                };
                for (mesh, h) in topol.loop_ccw_iter_mut(h0) {
                    mesh.halfedge_mut(h).face = Some(fnew);
                }
            }
        }
        // Adjust vertex halfedges.
        topol.adjust_outgoing_halfedge(v0);
        topol.adjust_outgoing_halfedge(v1);
        Ok(enew)
    }
}

impl<const DIM: usize, A> EditableTopology for PolyMeshT<DIM, A> where A: Adaptor<DIM> {}

impl EditableTopology for Topology {}

#[cfg(test)]
mod test {
    use crate::{
        HH,
        edit::EditableTopology,
        element::Handle,
        iterator::HasIterators,
        topol::{
            HasTopology, TopolCache,
            test::{loop_mesh, quad_box},
        },
        use_glam::PolyMeshF32,
    };

    #[test]
    fn t_box_check_edge_collapse() {
        let qbox = quad_box();
        for h in qbox.halfedges() {
            assert!(
                qbox.try_check_edge_collapse(h)
                    .expect("Cannot check halfedge collapse")
            );
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
        let mut hcache = Vec::new();
        let h = qbox
            .find_halfedge(5.into(), 6.into())
            .expect("Cannot find halfedge");
        qbox.try_collapse_edge(h, &mut hcache)
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
                    if estatus[e].deleted() {
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
                    if hstatus[h].deleted() {
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
        qbox.try_collapse_edge(h, &mut cache.halfedges)
            .expect("Cannot collapse edges");
        let h = qbox
            .find_halfedge(4.into(), 7.into())
            .expect("Cannot find halfedge");
        qbox.try_collapse_edge(h, &mut cache.halfedges)
            .expect("Cannot collapse edge");
        {
            let estatus = qbox
                .estatus
                .try_borrow()
                .expect("Cannot borrow the edge status property");
            assert_eq!(
                (3, 9),
                qbox.edges().fold((0usize, 0usize), |(del, ndel), e| {
                    if estatus[e].deleted() {
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
                    if hstatus[h].deleted() {
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
                    if fstatus[h].deleted() {
                        (del + 1, ndel)
                    } else {
                        (del, 1 + ndel)
                    }
                })
            );
            assert_eq!(
                (2, 3),
                qbox.faces().fold((0usize, 0usize), |(t, q), f| {
                    if fstatus[f].deleted() {
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
        let hnew = qbox.split_edge(h, v, false).expect("Cannot split edge");
        let ohnew = hnew.opposite();
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
                .filter(|e| eprop.get_cloned(*e).expect("Cannot read property") != usize::default())
                .count()
        );
        assert_eq!(
            2,
            qbox.halfedges()
                .filter(|h| hprop.get_cloned(*h).expect("Cannot read property") != usize::default())
                .count()
        );
        // Do the split and check topology.
        let h = e.halfedge(false);
        let oh = e.halfedge(true);
        let v = qbox.add_vertex().expect("Cannotr add vertex");
        let hnew = qbox.split_edge(h, v, true).expect("Cannot split edge");
        let ohnew = hnew.opposite();
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
                .filter(|e| eprop.get_cloned(*e).expect("Cannot read property") != usize::default())
                .count()
        );
        assert_eq!(
            4,
            qbox.halfedges()
                .filter(|h| hprop.get_cloned(*h).expect("Cannot read property") != usize::default())
                .count()
        );
        let eprop = eprop.try_borrow().expect("Cannot borrow edge");
        assert_eq!(eprop[e], eprop[hnew.edge()]);
        let hprop = hprop.try_borrow().expect("Cannot borrow halfedge");
        assert_eq!(hprop[h], hprop[hnew]);
        assert_eq!(hprop[oh], hprop[ohnew]);
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
        qbox.swap_edge_ccw(e).expect("Cannot swap edge");
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
        qbox.try_remove_edge(
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
        qbox.try_remove_edge(
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
                None,
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
                None,
            )
            .expect("Cannot insert edge");
        let (h, oh) = e.halfedges();
        mesh.check().expect("Topological errors found");
        assert!(!h.is_boundary(&mesh));
        assert!(!oh.is_boundary(&mesh));
        assert_eq!(3, mesh.loop_ccw_iter(h).count());
        assert_eq!(3, mesh.loop_ccw_iter(oh).count());
    }

    #[test]
    fn t_open_box_insert_edge_reuse_face() {
        let mut mesh = PolyMeshF32::unit_box().unwrap();
        mesh.delete_face(5.into(), false).unwrap();
        {
            let fstatus = mesh.face_status_prop();
            let fstatus = fstatus.try_borrow().unwrap();
            assert_eq!(5, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            assert_eq!(6, mesh.num_faces());
        }
        let h = mesh.halfedges().find(|h| h.is_boundary(&mesh)).unwrap();
        let e = mesh.insert_edge(h, h.prev(&mesh), Some(5.into())).unwrap();
        let (h, oh) = e.halfedges();
        assert!(!h.is_boundary(&mesh));
        assert!(oh.is_boundary(&mesh));
        assert_eq!(3, mesh.loop_ccw_iter(oh).count());
        assert_eq!(3, mesh.loop_ccw_iter(h).count());
        {
            let fstatus = mesh.face_status_prop();
            let fstatus = fstatus.try_borrow().unwrap();
            assert_eq!(6, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            assert_eq!(6, mesh.num_faces());
        }
    }

    #[test]
    fn t_box_insert_edge_reuse_face() {
        let mut mesh = PolyMeshF32::unit_box().unwrap();
        mesh.delete_face(5.into(), false).unwrap();
        let fstatus = mesh.face_status_prop();
        {
            let fstatus = fstatus.try_borrow().unwrap();
            assert_eq!(5, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            assert_eq!(6, mesh.num_faces());
        }
        mesh.add_quad_face(4.into(), 5.into(), 6.into(), 7.into())
            .unwrap();
        {
            let fstatus = fstatus.try_borrow().unwrap();
            assert_eq!(6, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
        }
        assert_eq!(7, mesh.num_faces());
        let h = mesh.find_halfedge(4.into(), 5.into()).unwrap();
        let e = mesh.insert_edge(h, h.prev(&mesh), Some(5.into())).unwrap();
        let (h, oh) = e.halfedges();
        assert!(!h.is_boundary(&mesh));
        assert!(!oh.is_boundary(&mesh));
        assert_eq!(3, mesh.loop_ccw_iter(oh).count());
        assert_eq!(3, mesh.loop_ccw_iter(h).count());
        {
            let fstatus = mesh.face_status_prop();
            let fstatus = fstatus.try_borrow().unwrap();
            assert_eq!(7, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            assert_eq!(7, mesh.num_faces());
        }
    }

    #[test]
    fn t_box_split_face() {
        let mut mesh = quad_box();
        const DEFAULT_PROP: u8 = 42u8;
        let mut fprop = mesh.create_face_prop(DEFAULT_PROP);
        {
            // Set face indices as property values.
            let mut fprop = fprop.try_borrow_mut().expect("Cannot borrow property");
            for f in mesh.faces() {
                fprop[f] = f.index() as u8;
            }
        }
        let vnew = mesh.add_vertex().expect("Cannot add vertex");
        mesh.split_face(5.into(), vnew, false)
            .expect("Cannot split face");
        let mesh = &mesh; // Immutable.
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(
            (4, 5),
            mesh.faces()
                .fold((0usize, 0usize), |(tris, quads), f| match f.valence(mesh) {
                    3 => (tris + 1, quads),
                    4 => (tris, quads + 1),
                    _ => (tris, quads),
                })
        );
        // Properties should not be copied to the new faces. They should get the default value.
        let fprop = fprop.try_borrow().expect("Cannot borrow property");
        assert_eq!(
            (6, 3),
            mesh.faces().fold((0usize, 0usize), |(old, new), f| {
                if fprop[f] == DEFAULT_PROP {
                    (old, 1 + new)
                } else {
                    (1 + old, new)
                }
            })
        );
    }

    #[test]
    fn t_box_split_face_copy_props() {
        let mut mesh = quad_box();
        const DEFAULT_PROP: u8 = 42u8;
        let mut fprop = mesh.create_face_prop(DEFAULT_PROP);
        {
            // Set face indices as property values.
            let mut fprop = fprop.try_borrow_mut().expect("Cannot borrow property");
            for f in mesh.faces() {
                fprop[f] = f.index() as u8;
            }
        }
        let vnew = mesh.add_vertex().expect("Cannot add vertex");
        mesh.split_face(5.into(), vnew, true)
            .expect("Cannot split face");
        let mesh = &mesh; // Immutable.
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(
            (4, 5),
            mesh.faces()
                .fold((0usize, 0usize), |(tris, quads), f| match f.valence(mesh) {
                    3 => (tris + 1, quads),
                    4 => (tris, quads + 1),
                    _ => (tris, quads),
                })
        );
        // Properties should not be copied to the new faces. They should get the default value.
        let fprop = fprop.try_borrow().expect("Cannot borrow property");
        assert_eq!(
            (9, 0),
            mesh.faces().fold((0usize, 0usize), |(old, new), f| {
                if fprop[f] == DEFAULT_PROP {
                    (old, 1 + new)
                } else {
                    (1 + old, new)
                }
            })
        );
    }

    #[test]
    fn t_box_insert_collapse() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot make a box");
        let h: HH = 6.into();
        let (hnew, _) = mesh
            .insert_edge(h.next(&mesh), h, None)
            .unwrap()
            .halfedges();
        mesh.check_topology().unwrap();
        assert!(mesh.try_check_edge_collapse(hnew).unwrap());
        let mut hcache = Vec::new();
        mesh.try_collapse_edge(hnew, &mut hcache).unwrap();
        mesh.check_topology().unwrap();
        mesh.garbage_collection().unwrap();
        assert_eq!(5, mesh.num_faces());
        assert_eq!(7, mesh.num_vertices());
        assert_eq!(10, mesh.num_edges());
    }
}
