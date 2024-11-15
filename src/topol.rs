use crate::{
    element::{Edge, Face, Halfedge, Handle, Vertex, EH, FH, HH, VH},
    error::Error,
    iterator,
    property::{EProperty, FProperty, HProperty, PropertyContainer, VProperty},
    status::Status,
};
use std::{cell::RefMut, ops::Range};

enum TentativeEdge {
    Old(HH),
    New {
        index: u32,
        from: VH,
        to: VH,
        prev: Option<HH>,
        next: Option<HH>,
        opp_prev: Option<HH>,
        opp_next: Option<HH>,
    },
}

#[derive(Default)]
pub struct TopolCache {
    loop_halfedges: Vec<Option<HH>>,
    needs_adjust: Vec<bool>,
    next_cache: Vec<(HH, HH)>,
    tentative: Vec<TentativeEdge>,
    pub(crate) halfedges: Vec<HH>,
    pub(crate) faces: Vec<FH>,
    pub(crate) edges: Vec<EH>,
    pub(crate) vertices: Vec<VH>,
}

impl TopolCache {
    fn clear(&mut self) {
        self.loop_halfedges.clear();
        self.needs_adjust.clear();
        self.next_cache.clear();
        self.tentative.clear();
        self.halfedges.clear();
        self.faces.clear();
        self.edges.clear();
        self.vertices.clear();
    }
}

pub struct Topology {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    pub(crate) vstatus: VProperty<Status>,
    pub(crate) hstatus: HProperty<Status>,
    pub(crate) estatus: EProperty<Status>,
    pub(crate) fstatus: FProperty<Status>,
    pub(crate) vprops: PropertyContainer<VH>,
    pub(crate) hprops: PropertyContainer<HH>,
    pub(crate) eprops: PropertyContainer<EH>,
    pub(crate) fprops: PropertyContainer<FH>,
}

impl Topology {
    pub fn new() -> Self {
        let (mut vprops, mut hprops, mut eprops, mut fprops) = (
            PropertyContainer::default(),
            PropertyContainer::default(),
            PropertyContainer::default(),
            PropertyContainer::default(),
        );
        let (vstatus, hstatus, estatus, fstatus) = (
            VProperty::new(&mut vprops, Default::default()),
            HProperty::new(&mut hprops, Default::default()),
            EProperty::new(&mut eprops, Default::default()),
            FProperty::new(&mut fprops, Default::default()),
        );
        Topology {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            vstatus,
            hstatus,
            estatus,
            fstatus,
            vprops,
            hprops,
            eprops,
            fprops,
        }
    }

    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        let (mut vprops, mut hprops, mut eprops, mut fprops) = (
            PropertyContainer::default(),
            PropertyContainer::default(),
            PropertyContainer::default(),
            PropertyContainer::default(),
        );
        let (vstatus, hstatus, estatus, fstatus) = (
            VProperty::with_capacity(nverts, &mut vprops, Default::default()),
            HProperty::with_capacity(nedges * 2, &mut hprops, Default::default()),
            EProperty::with_capacity(nedges, &mut eprops, Default::default()),
            FProperty::with_capacity(nfaces, &mut fprops, Default::default()),
        );
        Topology {
            vertices: Vec::with_capacity(nverts),
            edges: Vec::with_capacity(nedges),
            faces: Vec::with_capacity(nfaces),
            vstatus,
            hstatus,
            estatus,
            fstatus,
            vprops,
            hprops,
            eprops,
            fprops,
        }
    }

    pub fn reserve(&mut self, nverts: usize, nedges: usize, nfaces: usize) -> Result<(), Error> {
        // Elements.
        self.vertices.reserve(nverts);
        self.edges.reserve(nedges);
        self.faces.reserve(nfaces);
        // Properties.
        self.vprops.reserve(nverts)?;
        self.hprops.reserve(nedges * 2)?;
        self.eprops.reserve(nedges)?;
        self.fprops.reserve(nfaces)?;
        Ok(())
    }

    pub fn clear(&mut self) -> Result<(), Error> {
        // Elements.
        self.vertices.clear();
        self.edges.clear();
        self.faces.clear();
        // Properties.
        self.vprops.clear()?;
        self.hprops.clear()?;
        self.eprops.clear()?;
        self.fprops.clear()?;
        Ok(())
    }

    pub fn vertex_status(&self, v: VH) -> Result<Status, Error> {
        self.vstatus.get(v)
    }

    pub fn vertex_status_mut(&mut self, v: VH) -> Result<RefMut<'_, Status>, Error> {
        self.vstatus.get_mut(v)
    }

    pub fn halfedge_status(&self, h: HH) -> Result<Status, Error> {
        self.hstatus.get(h)
    }

    pub fn edge_status(&self, e: EH) -> Result<Status, Error> {
        self.estatus.get(e)
    }

    pub fn edge_status_mut(&mut self, e: EH) -> Result<RefMut<'_, Status>, Error> {
        self.estatus.get_mut(e)
    }

    pub fn face_status(&self, f: FH) -> Result<Status, Error> {
        self.fstatus.get(f)
    }

    pub fn face_status_mut(&mut self, f: FH) -> Result<RefMut<'_, Status>, Error> {
        self.fstatus.get_mut(f)
    }

    pub fn create_vertex_prop<T>(&mut self, default: T) -> VProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        VProperty::<T>::new(&mut self.vprops, default)
    }

    pub fn create_halfedge_prop<T>(&mut self, default: T) -> HProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        HProperty::<T>::new(&mut self.hprops, default)
    }

    pub fn create_edge_prop<T>(&mut self, default: T) -> EProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        EProperty::<T>::new(&mut self.eprops, default)
    }

    pub fn create_face_prop<T>(&mut self, default: T) -> FProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        FProperty::<T>::new(&mut self.fprops, default)
    }

    pub fn new_vprop_with_capacity<T>(&mut self, n: usize, default: T) -> VProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        VProperty::<T>::with_capacity(n, &mut self.vprops, default)
    }

    pub fn is_valid_vertex(&self, v: VH) -> bool {
        (v.index() as usize) < self.num_vertices()
    }

    pub fn is_valid_halfedge(&self, h: HH) -> bool {
        (h.index() as usize) < self.num_halfedges()
    }

    pub fn is_valid_edge(&self, e: EH) -> bool {
        (e.index() as usize) < self.num_edges()
    }

    pub fn is_valid_face(&self, f: FH) -> bool {
        (f.index() as usize) < self.num_faces()
    }

    fn vertex(&self, v: VH) -> &Vertex {
        &self.vertices[v.index() as usize]
    }

    pub(crate) fn vertex_mut(&mut self, v: VH) -> &mut Vertex {
        &mut self.vertices[v.index() as usize]
    }

    pub(crate) fn halfedge(&self, h: HH) -> &Halfedge {
        &self.edges[(h.index() >> 1) as usize].halfedges[(h.index() & 1) as usize]
    }

    pub(crate) fn halfedge_mut(&mut self, h: HH) -> &mut Halfedge {
        &mut self.edges[(h.index() >> 1) as usize].halfedges[(h.index() & 1) as usize]
    }

    pub(crate) fn face_mut(&mut self, f: FH) -> &mut Face {
        &mut self.faces[f.index() as usize]
    }

    pub fn vertex_halfedge(&self, v: VH) -> Option<HH> {
        self.vertex(v).halfedge
    }

    pub fn head_vertex(&self, h: HH) -> VH {
        self.halfedge(h).vertex
    }

    pub fn tail_vertex(&self, h: HH) -> VH {
        self.halfedge(self.opposite_halfedge(h)).vertex
    }

    pub fn prev_halfedge(&self, h: HH) -> HH {
        self.halfedge(h).prev
    }

    pub fn next_halfedge(&self, h: HH) -> HH {
        self.halfedge(h).next
    }

    pub fn halfedge_face(&self, h: HH) -> Option<FH> {
        self.halfedge(h).face
    }

    pub fn halfedge_edge(&self, h: HH) -> EH {
        (h.index() >> 1).into()
    }

    pub fn halfedge_pair(&self, e: EH) -> (HH, HH) {
        e.halfedges()
    }

    pub fn edge_halfedge(&self, e: EH, flag: bool) -> HH {
        e.halfedge(flag)
    }

    pub fn face_halfedge(&self, f: FH) -> HH {
        self.faces[f.index() as usize].halfedge
    }

    pub fn is_boundary_halfedge(&self, h: HH) -> bool {
        self.halfedge(h).face.is_none()
    }

    pub fn is_boundary_edge(&self, e: EH) -> bool {
        e.is_boundary(self)
    }

    pub fn is_boundary_vertex(&self, v: VH) -> bool {
        match self.vertices[v.index() as usize].halfedge {
            Some(h) => self.is_boundary_halfedge(h),
            None => true,
        }
    }

    pub fn opposite_halfedge(&self, h: HH) -> HH {
        (h.index() ^ 1).into()
    }

    pub fn cw_rotated_halfedge(&self, h: HH) -> HH {
        self.halfedge(self.opposite_halfedge(h)).next
    }

    pub fn ccw_rotated_halfedge(&self, h: HH) -> HH {
        self.opposite_halfedge(self.halfedge(h).prev)
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    pub fn num_halfedges(&self) -> usize {
        self.num_edges() * 2
    }

    pub fn vertices(&self) -> impl Iterator<Item = VH> {
        (0..(self.num_vertices() as u32)).map(|i| i.into())
    }

    pub fn halfedges(&self) -> impl Iterator<Item = HH> {
        (0..(self.num_halfedges() as u32)).map(|i| i.into())
    }

    pub fn edges(&self) -> impl Iterator<Item = EH> {
        (0..(self.num_edges() as u32)).map(|i| i.into())
    }

    pub fn faces(&self) -> impl Iterator<Item = FH> {
        (0..(self.num_faces() as u32)).map(|i| i.into())
    }

    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    pub fn find_halfedge(&self, from: VH, to: VH) -> Option<HH> {
        iterator::voh_ccw_iter(self, from).find(|h| self.head_vertex(*h) == to)
    }

    pub fn is_manifold_vertex(&self, v: VH) -> bool {
        /* If just the first outgoing halfedge is on the boundary, it just means
         * the vertex is on the boundary. If the first outgoing halfedge is not
         * on the boundary, it implies the vertex is in the interior. In both
         * cases the vertex is manifold. If any outgoing halfedge apart from the
         * first is on the boundary, it implies there are more than one gaps
         * when circulating around the vertex, making it non-manifold. For this
         * reason, we skip the first halfedge and check the rest. */
        iterator::voh_ccw_iter(self, v)
            .skip(1)
            .all(|h| !self.is_boundary_halfedge(h))
    }

    pub(crate) fn adjust_outgoing_halfedge(&mut self, v: VH) {
        let h = iterator::voh_ccw_iter(self, v).find(|h| self.is_boundary_halfedge(*h));
        if let Some(h) = h {
            self.vertex_mut(v).halfedge = Some(h);
        }
    }

    pub fn add_vertex(&mut self) -> Result<VH, Error> {
        let vi = self.vertices.len() as u32;
        self.vprops.push_value()?;
        self.vertices.push(Vertex { halfedge: None });
        Ok(vi.into())
    }

    pub fn add_vertices(&mut self, n: usize) -> Result<Range<u32>, Error> {
        let nverts = self.vertices.len() as u32;
        self.vprops.push_values(n)?;
        self.vertices
            .resize(self.vertices.len() + n, Vertex { halfedge: None });
        Ok(nverts..(nverts + n as u32))
    }

    pub(crate) fn new_edge(
        &mut self,
        from: VH,
        to: VH,
        prev: HH,
        next: HH,
        opp_prev: HH,
        opp_next: HH,
    ) -> Result<EH, Error> {
        let ei = self.edges.len() as u32;
        self.eprops.push_value()?;
        self.hprops.push_values(2)?;
        self.edges.push(Edge {
            halfedges: [
                Halfedge {
                    face: None,
                    vertex: to,
                    next,
                    prev,
                },
                Halfedge {
                    face: None,
                    vertex: from,
                    next: opp_next,
                    prev: opp_prev,
                },
            ],
        });
        Ok(ei.into())
    }

    pub(crate) fn new_face(&mut self, halfedge: HH) -> Result<FH, Error> {
        let fi = self.faces.len() as u32;
        self.fprops.push_value()?;
        self.faces.push(Face { halfedge });
        Ok(fi.into())
    }

    pub fn link_halfedges(&mut self, hprev: HH, hnext: HH) {
        self.halfedge_mut(hprev).next = hnext;
        self.halfedge_mut(hnext).prev = hprev;
    }

    pub fn add_face(&mut self, verts: &[VH], cache: &mut TopolCache) -> Result<FH, Error> {
        cache.clear();
        cache.loop_halfedges.reserve(verts.len());
        cache.needs_adjust.reserve(verts.len());
        cache.next_cache.reserve(verts.len() * 6);
        // Check for topological errors.
        for i in 0..verts.len() {
            if !self.is_boundary_vertex(verts[i]) {
                // Ensure vertex is manifold.
                return Err(Error::ComplexVertex(verts[i]));
            }
            // Ensure edge is manifold.
            let h = self.find_halfedge(verts[i], verts[(i + 1) % verts.len()]);
            match h {
                Some(h) if !self.is_boundary_halfedge(h) => return Err(Error::ComplexHalfedge(h)),
                _ => {} // Do nothing.
            }
            cache.loop_halfedges.push(h);
            cache.needs_adjust.push(false);
        }
        // If any vertex has more than two incident boundary edges, relinking might be necessary.
        for (prev, next) in (0..verts.len()).filter_map(|i| {
            match (
                cache.loop_halfedges[i],
                cache.loop_halfedges[(i + 1) % verts.len()],
            ) {
                (Some(prev), Some(next)) if self.next_halfedge(prev) != next => Some((prev, next)),
                _ => None,
            }
        }) {
            // Relink the patch.
            let boundprev = {
                let mut out = self.opposite_halfedge(next);
                loop {
                    out = self.opposite_halfedge(self.next_halfedge(out));
                    if self.is_boundary_halfedge(out) {
                        break;
                    }
                }
                out
            };
            let boundnext = self.next_halfedge(boundprev);
            // Ok ?
            if boundprev == prev {
                return Err(Error::PatchRelinkingFailed);
            }
            debug_assert!(
                self.is_boundary_halfedge(boundprev) && self.is_boundary_halfedge(boundnext)
            );
            // other halfedges.
            let pstart = self.next_halfedge(prev);
            let pend = self.prev_halfedge(next);
            // relink.
            cache.next_cache.extend_from_slice(&[
                (boundprev, pstart),
                (pend, boundnext),
                (prev, next),
            ]);
        }
        // Create boundary loop. No more errors allowed from this point.
        // If anything goes wrong, we panic.
        cache.tentative.clear();
        cache.tentative.reserve(verts.len());
        {
            let mut ei = self.edges.len() as u32;
            cache
                .tentative
                .extend((0..verts.len()).map(|i| match cache.loop_halfedges[i] {
                    Some(h) => TentativeEdge::Old(h),
                    None => TentativeEdge::New {
                        index: {
                            let current = ei;
                            ei += 1;
                            current << 1
                        },
                        from: verts[i],
                        to: verts[(i + 1) % verts.len()],
                        prev: None,
                        next: None,
                        opp_prev: None,
                        opp_next: None,
                    },
                }));
        }
        for (i, j) in (0..verts.len()).map(|i| (i, (i + 1) % verts.len())) {
            let (e0, e1) = if j == 0 {
                let (right, left) = cache.tentative.split_at_mut(i);
                (&mut left[0], &mut right[0])
            } else {
                let (left, right) = cache.tentative.split_at_mut(j);
                (&mut left[left.len() - 1], &mut right[0])
            };
            let v = verts[j];
            match (e0, e1) {
                (TentativeEdge::Old(_), TentativeEdge::Old(innernext)) => {
                    cache.needs_adjust[j] = self.vertex_halfedge(v) == Some(*innernext);
                }
                (
                    TentativeEdge::New {
                        index: innerprev,
                        opp_prev,
                        next,
                        ..
                    },
                    TentativeEdge::Old(innernext),
                ) => {
                    let innernext = *innernext;
                    let innerprev = innerprev.into();
                    let outernext = self.opposite_halfedge(innerprev);
                    let boundprev = self.prev_halfedge(innernext);
                    cache.next_cache.push((boundprev, outernext));
                    *opp_prev = Some(boundprev);
                    cache.next_cache.push((innerprev, innernext));
                    *next = Some(innernext);
                    self.vertex_mut(v).halfedge = Some(outernext);
                }
                (
                    TentativeEdge::Old(innerprev),
                    TentativeEdge::New {
                        index: innernext,
                        prev,
                        opp_next,
                        ..
                    },
                ) => {
                    let innerprev = *innerprev;
                    let innernext = innernext.into();
                    let outerprev = self.opposite_halfedge(innernext);
                    let boundnext = self.next_halfedge(innerprev);
                    cache.next_cache.push((outerprev, boundnext));
                    *opp_next = Some(boundnext);
                    cache.next_cache.push((innerprev, innernext));
                    *prev = Some(innerprev);
                    self.vertex_mut(v).halfedge = Some(boundnext);
                }
                (
                    TentativeEdge::New {
                        index: innerprev,
                        next,
                        opp_prev,
                        ..
                    },
                    TentativeEdge::New {
                        index: innernext,
                        prev,
                        opp_next,
                        ..
                    },
                ) => {
                    let innerprev = innerprev.into();
                    let innernext = innernext.into();
                    let outernext = self.opposite_halfedge(innerprev);
                    let outerprev = self.opposite_halfedge(innernext);
                    if let Some(boundnext) = self.vertex_halfedge(v) {
                        let boundprev = self.prev_halfedge(boundnext);
                        cache
                            .next_cache
                            .extend(&[(boundprev, outernext), (outerprev, boundnext)]);
                        *next = Some(innernext);
                        *opp_prev = Some(boundprev);
                        *prev = Some(innerprev);
                        *opp_next = Some(boundnext);
                    } else {
                        self.vertex_mut(v).halfedge = Some(outernext);
                        *next = Some(innernext);
                        *opp_prev = Some(outerprev);
                        *prev = Some(innerprev);
                        *opp_next = Some(outernext);
                    }
                }
            };
        }
        // Convert tentative edges into real edges.
        cache.halfedges.reserve(cache.tentative.len());
        const ERR: &str = "Unable to create edge loop";
        for tedge in &cache.tentative {
            let h = match tedge {
                TentativeEdge::Old(h) => *h,
                TentativeEdge::New {
                    index,
                    from,
                    to,
                    prev,
                    next,
                    opp_prev,
                    opp_next,
                } => {
                    let ei = self.new_edge(
                        *from,
                        *to,
                        prev.expect(ERR),
                        next.expect(ERR),
                        opp_prev.expect(ERR),
                        opp_next.expect(ERR),
                    )?;
                    assert_eq!(*index >> 1, ei.index(), "Failed to create an edge loop");
                    index.into()
                }
            };
            cache.halfedges.push(h);
        }
        // Create the face.
        let fnew = self.new_face(match cache.tentative.last().expect(ERR) {
            TentativeEdge::Old(index) => *index,
            TentativeEdge::New { index, .. } => index.into(),
        })?;
        for h in &cache.halfedges {
            self.halfedge_mut(*h).face = Some(fnew);
        }
        // Process next halfedge cache.
        for (prev, next) in cache.next_cache.drain(..) {
            self.link_halfedges(prev, next);
        }
        // Adjust vertices' halfedge handles.
        for (i, vert) in verts.iter().enumerate() {
            if cache.needs_adjust[i] {
                self.adjust_outgoing_halfedge(*vert);
            }
        }
        Ok(fnew)
    }

    pub fn delete_vertex(
        &mut self,
        v: VH,
        delete_isolated_vertices: bool,
        cache: &mut TopolCache,
    ) -> Result<(), Error> {
        /* Deleting faces changes the local topology. So we cannot delete them
         * as we iterate over them. Instead we collect them into cache and then
         * delete them. */
        cache.faces.clear();
        cache.faces.extend(iterator::vf_ccw_iter(self, v));
        for f in &cache.faces {
            self.delete_face(
                *f,
                delete_isolated_vertices,
                &mut cache.edges,
                &mut cache.vertices,
            )?;
        }
        self.vertex_status_mut(v)?.set_deleted(true);
        Ok(())
    }

    pub fn delete_edge(
        &mut self,
        e: EH,
        delete_isolated_vertices: bool,
        ecache: &mut Vec<EH>,
        vcache: &mut Vec<VH>,
    ) -> Result<(), Error> {
        let h0 = self.edge_halfedge(e, false);
        let h1 = self.edge_halfedge(e, true);
        let f0 = self.halfedge_face(h0);
        let f1 = self.halfedge_face(h1);
        if let Some(f) = f0 {
            self.delete_face(f, delete_isolated_vertices, ecache, vcache)?;
        }
        if let Some(f) = f1 {
            self.delete_face(f, delete_isolated_vertices, ecache, vcache)?;
        }
        /* If either face was valid, the edge is deleted inside the call to
         * delete_face. Otherwise we have to mark them deleted here. */
        if f0.is_none() && f1.is_none() {
            self.edge_status_mut(e)?.set_deleted(true);
            {
                let mut hstatus = self.hstatus.try_borrow_mut()?;
                hstatus[h0.index() as usize].set_deleted(true);
                hstatus[h1.index() as usize].set_deleted(true);
            }
        }
        Ok(())
    }

    pub fn delete_face(
        &mut self,
        f: FH,
        delete_isolated_vertices: bool,
        ecache: &mut Vec<EH>,
        vcache: &mut Vec<VH>,
    ) -> Result<(), Error> {
        {
            let mut fs = self.face_status_mut(f)?;
            // Check if the face was already deleted.
            if fs.deleted() {
                return Err(Error::DeletedFace(f));
            }
            fs.set_deleted(true); // Mark face deleted.
        }
        // Collect neighborhood topology.
        ecache.clear();
        vcache.clear();
        for (mesh, h) in iterator::fh_ccw_iter_mut(self, f) {
            mesh.halfedge_mut(h).face = None; // Disconnect from face.
            if mesh.is_boundary_halfedge(mesh.opposite_halfedge(h)) {
                ecache.push(mesh.halfedge_edge(h));
            }
            vcache.push(mesh.head_vertex(h));
        }
        // Delete collected topology.
        for e in ecache.drain(..) {
            let h0 = self.edge_halfedge(e, false);
            let v0 = self.head_vertex(h0);
            let next0 = self.next_halfedge(h0);
            let prev0 = self.prev_halfedge(h0);
            let h1 = self.edge_halfedge(e, true);
            let v1 = self.head_vertex(h1);
            let next1 = self.next_halfedge(h1);
            let prev1 = self.prev_halfedge(h1);
            // Adjust halfedge links and mark edge and halfedges deleted.
            self.link_halfedges(prev0, next1);
            self.link_halfedges(prev1, next0);
            self.edge_status_mut(e)?.set_deleted(true);
            {
                let mut hstatus = self.hstatus.try_borrow_mut()?;
                hstatus[h0.index() as usize].set_deleted(true);
                hstatus[h1.index() as usize].set_deleted(true);
            }
            // Update vertices.
            if self.vertex_halfedge(v0) == Some(h1) {
                // Isolated?
                if next0 == h1 {
                    if delete_isolated_vertices {
                        self.vertex_status_mut(v0)?.set_deleted(true);
                    }
                    self.vertex_mut(v0).halfedge = None;
                } else {
                    self.vertex_mut(v0).halfedge = Some(next0);
                }
            }
            if self.vertex_halfedge(v1) == Some(h0) {
                // Isolated?
                if next1 == h0 {
                    if delete_isolated_vertices {
                        self.vertex_status_mut(v1)?.set_deleted(true);
                    }
                    self.vertex_mut(v1).halfedge = None;
                } else {
                    self.vertex_mut(v1).halfedge = Some(next1);
                }
            }
        }
        // Update outgoing halfedges.
        for v in vcache.drain(..) {
            self.adjust_outgoing_halfedge(v);
        }
        Ok(())
    }

    pub fn garbage_collection(&mut self, cache: &mut TopolCache) -> Result<(), Error> {
        // Setup mapping.
        let vmap = &mut cache.vertices;
        vmap.clear();
        vmap.reserve(self.num_vertices());
        vmap.extend(self.vertices());
        let hmap = &mut cache.halfedges;
        hmap.clear();
        hmap.reserve(self.num_halfedges());
        hmap.extend(self.halfedges());
        let fmap = &mut cache.faces;
        fmap.clear();
        fmap.reserve(self.num_faces());
        fmap.extend(self.faces());
        // Remove vertices.
        if self.num_vertices() > 0 {
            let mut left = 0usize;
            let mut right = self.num_vertices() - 1;
            let newlen = loop {
                // Find first deleted and last un-deleted.
                // Use scope to borrow the status vector.
                {
                    let status = self.vstatus.try_borrow()?;
                    while !status[left].deleted() && left < right {
                        left += 1;
                    }
                    while status[right].deleted() && left < right {
                        right -= 1;
                    }
                    if left >= right {
                        break left + if status[left].deleted() { 0 } else { 1 };
                    }
                }
                // Swap.
                self.vertices.swap(left, right);
                vmap.swap(left, right);
                self.vprops.swap(left, right)?;
            };
            self.vertices.truncate(newlen);
            self.vprops.resize(self.num_vertices())?;
        }
        // Remove edges.
        if self.num_edges() > 0 {
            let mut left = 0usize;
            let mut right = self.num_edges() - 1;
            let newlen = loop {
                // Find first deleted and last un-deleted.
                // Use scope to borrow the status vector.
                {
                    let status = self.estatus.try_borrow()?;
                    while !status[left].deleted() && left < right {
                        left += 1;
                    }
                    while status[right].deleted() && left < right {
                        right -= 1;
                    }
                    if left >= right {
                        break left + if status[left].deleted() { 0 } else { 1 };
                    }
                }
                // Swap.
                self.edges.swap(left, right);
                hmap.swap(2 * left, 2 * right);
                hmap.swap(2 * left + 1, 2 * right + 1);
                self.eprops.swap(left, right)?;
                self.hprops.swap(2 * left, 2 * right)?;
                self.hprops.swap(2 * left + 1, 2 * right + 1)?;
            };
            self.edges.truncate(newlen);
            self.eprops.resize(self.num_edges())?;
            self.hprops.resize(self.num_halfedges())?;
        }
        // Remove faces.
        if self.num_faces() > 0 {
            let mut left = 0usize;
            let mut right = self.num_faces() - 1;
            let newlen = loop {
                // Find first deleted and last un-deleted.
                // Use scope to borrow the status vector.
                {
                    let status = self.fstatus.try_borrow()?;
                    while !status[left].deleted() && left < right {
                        left += 1;
                    }
                    while status[right].deleted() && left < right {
                        right -= 1;
                    }
                    if left >= right {
                        break left + if status[left].deleted() { 0 } else { 1 };
                    }
                }
                // Swap.
                self.faces.swap(left, right);
                fmap.swap(left, right);
                self.fprops.swap(left, right)?;
            };
            self.faces.truncate(newlen);
            self.fprops.resize(self.num_faces())?;
        }
        // Update handles of vertices.
        for v in self.vertices.iter_mut() {
            if let Some(h) = v.halfedge {
                v.halfedge = Some(hmap[h.index() as usize]);
            }
        }
        // Update edges vertices.
        for e in self.edges.iter_mut() {
            for h in e.halfedges.iter_mut() {
                h.vertex = vmap[h.vertex.index() as usize];
            }
        }
        // Update halfedge connectivity.
        for h in self.halfedges() {
            self.link_halfedges(h, hmap[self.next_halfedge(h).index() as usize]);
            if let Some(f) = self.halfedge_face(h) {
                self.halfedge_mut(h).face = Some(fmap[f.index() as usize]);
            }
        }
        // Update handles of faces.
        for f in self.faces.iter_mut() {
            f.halfedge = hmap[f.halfedge.index() as usize];
        }
        // Clean up properties.
        self.vprops.garbage_collection();
        self.hprops.garbage_collection();
        self.eprops.garbage_collection();
        self.fprops.garbage_collection();
        Ok(())
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub(crate) mod test {
    use arrayvec::ArrayVec;

    use crate::{
        alum_glam::PolyMeshF32,
        iterator,
        macros::assert_f32_eq,
        topol::{Handle, VH},
    };

    use super::{TopolCache, Topology};

    /**
     * Makes a box with the following topology.
     * ```text
     *
     *      7-----------6
     *     /|          /|
     *    / |         / |
     *   4-----------5  |
     *   |  |        |  |
     *   |  3--------|--2
     *   | /         | /
     *   |/          |/
     *   0-----------1
     * ```
     */
    pub fn quad_box() -> Topology {
        let mut topol = Topology::with_capacity(8, 12, 6);
        let verts: Vec<_> = (0..8)
            .map(|_| topol.add_vertex().expect("Unable to add a vertex").index())
            .collect();
        assert_eq!(verts, (0u32..8).collect::<Vec<_>>());
        let mut cache = TopolCache::default();
        let faces: Vec<_> = [
            [0u32, 3, 2, 1],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
            [4, 5, 6, 7],
        ]
        .iter()
        .map(|indices| {
            topol
                .add_face(&indices.map(|i| i.into()), &mut cache)
                .expect("Unable to add a face")
        })
        .collect();
        assert_eq!(faces, (0..6).map(|i| i.into()).collect::<Vec<_>>());
        assert_eq!(topol.num_vertices(), 8);
        assert_eq!(topol.num_halfedges(), 24);
        assert_eq!(topol.num_edges(), 12);
        assert_eq!(topol.num_faces(), 6);
        topol
    }

    /**

                        12---------13---------14---------15
                       /          /          /          /
                      /   f5     /   f6     /    f7    /
                     /          /          /          /
                    /          /          /          /
                   8----------9----------10---------11
                  /          /          /          /
                 /    f3    /          /    f4    /
                /          /          /          /
               /          /          /          /
              4----------5----------6----------7
             /          /          /          /
            /   f0     /    f1    /    f2    /
           /          /          /          /
          /          /          /          /
         0----------1----------2----------3
    */
    pub fn loop_mesh() -> Topology {
        let mut topol = Topology::with_capacity(16, 24, 8);
        let mut cache = TopolCache::default();
        for _ in 0u32..16 {
            topol.add_vertex().expect("Unable to add vertex");
        }
        for fvi in [
            [0u32, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [4, 5, 9, 8],
            [6, 7, 11, 10],
            [8, 9, 13, 12],
            [9, 10, 14, 13],
            [10, 11, 15, 14],
        ] {
            let vs = fvi.iter().map(|i| i.into()).collect::<ArrayVec<VH, 4>>();
            topol.add_face(&vs, &mut cache).expect("Unable to add face");
        }
        topol
    }

    #[test]
    fn t_triangle() {
        let mut topol = Topology::default();
        let mut cache = TopolCache::default();
        let verts: Vec<_> = (0..3).flat_map(|_| topol.add_vertex()).collect();
        assert_eq!(verts, (0..3u32).map(|idx| idx.into()).collect::<Vec<_>>());
        let face = topol.add_face(&verts, &mut cache).unwrap();
        assert_eq!(topol.num_faces(), 1);
        assert_eq!(topol.num_edges(), 3);
        assert_eq!(topol.num_halfedges(), 6);
        assert_eq!(topol.num_vertices(), 3);
        assert_eq!(face.index(), 0);
        for v in topol.vertices() {
            let h = topol
                .vertex_halfedge(v)
                .expect("Vertex must have an incident halfedge");
            assert!(topol.is_boundary_halfedge(h));
            let oh = topol.opposite_halfedge(h);
            assert!(!topol.is_boundary_halfedge(oh));
            assert_eq!(
                topol
                    .halfedge_face(oh)
                    .expect("Halfedge must have an incident face"),
                face
            );
        }
        assert_eq!(
            topol
                .halfedges()
                .filter(|h| topol.is_boundary_halfedge(*h))
                .count(),
            3
        );
        assert_eq!(
            topol
                .halfedges()
                .filter(|h| !topol.is_boundary_halfedge(*h))
                .count(),
            3
        );
        for (i, j) in (0u32..3).map(|i| (i, (i + 1) % 3)) {
            let h = topol.find_halfedge(i.into(), j.into()).unwrap();
            assert!(!topol.is_boundary_halfedge(h));
        }
    }

    #[test]
    fn t_two_triangles() {
        let mut topol = Topology::default();
        let mut cache = TopolCache::default();
        let verts: Vec<_> = (0..4)
            .map(|_| topol.add_vertex().expect("Cannot add vertex"))
            .collect();
        assert_eq!(verts.len(), 4);
        let faces = vec![
            topol
                .add_face(&[verts[0], verts[1], verts[2]], &mut cache)
                .expect("Cannot add face"),
            topol
                .add_face(&[verts[0], verts[2], verts[3]], &mut cache)
                .expect("Cannot add face"),
        ];
        assert_eq!(
            faces,
            [0u32, 1].iter().map(|idx| idx.into()).collect::<Vec<_>>()
        );
        assert_eq!(topol.num_vertices(), 4);
        assert_eq!(topol.num_halfedges(), 10);
        assert_eq!(topol.num_edges(), 5);
        assert_eq!(topol.num_faces(), 2);
        assert_eq!(
            topol.edges().filter(|e| topol.is_boundary_edge(*e)).count(),
            4
        );
        assert_eq!(
            topol
                .edges()
                .filter(|e| !topol.is_boundary_edge(*e))
                .count(),
            1
        );
    }

    #[test]
    fn t_quad() {
        let mut topol = Topology::default();
        let mut cache = TopolCache::default();
        let verts: Vec<_> = (0..4).flat_map(|_| topol.add_vertex()).collect();
        assert_eq!(verts, vec![0.into(), 1.into(), 2.into(), 3.into()]);
        let face = topol.add_face(&verts, &mut cache).unwrap();
        assert_eq!(topol.num_faces(), 1);
        assert_eq!(topol.num_edges(), 4);
        assert_eq!(topol.num_halfedges(), 8);
        assert_eq!(topol.num_vertices(), 4);
        assert_eq!(face, 0.into());
        for v in topol.vertices() {
            let h = topol
                .vertex_halfedge(v)
                .expect("Vertex must have an incident halfedge");
            assert!(topol.is_boundary_halfedge(h));
            let oh = topol.opposite_halfedge(h);
            assert!(!topol.is_boundary_halfedge(oh));
            assert_eq!(
                topol
                    .halfedge_face(oh)
                    .expect("Halfedge must have an incident face"),
                face
            );
        }
        assert_eq!(
            topol
                .halfedges()
                .filter(|h| topol.is_boundary_halfedge(*h))
                .count(),
            4
        );
        assert_eq!(
            topol
                .halfedges()
                .filter(|h| !topol.is_boundary_halfedge(*h))
                .count(),
            4
        );
    }

    #[test]
    fn t_box_manifold() {
        let qbox = quad_box();
        assert!(
            qbox.halfedges().all(|h| !qbox.is_boundary_halfedge(h)),
            "Not expecting any boundary edges"
        );
    }

    #[test]
    fn t_box_vertex_valences() {
        let qbox = quad_box();
        for v in qbox.vertices() {
            assert_eq!(v.valence(&qbox), 3);
        }
    }

    #[test]
    fn t_box_face_valence() {
        let qbox = quad_box();
        for f in qbox.faces() {
            assert_eq!(f.valence(&qbox), 4);
        }
    }

    #[test]
    fn t_loop_mesh_add_face() {
        /*

                            12---------13---------14---------15
                           /          /          /          /
                          /   f5     /   f6     /    f7    /
                         /          /          /          /
                        /          /          /          /
                       8----------9----------10---------11
                      /          / v1--v0   /          /
                     /    f3    /    \  |  /    f4    /
                    /          /      \ | /          /
                   /          /        \|/          /
                  4----------5----------6----------7
                 /          /          /          /
                /   f0     /    f1    /    f2    /
               /          /          /          /
              /          /          /          /
             0----------1----------2----------3
        */
        let mut mesh = loop_mesh();
        assert_eq!(
            mesh.edges().filter(|e| mesh.is_boundary_edge(*e)).count(),
            16
        );
        // Add a floating triangle at v6.
        let v0 = mesh.add_vertex().expect("Unable to add vertex");
        let v1 = mesh.add_vertex().expect("Unable to add vertex");
        let mut cache = TopolCache::default();
        let f0 = mesh
            .add_face(&[6.into(), v0, v1], &mut cache)
            .expect("Unable to add a face");
        assert_eq!(f0.index(), 8);
        assert_eq!(
            iterator::vf_ccw_iter(&mesh, 6.into())
                .map(|i| i.index())
                .collect::<Vec<_>>(),
            [8, 1, 2, 4]
        );
        assert_eq!(
            iterator::vv_ccw_iter(&mesh, 6.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            [10, v0.index(), v1.index(), 5, 2, 7]
        );
        assert_eq!(iterator::ve_ccw_iter(&mesh, 6.into()).count(), 6);
        assert_eq!(
            mesh.vertices()
                .filter(|v| mesh.is_manifold_vertex(*v))
                .count(),
            17
        );
        assert!(!mesh.is_manifold_vertex(6.into()));
        // Check halfedge connectivity at v6.
        {
            let h = mesh
                .find_halfedge(5.into(), 6.into())
                .expect("Cannot find halfedge");
            assert_eq!(mesh.tail_vertex(h), 5.into());
            assert_eq!(mesh.head_vertex(h), 6.into());
            let h1 = mesh.next_halfedge(h);
            assert_eq!(
                h1,
                mesh.find_halfedge(6.into(), v1)
                    .expect("Cannot find halfedge")
            );
            let h = mesh
                .find_halfedge(6.into(), 10.into())
                .expect("Cannot find halfedge");
            assert_eq!(mesh.tail_vertex(h), 6.into());
            assert_eq!(mesh.head_vertex(h), 10.into());
            let h = mesh.prev_halfedge(h);
            assert_eq!(
                h,
                mesh.find_halfedge(v0, 6.into())
                    .expect("Cannot find halfedge")
            );
        }
        assert_eq!(
            mesh.edges().filter(|e| mesh.is_boundary_edge(*e)).count(),
            19
        );
        // Add another face.
        let f1 = mesh
            .add_face(&[5.into(), 6.into(), v1], &mut cache)
            .expect("Unable to add face");
        assert_eq!(f1.index(), 9);
        assert_eq!(mesh.num_edges(), 28);
        assert_eq!(
            mesh.edges().filter(|e| mesh.is_boundary_edge(*e)).count(),
            18
        );
        assert_eq!(
            iterator::vf_ccw_iter(&mesh, 6.into())
                .map(|f| f.index())
                .collect::<Vec<_>>(),
            [8, 9, 1, 2, 4]
        );
        assert_eq!(
            iterator::vf_ccw_iter(&mesh, 5.into())
                .map(|f| f.index())
                .collect::<Vec<_>>(),
            [3, 0, 1, 9]
        );
        assert_eq!(
            iterator::vv_ccw_iter(&mesh, 6.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            [10, v0.index(), v1.index(), 5, 2, 7]
        );
        assert_eq!(
            iterator::vv_ccw_iter(&mesh, 5.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            [v1.index(), 9, 4, 1, 6]
        );
        assert!(mesh.is_manifold_vertex(6.into()));
        assert!(mesh.is_manifold_vertex(5.into()));
    }

    #[test]
    fn t_quad_box_prop_len() {
        let qbox = quad_box();
        assert_eq!(qbox.num_vertices(), qbox.vprops.len());
        assert_eq!(qbox.num_halfedges(), qbox.hprops.len());
        assert_eq!(qbox.num_edges(), qbox.eprops.len());
        assert_eq!(qbox.num_faces(), qbox.fprops.len());
    }

    #[test]
    fn t_quad_box_delete_face() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        qbox.delete_face(5.into(), true, &mut cache.edges, &mut cache.vertices)
            .expect("Cannot delete face");
        assert!(qbox
            .face_status(5.into())
            .expect("Cannot read face status")
            .deleted());
        qbox.garbage_collection(&mut cache)
            .expect("Cannot garbage collect mesh");
        assert_eq!(qbox.num_faces(), 5);
        assert_eq!(qbox.num_edges(), 12);
        assert_eq!(qbox.num_vertices(), 8);
        for v in qbox.vertices() {
            assert_eq!(qbox.is_boundary_vertex(v), v.index() > 3);
        }
        assert_eq!(
            4,
            qbox.edges().filter(|e| qbox.is_boundary_edge(*e)).count()
        );
        assert_eq!(
            4,
            qbox.vertices()
                .filter(|v| qbox.is_boundary_vertex(*v))
                .count()
        );
    }

    #[test]
    fn t_quad_box_delete_edge() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        qbox.delete_edge(
            qbox.halfedge_edge(
                qbox.find_halfedge(5.into(), 6.into())
                    .expect("Cannot find halfedge"),
            ),
            true,
            &mut cache.edges,
            &mut cache.vertices,
        )
        .expect("Cannot delete edge");
        assert_eq!(
            11,
            qbox.edges()
                .filter(|e| !qbox
                    .edge_status(*e)
                    .expect("Cannot read edge status")
                    .deleted())
                .count()
        );
        assert_eq!(
            4,
            qbox.faces()
                .filter(|f| !qbox
                    .face_status(*f)
                    .expect("Cannot read face status")
                    .deleted())
                .count()
        );
        qbox.garbage_collection(&mut cache)
            .expect("Failed to garbage collect");
        assert_eq!(qbox.num_faces(), 4);
        assert_eq!(qbox.num_edges(), 11);
        assert_eq!(qbox.num_vertices(), 8);
        assert_eq!(
            6,
            qbox.edges().filter(|e| qbox.is_boundary_edge(*e)).count()
        );
        assert_eq!(
            6,
            qbox.vertices()
                .filter(|v| qbox.is_boundary_vertex(*v))
                .count()
        );
    }

    #[test]
    fn t_quad_box_delete_vertex() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        qbox.delete_vertex(5.into(), true, &mut cache)
            .expect("Cannot delete face");
        for v in qbox.vertices() {
            assert_eq!(
                qbox.vertex_status(v)
                    .expect("Cannot read vertex status")
                    .deleted(),
                v.index() == 5
            );
        }
        qbox.garbage_collection(&mut cache)
            .expect("Garbage collection failed");
        assert_eq!(qbox.num_vertices(), 7);
        assert_eq!(qbox.num_edges(), 9);
        assert_eq!(qbox.num_halfedges(), 18);
        assert_eq!(qbox.num_faces(), 3);
        assert_eq!(
            6,
            qbox.edges().filter(|e| qbox.is_boundary_edge(*e)).count()
        );
        assert_eq!(
            6,
            qbox.vertices()
                .filter(|v| qbox.is_boundary_vertex(*v))
                .count()
        );
    }

    #[test]
    fn t_quad_box_delete_three_faces() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        for fi in 0..3u32 {
            qbox.delete_face(fi.into(), true, &mut cache.edges, &mut cache.vertices)
                .expect("Cannot delete faces");
        }
        qbox.garbage_collection(&mut cache)
            .expect("Failed to garbage collect");
        assert_eq!(3, qbox.num_faces());
        assert_eq!(9, qbox.num_edges());
        assert_eq!(18, qbox.num_halfedges());
        assert_eq!(7, qbox.num_vertices());
    }

    #[test]
    fn t_primitive_check_topology() {
        PolyMeshF32::tetrahedron(1.0)
            .expect("Cannot crate mesh")
            .check_topology()
            .expect("Topological check failed");
        PolyMeshF32::hexahedron(1.0)
            .expect("Cannot crate mesh")
            .check_topology()
            .expect("Topological check failed");
        PolyMeshF32::octahedron(1.0)
            .expect("Cannot crate mesh")
            .check_topology()
            .expect("Topological check failed");
        PolyMeshF32::icosahedron(1.0)
            .expect("Cannot crate mesh")
            .check_topology()
            .expect("Topological check failed");
        PolyMeshF32::dodecahedron(1.0)
            .expect("Cannot crate mesh")
            .check_topology()
            .expect("Topological check failed");
    }

    #[test]
    fn t_hexahedron_delete_faces_iter_mut() {
        let mut mesh = PolyMeshF32::hexahedron(1.0).expect("Cannot create hexahedron");
        let area = mesh.try_calc_area().expect("Cannot compute area");
        assert_eq!(8, mesh.num_vertices());
        assert_eq!(12, mesh.num_edges());
        assert_eq!(24, mesh.num_halfedges());
        assert_eq!(6, mesh.num_faces());
        for v in mesh.vertices() {
            // using a random condition to delete some faces, to make sure I can
            // modify the mesh. If you try changing any of the `m` inside this
            // loop to `mesh`, the borrow checker should complain and the code
            // should not compile.
            for (m, h) in mesh.voh_ccw_iter_mut(v) {
                if let Some(f) = m.halfedge_face(h) {
                    if (f.index() + h.index()) % 2 != 0 {
                        m.delete_face(f, true).expect("Cannot delete face");
                    }
                }
            }
        }
        mesh.garbage_collection()
            .expect("Garbage collection failed");
        mesh.check_topology().expect("Topological errors found");
        // It just so happens, we should be left with two adjacent faces meeting
        // at 90 degrees. This is because of the specific topology of the
        // hexahedron. I am just asserting the stuff below to ensure
        // determinimsm.
        assert_eq!(6, mesh.num_vertices());
        assert_eq!(7, mesh.num_edges());
        assert_eq!(14, mesh.num_halfedges());
        assert_eq!(2, mesh.num_faces());
        assert_f32_eq!(
            area / 3.0,
            mesh.try_calc_area().expect("Cannot compute area"),
            1e-6
        );
        assert_eq!(
            (6, 1),
            mesh.edges()
                .fold((0usize, 0usize), |(b, nb), e| if e.is_boundary(&mesh) {
                    (b + 1, nb)
                } else {
                    (b, nb + 1)
                })
        );
    }
}
