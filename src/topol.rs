use crate::{
    element::{Edge, Face, Halfedge, Handle, Vertex, EH, FH, HH, VH},
    error::Error,
    iterator,
    property::{EProperty, FProperty, HProperty, PropertyContainer, TPropData, VProperty},
};

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
pub(crate) struct TopolCache {
    loop_halfedges: Vec<Option<HH>>,
    needs_adjust: Vec<bool>,
    next_cache: Vec<(HH, HH)>,
    tentative: Vec<TentativeEdge>,
    halfedges: Vec<HH>,
}

impl TopolCache {
    fn clear(&mut self) {
        self.loop_halfedges.clear();
        self.needs_adjust.clear();
        self.next_cache.clear();
        self.tentative.clear();
        self.halfedges.clear();
    }
}

pub(crate) struct Topology {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    vprops: PropertyContainer,
    hprops: PropertyContainer,
    eprops: PropertyContainer,
    fprops: PropertyContainer,
}

impl Topology {
    pub fn new() -> Self {
        Topology {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            vprops: PropertyContainer::new(),
            hprops: PropertyContainer::new(),
            eprops: PropertyContainer::new(),
            fprops: PropertyContainer::new(),
        }
    }

    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        Topology {
            vertices: Vec::with_capacity(nverts),
            edges: Vec::with_capacity(nedges),
            faces: Vec::with_capacity(nfaces),
            vprops: PropertyContainer::new(),
            hprops: PropertyContainer::new(),
            eprops: PropertyContainer::new(),
            fprops: PropertyContainer::new(),
        }
    }

    pub fn create_vertex_prop<T: TPropData>(&mut self) -> VProperty<T> {
        VProperty::<T>::new(&mut self.vprops)
    }

    pub fn create_halfedge_prop<T: TPropData>(&mut self) -> HProperty<T> {
        HProperty::<T>::new(&mut self.hprops)
    }

    pub fn create_edge_prop<T: TPropData>(&mut self) -> EProperty<T> {
        EProperty::<T>::new(&mut self.eprops)
    }

    pub fn create_face_prop<T: TPropData>(&mut self) -> FProperty<T> {
        FProperty::<T>::new(&mut self.fprops)
    }

    fn vertex(&self, v: VH) -> &Vertex {
        &self.vertices[v.index() as usize]
    }

    fn halfedge(&self, h: HH) -> &Halfedge {
        &self.edges[(h.index() >> 1) as usize].halfedges[(h.index() & 1) as usize]
    }

    fn halfedge_mut(&mut self, h: HH) -> &mut Halfedge {
        &mut self.edges[(h.index() >> 1) as usize].halfedges[(h.index() & 1) as usize]
    }

    pub fn vertex_halfedge(&self, v: VH) -> Option<HH> {
        self.vertex(v).halfedge
    }

    pub fn to_vertex(&self, h: HH) -> VH {
        self.halfedge(h).vertex
    }

    pub fn from_vertex(&self, h: HH) -> VH {
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

    pub fn edge_halfedge(&self, e: EH, flag: bool) -> HH {
        ((e.index() << 1) & if flag { 1 } else { 0 }).into()
    }

    pub fn face_halfedge(&self, f: FH) -> HH {
        self.faces[f.index() as usize].halfedge
    }

    pub fn is_boundary_halfedge(&self, h: HH) -> bool {
        self.halfedge(h).face.is_none()
    }

    pub fn is_boundary_edge(&self, e: EH) -> bool {
        let h = (e.index() << 1).into();
        self.is_boundary_halfedge(h) || self.is_boundary_halfedge(self.opposite_halfedge(h))
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
        iterator::voh_ccw_iter(self, from).find(|h| self.to_vertex(*h) == to)
    }

    pub fn is_manifold_vertex(&self, v: VH) -> bool {
        /* If just the first outgoing halfedge is on the boundary, it just means
         * the vertex is on the boundary. If the first outgoing halfedge is not
         * on the boundary, it implies the vertex is in the interior. In both
         * cases the vertex is manifold. If any outgoing halfedge apart from the
         * first is on the boundary, it implies there are more than one gaps
         * when circulating around the vertex, making it non-manifold. For this
         * reason, we skip the first halfedge and check the rest.
         */
        iterator::voh_ccw_iter(self, v)
            .skip(1)
            .all(|h| !self.is_boundary_halfedge(h))
    }

    fn adjust_outgoing_halfedge(&mut self, v: VH) {
        let h = iterator::voh_ccw_iter(self, v).find(|h| self.is_boundary_halfedge(*h));
        if let Some(h) = h {
            self.set_vertex_halfedge(v, h)
        }
    }

    pub fn add_vertex(&mut self) -> Result<VH, Error> {
        let vi = self.vertices.len() as u32;
        self.vprops.push_value()?;
        self.vertices.push(Vertex { halfedge: None });
        Ok(vi.into())
    }

    fn new_edge(
        &mut self,
        from: VH,
        to: VH,
        prev: HH,
        next: HH,
        opp_prev: HH,
        opp_next: HH,
    ) -> Result<u32, Error> {
        let ei = self.edges.len() as u32;
        self.eprops.push_value()?;
        for _ in 0..2 {
            self.hprops.push_value()?;
        }
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
        Ok(ei)
    }

    fn new_face(&mut self, halfedge: HH) -> Result<FH, Error> {
        let fi = self.faces.len() as u32;
        self.fprops.push_value()?;
        self.faces.push(Face { halfedge });
        Ok(fi.into())
    }

    pub fn set_vertex_halfedge(&mut self, v: VH, h: HH) {
        self.vertices[v.index() as usize].halfedge = Some(h);
    }

    pub fn set_next_halfedge(&mut self, hprev: HH, hnext: HH) {
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
                    self.set_vertex_halfedge(v, outernext);
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
                    self.set_vertex_halfedge(v, boundnext);
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
                        self.set_vertex_halfedge(v, outernext);
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
                    assert_eq!(*index >> 1, ei, "Failed to create an edge loop");
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
            self.set_next_halfedge(prev, next);
        }
        // Adjust vertices' halfedge handles.
        for i in 0..verts.len() {
            if cache.needs_adjust[i] {
                self.adjust_outgoing_halfedge(verts[i]);
            }
        }
        Ok(fnew)
    }

    pub fn vertex_valence(&self, v: VH) -> usize {
        iterator::voh_ccw_iter(self, v).count()
    }

    pub fn face_valence(&self, f: FH) -> usize {
        iterator::fv_ccw_iter(self, f).count()
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod test {
    use arrayvec::ArrayVec;

    use crate::{
        iterator,
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
    fn quad_box() -> Topology {
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

    fn loop_mesh() -> Topology {
        /*

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
            assert_eq!(qbox.vertex_valence(v), 3);
        }
    }

    #[test]
    fn t_box_face_valence() {
        let qbox = quad_box();
        for f in qbox.faces() {
            assert_eq!(qbox.face_valence(f), 4);
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
            assert_eq!(mesh.from_vertex(h), 5.into());
            assert_eq!(mesh.to_vertex(h), 6.into());
            let h1 = mesh.next_halfedge(h);
            assert_eq!(
                h1,
                mesh.find_halfedge(6.into(), v1)
                    .expect("Cannot find halfedge")
            );
            let h = mesh
                .find_halfedge(6.into(), 10.into())
                .expect("Cannot find halfedge");
            assert_eq!(mesh.from_vertex(h), 6.into());
            assert_eq!(mesh.to_vertex(h), 10.into());
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
}
