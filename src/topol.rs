use crate::iterator::HasIterators;
use crate::property::EProperty;
use crate::{
    element::{
        ERange, Edge, FRange, Face, HRange, Halfedge, Handle, VRange, Vertex, EH, FH, HH, VH,
    },
    error::Error,
    property::{FProperty, HProperty, PropertyContainer, VProperty},
    status::Status,
};
use std::cell::RefMut;

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

/// This trait defines functions to access the topology of a mesh.
///
/// It also has other related functions to create properties, access the status
/// of an element etc.
pub trait HasTopology: Sized {
    fn topology(&self) -> &Topology;

    fn topology_mut(&mut self) -> &mut Topology;

    /// Number of vertices.
    fn num_vertices(&self) -> usize {
        self.topology().vertices.len()
    }

    /// Number of halfedges.
    fn num_halfedges(&self) -> usize {
        self.num_edges() * 2
    }

    /// Number of edges.
    fn num_edges(&self) -> usize {
        self.topology().edges.len()
    }

    /// Number of faces.
    fn num_faces(&self) -> usize {
        self.topology().faces.len()
    }

    /// Create a new vertex property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_vertex_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_vertices(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    fn create_vertex_prop<T>(&mut self, default: T) -> VProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        VProperty::<T>::new(&mut self.topology_mut().vprops, default)
    }

    /// Create a new halfedge property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_halfedge_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_halfedges(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    fn create_halfedge_prop<T>(&mut self, default: T) -> HProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        HProperty::<T>::new(&mut self.topology_mut().hprops, default)
    }

    /// Create a new edge property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_edge_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_edges(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    fn create_edge_prop<T>(&mut self, default: T) -> EProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        EProperty::<T>::new(&mut self.topology_mut().eprops, default)
    }

    /// Create a new face property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_face_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_faces(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    fn create_face_prop<T>(&mut self, default: T) -> FProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        FProperty::<T>::new(&mut self.topology_mut().fprops, default)
    }

    /// Get the number of properties defined on vertices.
    ///
    /// This may include empty slots corresponding to properties that have been
    /// dropped. Performing a garbage collection can clean up all empty slots
    /// and provide an accurate count of the properties.
    fn num_vertex_props(&self) -> usize {
        self.topology().vprops.num_properties()
    }

    /// Get the number of properties defined on halfedges.
    ///
    /// This may include empty slots corresponding to properties that have been
    /// dropped. Performing a garbage collection can clean up all empty slots
    /// and provide an accurate count of the properties.
    fn num_halfedge_props(&self) -> usize {
        self.topology().hprops.num_properties()
    }

    /// Get the number of properties defined on edges.
    ///
    /// This may include empty slots corresponding to properties that have been
    /// dropped. Performing a garbage collection can clean up all empty slots
    /// and provide an accurate count of the properties.
    fn num_edge_props(&self) -> usize {
        self.topology().eprops.num_properties()
    }

    /// Get the number of properties defined on faces.
    ///
    /// This may include empty slots corresponding to properties that have been
    /// dropped. Performing a garbage collection can clean up all empty slots
    /// and provide an accurate count of the properties.
    fn num_face_props(&self) -> usize {
        self.topology().fprops.num_properties()
    }

    /// Reserve memory for the given number of elements.
    ///
    /// The memory is also reserved for all properties.
    fn reserve(&mut self, nverts: usize, nedges: usize, nfaces: usize) -> Result<(), Error> {
        // Elements.
        self.topology_mut().vertices.reserve(nverts);
        self.topology_mut().edges.reserve(nedges);
        self.topology_mut().faces.reserve(nfaces);
        // Properties.
        self.topology_mut().vprops.reserve(nverts)?;
        self.topology_mut().hprops.reserve(nedges * 2)?;
        self.topology_mut().eprops.reserve(nedges)?;
        self.topology_mut().fprops.reserve(nfaces)?;
        Ok(())
    }

    /// Delete all elements and their properties.
    fn clear(&mut self) -> Result<(), Error> {
        // Elements.
        self.topology_mut().vertices.clear();
        self.topology_mut().edges.clear();
        self.topology_mut().faces.clear();
        // Properties.
        self.topology_mut().vprops.clear()?;
        self.topology_mut().hprops.clear()?;
        self.topology_mut().eprops.clear()?;
        self.topology_mut().fprops.clear()?;
        Ok(())
    }

    /// Iterator over the vertices of the mesh.
    fn vertices(&self) -> VRange {
        (0..(self.num_vertices() as u32)).into()
    }

    /// Iterator over the halfedges of the mesh.
    fn halfedges(&self) -> HRange {
        (0..(self.num_halfedges() as u32)).into()
    }

    /// Iterator over the edges of the mesh.
    fn edges(&self) -> ERange {
        (0..(self.num_edges() as u32)).into()
    }

    /// Iterator over the faces of the mesh.
    fn faces(&self) -> FRange {
        (0..(self.num_faces() as u32)).into()
    }

    /// Get the vertex status property.
    fn vertex_status_prop(&self) -> VProperty<Status> {
        self.topology().vstatus.clone()
    }

    /// Get the halfedge status property.
    fn halfedge_status_prop(&self) -> HProperty<Status> {
        self.topology().hstatus.clone()
    }

    /// Get the edge status property.
    fn edge_status_prop(&self) -> EProperty<Status> {
        self.topology().estatus.clone()
    }

    /// Get the face status property.
    fn face_status_prop(&self) -> FProperty<Status> {
        self.topology().fstatus.clone()
    }

    /// The status of a vertex.
    fn vertex_status(&self, v: VH) -> Result<Status, Error> {
        self.topology().vstatus.get_cloned(v)
    }

    /// The status of a halfedge.
    fn halfedge_status(&self, h: HH) -> Result<Status, Error> {
        self.topology().hstatus.get_cloned(h)
    }

    /// The status of an edge.
    fn edge_status(&self, e: EH) -> Result<Status, Error> {
        self.topology().estatus.get_cloned(e)
    }

    /// The status of a face.
    fn face_status(&self, f: FH) -> Result<Status, Error> {
        self.topology().fstatus.get_cloned(f)
    }

    /// The status of a vertex as mutable.
    fn vertex_status_mut(&mut self, v: VH) -> Result<RefMut<'_, Status>, Error> {
        self.topology_mut().vstatus.get_mut(v)
    }

    /// The status of a halfedge as mutable.
    fn halfedge_status_mut(&mut self, h: HH) -> Result<RefMut<'_, Status>, Error> {
        self.topology_mut().hstatus.get_mut(h)
    }

    /// The status of an edge as mutable.
    fn edge_status_mut(&mut self, e: EH) -> Result<RefMut<'_, Status>, Error> {
        self.topology_mut().estatus.get_mut(e)
    }

    /// The status of a face as mutable.
    fn face_status_mut(&mut self, f: FH) -> Result<RefMut<'_, Status>, Error> {
        self.topology_mut().fstatus.get_mut(f)
    }

    /// Check the topology of the mesh.
    ///
    /// This function will return an error if any errors are found in the topolgy.
    fn check_topology(&self) -> Result<(), Error> {
        self.topology().check()
    }

    /// Delete the isolated vertices of this mesh.
    ///
    /// A vertex is isolated if it has no incident edges. If successful, the
    /// number of vertices deleted is returned. This function may fail when
    /// trying to borrow properties.
    fn delete_isolated_vertices(&mut self) -> Result<usize, Error> {
        let mut vstatus = self.vertex_status_prop();
        let mut vstatus = vstatus.try_borrow_mut()?;
        Ok(self
            .vertices()
            .fold(0usize, |count, v| match v.halfedge(self) {
                Some(_) => count,
                None => {
                    vstatus[v].set_deleted(true);
                    count + 1
                }
            }))
    }

    /// Check for deleted elements, and return an error if any are found.
    ///
    /// This is useful for making sure there are no deleted elements. If any are
    /// found, calling garbage collection will clean them up.
    fn check_for_deleted(&self) -> Result<(), Error> {
        // Ensure no deleted topology.
        let fstatus = self.face_status_prop();
        let fstatus = fstatus.try_borrow()?;
        if let Some(f) = self.faces().find(|f| fstatus[*f].deleted()) {
            return Err(Error::DeletedFace(f));
        }
        let estatus = self.edge_status_prop();
        let estatus = estatus.try_borrow()?;
        if let Some(e) = self.edges().find(|e| estatus[*e].deleted()) {
            return Err(Error::DeletedEdge(e));
        }
        let vstatus = self.vertex_status_prop();
        let vstatus = vstatus.try_borrow()?;
        if let Some(v) = self.vertices().find(|v| vstatus[*v].deleted()) {
            return Err(Error::DeletedVertex(v));
        }
        Ok(())
    }

    /// Check if the mesh is closed, i.e. has no boundary edges.
    ///
    /// A mesh is closed edges have exactly two incident faces.
    fn is_closed(&self) -> bool {
        self.halfedges().all(|h| h.is_boundary(self))
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

    pub(crate) fn vertex(&self, v: VH) -> &Vertex {
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

    pub(crate) fn face(&self, f: FH) -> &Face {
        &self.faces[f.index() as usize]
    }

    pub(crate) fn face_mut(&mut self, f: FH) -> &mut Face {
        &mut self.faces[f.index() as usize]
    }

    pub fn add_vertex(&mut self) -> Result<VH, Error> {
        let vi = self.vertices.len() as u32;
        self.vprops.push_value()?;
        self.vertices.push(Vertex { halfedge: None });
        Ok(vi.into())
    }

    pub fn add_vertices(&mut self, n: usize) -> Result<VRange, Error> {
        let nverts = self.vertices.len() as u32;
        self.vprops.push_values(n)?;
        self.vertices
            .resize(self.vertices.len() + n, Vertex { halfedge: None });
        Ok((nverts..(nverts + n as u32)).into())
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
            if !verts[i].is_boundary(self) {
                // Ensure vertex is manifold.
                return Err(Error::ComplexVertex(verts[i]));
            }
            // Ensure edge is manifold.
            let h = self.find_halfedge(verts[i], verts[(i + 1) % verts.len()]);
            match h {
                Some(h) if !h.is_boundary(self) => return Err(Error::ComplexHalfedge(h)),
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
                (Some(prev), Some(next)) if prev.next(self) != next => Some((prev, next)),
                _ => None,
            }
        }) {
            // Relink the patch.
            let boundprev = {
                let mut out = next.opposite();
                loop {
                    out = out.next(self).opposite();
                    if out.is_boundary(self) {
                        break;
                    }
                }
                out
            };
            let boundnext = boundprev.next(self);
            // Ok ?
            if boundprev == prev {
                return Err(Error::PatchRelinkingFailed);
            }
            debug_assert!(boundprev.is_boundary(self) && boundnext.is_boundary(self));
            // other halfedges.
            let pstart = prev.next(self);
            let pend = next.prev(self);
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
                    cache.needs_adjust[j] = v.halfedge(self) == Some(*innernext);
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
                    let innerprev: HH = innerprev.into();
                    let outernext = innerprev.opposite();
                    let boundprev = innernext.prev(self);
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
                    let innernext: HH = innernext.into();
                    let outerprev = innernext.opposite();
                    let boundnext = innerprev.next(self);
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
                    let innerprev: HH = innerprev.into();
                    let innernext: HH = innernext.into();
                    let outernext = innerprev.opposite();
                    let outerprev = innernext.opposite();
                    if let Some(boundnext) = v.halfedge(self) {
                        let boundprev = boundnext.prev(self);
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
        cache.faces.extend(self.vf_ccw_iter(v));
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
        let h0 = e.halfedge(false);
        let h1 = e.halfedge(true);
        let f0 = h0.face(self);
        let f1 = h1.face(self);
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
                hstatus[h0].set_deleted(true);
                hstatus[h1].set_deleted(true);
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
        for (mesh, h) in self.fh_ccw_iter_mut(f) {
            mesh.halfedge_mut(h).face = None; // Disconnect from face.
            if h.opposite().is_boundary(mesh) {
                ecache.push(h.edge());
            }
            vcache.push(h.head(mesh));
        }
        // Delete collected topology.
        for e in ecache.drain(..) {
            let h0 = e.halfedge(false);
            let v0 = h0.head(self);
            let next0 = h0.next(self);
            let prev0 = h0.prev(self);
            let h1 = e.halfedge(true);
            let v1 = h1.head(self);
            let next1 = h1.next(self);
            let prev1 = h1.prev(self);
            // Adjust halfedge links and mark edge and halfedges deleted.
            self.link_halfedges(prev0, next1);
            self.link_halfedges(prev1, next0);
            self.edge_status_mut(e)?.set_deleted(true);
            {
                let mut hstatus = self.hstatus.try_borrow_mut()?;
                hstatus[h0].set_deleted(true);
                hstatus[h1].set_deleted(true);
            }
            // Update vertices.
            if v0.halfedge(self) == Some(h1) {
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
            if v1.halfedge(self) == Some(h0) {
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
                    let status: &[Status] = &status;
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
                    let status: &[Status] = &status;
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
                    let status: &[Status] = &status;
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
            self.link_halfedges(h, hmap[h.next(self).index() as usize]);
            if let Some(f) = h.face(self) {
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

impl Clone for Topology {
    /// Clone the topology. This is not a full clone.
    ///
    /// This clones the topology, the elements and the built-in properties such
    /// as statuses for vertices, halfedges, edges, and faces. This doesn't
    /// clone any other properties because the topology doesn't fully own the
    /// other properties. It is upto the owners of the property (i.e. users of
    /// the API) to clone the other properties. Furthermore, the built-in
    /// properties are only cloned if they can be successfully borrowed. If the
    /// caller already borrowed the builtin properties, they may not be copied.
    fn clone(&self) -> Self {
        let (mut vprops, mut hprops, mut eprops, mut fprops) = (
            PropertyContainer::new_with_size(self.num_vertices()),
            PropertyContainer::new_with_size(self.num_halfedges()),
            PropertyContainer::new_with_size(self.num_edges()),
            PropertyContainer::new_with_size(self.num_faces()),
        );
        let (mut vstatus, mut hstatus, mut estatus, mut fstatus) = (
            VProperty::new(&mut vprops, Default::default()),
            HProperty::new(&mut hprops, Default::default()),
            EProperty::new(&mut eprops, Default::default()),
            FProperty::new(&mut fprops, Default::default()),
        );
        match (self.vstatus.try_borrow(), vstatus.try_borrow_mut()) {
            (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
            _ => {} //
        }
        match (self.hstatus.try_borrow(), hstatus.try_borrow_mut()) {
            (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
            _ => {} //
        }
        match (self.estatus.try_borrow(), estatus.try_borrow_mut()) {
            (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
            _ => {} //
        }
        match (self.fstatus.try_borrow(), fstatus.try_borrow_mut()) {
            (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
            _ => {} //
        }
        Self {
            vertices: self.vertices.clone(),
            edges: self.edges.clone(),
            faces: self.faces.clone(),
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
}

impl HasTopology for Topology {
    fn topology(&self) -> &Topology {
        self
    }

    fn topology_mut(&mut self) -> &mut Topology {
        self
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
        iterator::HasIterators,
        macros::assert_f32_eq,
        topol::{Handle, HasTopology, VH},
        use_glam::PolyMeshF32,
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
            let h = v
                .halfedge(&topol)
                .expect("Vertex must have an incident halfedge");
            assert!(h.is_boundary(&topol));
            let oh = h.opposite();
            assert!(!oh.is_boundary(&topol));
            assert_eq!(
                oh.face(&topol)
                    .expect("Halfedge must have an incident face"),
                face
            );
        }
        assert_eq!(
            topol.halfedges().filter(|h| h.is_boundary(&topol)).count(),
            3
        );
        assert_eq!(
            topol.halfedges().filter(|h| !h.is_boundary(&topol)).count(),
            3
        );
        for (i, j) in (0u32..3).map(|i| (i, (i + 1) % 3)) {
            let h = topol.find_halfedge(i.into(), j.into()).unwrap();
            assert!(!h.is_boundary(&topol));
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
        assert_eq!(topol.edges().filter(|e| e.is_boundary(&topol)).count(), 4);
        assert_eq!(topol.edges().filter(|e| !e.is_boundary(&topol)).count(), 1);
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
            let h = v
                .halfedge(&topol)
                .expect("Vertex must have an incident halfedge");
            assert!(h.is_boundary(&topol));
            let oh = h.opposite();
            assert!(!oh.is_boundary(&topol));
            assert_eq!(
                oh.face(&topol)
                    .expect("Halfedge must have an incident face"),
                face
            );
        }
        assert_eq!(
            topol.halfedges().filter(|h| h.is_boundary(&topol)).count(),
            4
        );
        assert_eq!(
            topol.halfedges().filter(|h| !h.is_boundary(&topol)).count(),
            4
        );
    }

    #[test]
    fn t_box_manifold() {
        let qbox = quad_box();
        assert!(
            qbox.halfedges().all(|h| !h.is_boundary(&qbox)),
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
        assert_eq!(mesh.edges().filter(|e| e.is_boundary(&mesh)).count(), 16);
        // Add a floating triangle at v6.
        let v0 = mesh.add_vertex().expect("Unable to add vertex");
        let v1 = mesh.add_vertex().expect("Unable to add vertex");
        let mut cache = TopolCache::default();
        let f0 = mesh
            .add_face(&[6.into(), v0, v1], &mut cache)
            .expect("Unable to add a face");
        assert_eq!(f0.index(), 8);
        assert_eq!(
            mesh.vf_ccw_iter(6.into())
                .map(|i| i.index())
                .collect::<Vec<_>>(),
            [8, 1, 2, 4]
        );
        assert_eq!(
            mesh.vv_ccw_iter(6.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            [10, v0.index(), v1.index(), 5, 2, 7]
        );
        assert_eq!(mesh.ve_ccw_iter(6.into()).count(), 6);
        assert_eq!(mesh.vertices().filter(|v| v.is_manifold(&mesh)).count(), 17);
        let v: VH = 6.into();
        assert!(!v.is_manifold(&mesh));
        // Check halfedge connectivity at v6.
        {
            let h = mesh
                .find_halfedge(5.into(), 6.into())
                .expect("Cannot find halfedge");
            assert_eq!(h.tail(&mesh), 5.into());
            assert_eq!(h.head(&mesh), 6.into());
            let h1 = h.next(&mesh);
            assert_eq!(
                h1,
                mesh.find_halfedge(6.into(), v1)
                    .expect("Cannot find halfedge")
            );
            let h = mesh
                .find_halfedge(6.into(), 10.into())
                .expect("Cannot find halfedge");
            assert_eq!(h.tail(&mesh), 6.into());
            assert_eq!(h.head(&mesh), 10.into());
            let h = h.prev(&mesh);
            assert_eq!(
                h,
                mesh.find_halfedge(v0, 6.into())
                    .expect("Cannot find halfedge")
            );
        }
        assert_eq!(mesh.edges().filter(|e| e.is_boundary(&mesh)).count(), 19);
        // Add another face.
        let f1 = mesh
            .add_face(&[5.into(), 6.into(), v1], &mut cache)
            .expect("Unable to add face");
        assert_eq!(f1.index(), 9);
        assert_eq!(mesh.num_edges(), 28);
        assert_eq!(mesh.edges().filter(|e| e.is_boundary(&mesh)).count(), 18);
        assert_eq!(
            mesh.vf_ccw_iter(6.into())
                .map(|f| f.index())
                .collect::<Vec<_>>(),
            [8, 9, 1, 2, 4]
        );
        assert_eq!(
            mesh.vf_ccw_iter(5.into())
                .map(|f| f.index())
                .collect::<Vec<_>>(),
            [3, 0, 1, 9]
        );
        assert_eq!(
            mesh.vv_ccw_iter(6.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            [10, v0.index(), v1.index(), 5, 2, 7]
        );
        assert_eq!(
            mesh.vv_ccw_iter(5.into())
                .map(|v| v.index())
                .collect::<Vec<_>>(),
            [v1.index(), 9, 4, 1, 6]
        );
        let (v0, v1): (VH, VH) = (6.into(), 5.into());
        assert!(v0.is_manifold(&mesh));
        assert!(v1.is_manifold(&mesh));
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
            assert_eq!(v.is_boundary(&qbox), v.index() > 3);
        }
        assert_eq!(4, qbox.edges().filter(|e| e.is_boundary(&qbox)).count());
        assert_eq!(4, qbox.vertices().filter(|v| v.is_boundary(&qbox)).count());
    }

    #[test]
    fn t_quad_box_delete_edge() {
        let mut qbox = quad_box();
        let mut cache = TopolCache::default();
        qbox.delete_edge(
            qbox.find_halfedge(5.into(), 6.into())
                .expect("Cannot find halfedge")
                .edge(),
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
        assert_eq!(6, qbox.edges().filter(|e| e.is_boundary(&qbox)).count());
        assert_eq!(6, qbox.vertices().filter(|v| v.is_boundary(&qbox)).count());
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
        assert_eq!(6, qbox.edges().filter(|e| e.is_boundary(&qbox)).count());
        assert_eq!(6, qbox.vertices().filter(|v| v.is_boundary(&qbox)).count());
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
                if let Some(f) = h.face(m) {
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

    #[test]
    fn t_quad_box_clone() {
        let mut mesh = quad_box();
        let myprop = mesh.create_halfedge_prop(0u8); // Create custom halfedge property.
        {
            let myprop = myprop.try_borrow().expect("Cannot borrow property");
            assert_eq!(myprop.len(), mesh.num_halfedges());
        }
        // Tag the odd numbered elements.
        for v in mesh.vertices().filter(|v| v.index() % 2 != 0) {
            mesh.vertex_status_mut(v)
                .expect("Cannot access vertex status")
                .set_tagged(true);
        }
        for h in mesh.halfedges().filter(|h| h.index() % 2 != 0) {
            mesh.halfedge_status_mut(h)
                .expect("Cannot access halfedge status")
                .set_tagged(true);
        }
        for e in mesh.edges().filter(|e| e.index() % 2 != 0) {
            mesh.edge_status_mut(e)
                .expect("Cannot access edge status")
                .set_tagged(true);
        }
        for f in mesh.faces().filter(|f| f.index() % 2 != 0) {
            mesh.face_status_mut(f)
                .expect("Cannot access face status")
                .set_tagged(true);
        }
        let copy = mesh.clone();
        for v in copy.vertices() {
            assert_eq!(
                v.index() % 2 != 0,
                copy.vertex_status(v)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        for h in copy.halfedges() {
            assert_eq!(
                h.index() % 2 != 0,
                copy.halfedge_status(h)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        for e in copy.edges() {
            assert_eq!(
                e.index() % 2 != 0,
                copy.edge_status(e)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        for f in copy.faces() {
            assert_eq!(
                f.index() % 2 != 0,
                copy.face_status(f)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        // Caller owned properties are not cloned. There is exactly one caller owned property.
        assert_eq!(1 + copy.num_halfedge_props(), mesh.num_halfedge_props());
    }
}
