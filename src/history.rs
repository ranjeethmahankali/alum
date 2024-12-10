/*!
# Recording and Undoing Mesh Edits

Often algorithms require performing topological edits to a mesh. When something
goes wrong, it is quite useful to be able to go back to a previous state of the
mesh with a valid state. This module provides features that allow just that.

`TopolHistory` can be used to record topological edits being made to a mesh, and
that can be later used to revert said operations and restore the topology of the
mesh to an earlier checkpoint.
*/

use crate::{
    element::{Edge, Face, Halfedge, Vertex},
    topol::Topology,
    EPropBuf, EditableTopology, Error, FPropBuf, HPropBuf, Handle, PropBuf, Status, VPropBuf, EH,
    FH, HH, VH,
};

enum Element {
    Vertex(VH, Vertex, Status),
    Halfedge(HH, Halfedge, Status),
    EdgeStatus(EH, Status),
    Face(FH, Face, Status),
}

/// This can be used to remember and revert the topological edits being made to
/// a mesh.
///
/// ```rust
/// use alum::{use_glam::PolyMeshF32, TopolHistory, HasTopology, EditableTopology, HH};
///
/// let mut mesh = PolyMeshF32::unit_box().unwrap();
/// let mut history = TopolHistory::default();
/// let ckpt = history.check_point();
/// let h: HH = 0.into();
/// history.try_commit_edge_collapse(&mesh, h).unwrap();
/// let mut hcache = Vec::new();
/// mesh.try_collapse_edge(h, &mut hcache).unwrap();
/// {
///     let estatus = mesh.edge_status_prop();
///     let estatus = estatus.try_borrow().unwrap();
///     assert_eq!(11, mesh.edges().filter(|&e| !estatus[e].deleted()).count());
///     let vstatus = mesh.vertex_status_prop();
///     let vstatus = vstatus.try_borrow().unwrap();
///     assert_eq!(7, mesh.vertices().filter(|&v| !vstatus[v].deleted()).count());
/// }
/// assert!(history.try_restore(ckpt, &mut mesh).unwrap());
/// {
///     let estatus = mesh.edge_status_prop();
///     let estatus = estatus.try_borrow().unwrap();
///     assert_eq!(12, mesh.edges().filter(|&e| !estatus[e].deleted()).count());
///     let vstatus = mesh.vertex_status_prop();
///     let vstatus = vstatus.try_borrow().unwrap();
///     assert_eq!(8, mesh.vertices().filter(|&v| !vstatus[v].deleted()).count());
/// }
/// ```
#[derive(Default)]
pub struct TopolHistory {
    ops: Vec<Element>,
    vertices: Vec<VH>,
    halfedges: Vec<HH>,
    edges: Vec<EH>,
    faces: Vec<FH>,
}

impl TopolHistory {
    fn commit_vertices(&mut self, mesh: &Topology, vstatus: &VPropBuf<Status>) {
        self.vertices.sort();
        self.vertices.dedup();
        self.ops.extend(
            self.vertices
                .drain(..)
                .map(|v| Element::Vertex(v, *mesh.vertex(v), vstatus[v])),
        );
    }

    fn commit_edges(
        &mut self,
        mesh: &Topology,
        hstatus: &HPropBuf<Status>,
        estatus: &EPropBuf<Status>,
    ) {
        self.halfedges.sort();
        self.halfedges.dedup();
        self.edges.extend(self.halfedges.iter().map(|h| h.edge()));
        self.edges.sort();
        self.edges.dedup();
        self.ops.extend(
            self.halfedges
                .drain(..)
                .map(|h| Element::Halfedge(h, *mesh.halfedge(h), hstatus[h])),
        );
        self.ops.extend(
            self.edges
                .drain(..)
                .map(|e| Element::EdgeStatus(e, estatus[e])),
        );
    }

    fn commit_faces(&mut self, mesh: &Topology, fstatus: &FPropBuf<Status>) {
        self.faces.sort();
        self.faces.dedup();
        self.ops.extend(
            self.faces
                .drain(..)
                .map(|f| Element::Face(f, *mesh.face(f), fstatus[f])),
        );
    }

    fn push_degenerate_triangle_collapse(&mut self, mesh: &impl EditableTopology, h: HH) {
        let h1 = h.next(mesh);
        let o = h.opposite();
        let o1 = h1.opposite();
        let on = o.next(mesh);
        let op = o.prev(mesh);
        self.halfedges.extend_from_slice(&[h1, o, o1, on, op]);
        let v0 = h.head(mesh);
        let v1 = h1.head(mesh);
        self.vertices.extend_from_slice(&[v0, v1]);
        let fh = h.face(mesh);
        let fo = o.face(mesh);
        self.faces.extend([fh, fo].iter().filter_map(|f| *f));
    }

    fn commit(
        &mut self,
        mesh: &Topology,
        vstatus: &VPropBuf<Status>,
        hstatus: &HPropBuf<Status>,
        estatus: &EPropBuf<Status>,
        fstatus: &FPropBuf<Status>,
    ) {
        self.commit_vertices(mesh, vstatus);
        self.commit_edges(mesh, hstatus, estatus);
        self.commit_faces(mesh, fstatus);
    }

    /// This is similar to [`commit_edge_collapse`](Self::commit_edge_collapse),
    /// except it tries to borrow the required properties and returns an error
    /// if borrowing fails.
    pub fn try_commit_edge_collapse(
        &mut self,
        mesh: &impl EditableTopology,
        h: HH,
    ) -> Result<(), Error> {
        let (vstatus, hstatus, estatus, fstatus) = (
            mesh.vertex_status_prop(),
            mesh.halfedge_status_prop(),
            mesh.edge_status_prop(),
            mesh.face_status_prop(),
        );
        let (vstatus, hstatus, estatus, fstatus) = (
            vstatus.try_borrow()?,
            hstatus.try_borrow()?,
            estatus.try_borrow()?,
            fstatus.try_borrow()?,
        );
        self.commit_edge_collapse(mesh, h, &vstatus, &hstatus, &estatus, &fstatus);
        Ok(())
    }

    /// Commit an edge collapse to history.
    ///
    /// This encodes the necessary information into the history to later revert
    /// the edge collapse and restore the topology as it was before the edge
    /// collapse.
    pub fn commit_edge_collapse(
        &mut self,
        mesh: &impl EditableTopology,
        h: HH,
        vstatus: &VPropBuf<Status>,
        hstatus: &HPropBuf<Status>,
        estatus: &EPropBuf<Status>,
        fstatus: &FPropBuf<Status>,
    ) {
        // Cache neighborhood.
        let hn = h.next(mesh);
        let hp = h.prev(mesh);
        let o = h.opposite();
        let on = o.next(mesh);
        let op = o.prev(mesh);
        self.halfedges.extend_from_slice(&[h, hn, hp, o, on, op]);
        let fh = h.face(mesh);
        let fo = o.face(mesh);
        self.faces.extend([fh, fo].iter().filter_map(|f| *f));
        let vh = h.head(mesh);
        let vo = o.head(mesh);
        self.vertices.extend_from_slice(&[vh, vo]);
        self.halfedges.extend(mesh.vih_ccw_iter(vo));
        // Record degenerate triangle collapses.
        if hn.next(mesh) == hp {
            self.push_degenerate_triangle_collapse(mesh, hn);
        }
        if on.next(mesh) == op {
            self.push_degenerate_triangle_collapse(mesh, on);
        }
        // Commit to history.
        self.commit(mesh.topology(), vstatus, hstatus, estatus, fstatus);
    }

    /// This is similar to
    /// [`commit_edge_insertion`](Self::commit_edge_insertion), except it tries
    /// to borrow the necessary properties and returns an error if the borrowing
    /// fails.
    pub fn try_commit_edge_insertion(
        &mut self,
        mesh: &impl EditableTopology,
        prev: HH,
        next: HH,
    ) -> Result<(), Error> {
        let (vstatus, hstatus, estatus, fstatus) = (
            mesh.vertex_status_prop(),
            mesh.halfedge_status_prop(),
            mesh.edge_status_prop(),
            mesh.face_status_prop(),
        );
        let (vstatus, hstatus, estatus, fstatus) = (
            vstatus.try_borrow()?,
            hstatus.try_borrow()?,
            estatus.try_borrow()?,
            fstatus.try_borrow()?,
        );
        self.commit_edge_insertion(mesh, prev, next, (&vstatus, &hstatus, &estatus, &fstatus));
        Ok(())
    }

    /// Commit an edge insertion to history.
    ///
    /// This encodes the necessary information into the history to later revert
    /// the edge insertion and restore the previous topology of the mesh. Once
    /// reverted, the inserted edge and the corresponding halfedges will still
    /// be present in the mesh as a deleted elements.
    pub fn commit_edge_insertion(
        &mut self,
        mesh: &impl EditableTopology,
        prev: HH,
        next: HH,
        (vstatus, hstatus, estatus, fstatus): (
            &VPropBuf<Status>,
            &HPropBuf<Status>,
            &EPropBuf<Status>,
            &FPropBuf<Status>,
        ),
    ) {
        let (p1, n0) = (prev.next(mesh), next.prev(mesh));
        let f = prev.face(mesh);
        let v0 = prev.head(mesh);
        let v1 = next.tail(mesh);
        let ehnew: EH = (mesh.num_edges() as u32).into();
        let fhnew = f.map(|_| (mesh.num_faces() as u32).into());
        let new_edge = Edge {
            halfedges: [
                Halfedge {
                    face: f,
                    vertex: v1,
                    next,
                    prev,
                },
                Halfedge {
                    face: fhnew,
                    vertex: v0,
                    next: p1,
                    prev: n0,
                },
            ],
        };
        let (h0, h1) = ehnew.halfedges();
        let delstatus = {
            let mut s = Status::default();
            s.set_deleted(true);
            s
        };
        self.vertices.extend_from_slice(&[v0, v1]);
        self.halfedges.extend(mesh.loop_ccw_iter(prev));
        self.ops.push(Element::EdgeStatus(ehnew, delstatus));
        self.ops
            .push(Element::Halfedge(h0, new_edge.halfedges[0], delstatus));
        self.ops
            .push(Element::Halfedge(h1, new_edge.halfedges[1], delstatus));
        if let Some(fhnew) = fhnew {
            self.ops
                .push(Element::Face(fhnew, Face { halfedge: h1 }, delstatus));
        }
        if let Some(f) = f {
            self.faces.push(f);
        }
        self.commit_vertices(mesh.topology(), vstatus);
        self.commit_edges(mesh.topology(), hstatus, estatus);
        self.commit_faces(mesh.topology(), fstatus);
    }

    /// Get the checkpoint representing the current state of the history.
    pub fn check_point(&self) -> usize {
        self.ops.len()
    }

    /// This is similar to [`restore`](Self::restore) except it tries to borrow
    /// the necessary properties and returns an error if the borrowing fails.
    pub fn try_restore(
        &mut self,
        check_point: usize,
        mesh: &mut impl EditableTopology,
    ) -> Result<bool, Error> {
        let (mut vstatus, mut hstatus, mut estatus, mut fstatus) = (
            mesh.vertex_status_prop(),
            mesh.halfedge_status_prop(),
            mesh.edge_status_prop(),
            mesh.face_status_prop(),
        );
        let (mut vstatus, mut hstatus, mut estatus, mut fstatus) = (
            vstatus.try_borrow_mut().unwrap(),
            hstatus.try_borrow_mut().unwrap(),
            estatus.try_borrow_mut().unwrap(),
            fstatus.try_borrow_mut().unwrap(),
        );
        Ok(self.restore(
            check_point,
            mesh,
            &mut vstatus,
            &mut hstatus,
            &mut estatus,
            &mut fstatus,
        ))
    }

    /// Undo the commited operations on the given mesh and revert it's topology
    /// up to the provided checkpoint.
    ///
    /// If the checkpoint is valid, and the restoration occurs successfully
    /// `true` is returned. Otherwise `false` is returned.
    pub fn restore(
        &mut self,
        check_point: usize,
        mesh: &mut impl EditableTopology,
        vstatus: &mut VPropBuf<Status>,
        hstatus: &mut HPropBuf<Status>,
        estatus: &mut EPropBuf<Status>,
        fstatus: &mut FPropBuf<Status>,
    ) -> bool {
        if check_point >= self.ops.len() {
            return false;
        }
        let mesh = mesh.topology_mut();
        for op in self.ops[check_point..].iter().rev() {
            match op {
                Element::Vertex(vh, vertex, status) => {
                    *mesh.vertex_mut(*vh) = *vertex;
                    vstatus[*vh] = *status;
                }
                Element::Halfedge(hh, halfedge, status) => {
                    *mesh.halfedge_mut(*hh) = *halfedge;
                    hstatus[*hh] = *status;
                }
                Element::EdgeStatus(eh, status) => {
                    estatus[*eh] = *status;
                }
                Element::Face(fh, face, status) => {
                    *mesh.face_mut(*fh) = *face;
                    fstatus[*fh] = *status;
                }
            }
        }
        self.erase(check_point)
    }

    /// Erase the history up to the given checkpoint.
    ///
    /// This will remove all operations commited after the checkpoint from the
    /// history. This is meant to be used when a topological edit failed, and
    /// remembering that edit in the history is not necessary.
    pub fn erase(&mut self, checkpt: usize) -> bool {
        if checkpt < self.ops.len() {
            self.ops.truncate(checkpt);
            true
        } else {
            false
        }
    }
}

#[derive(Default)]
struct PropHistory<H, T>
where
    H: Handle,
    T: Copy + Clone + 'static,
{
    values: Vec<(H, T)>,
}

impl<H, T> PropHistory<H, T>
where
    H: Handle,
    T: Copy + Clone + 'static,
{
    pub fn commit(&mut self, handle: H, prop: &PropBuf<H, T>) {
        self.values.push((handle, prop[handle]));
    }

    pub fn check_point(&self) -> usize {
        self.values.len()
    }

    pub fn restore(&mut self, check_point: usize, prop: &mut PropBuf<H, T>) -> bool {
        if check_point >= self.values.len() {
            return false;
        }
        for (h, v) in self.values.drain(..).rev() {
            prop[h] = v;
        }
        true
    }
}

#[cfg(test)]
mod test {
    use super::TopolHistory;
    use crate::{use_glam::PolyMeshF32, EditableTopology, HasIterators, HasTopology, HH};
    use glam::vec3;

    fn compare_meshes(a: &PolyMeshF32, b: &PolyMeshF32) {
        // Compare vertices.
        let avs = a.vertex_status_prop();
        let avs = avs.try_borrow().unwrap();
        let bvs = b.vertex_status_prop();
        let bvs = bvs.try_borrow().unwrap();
        let apts = a.points();
        let apts = apts.try_borrow().unwrap();
        let bpts = b.points();
        let bpts = bpts.try_borrow().unwrap();
        let nv = usize::min(a.num_vertices(), b.num_vertices());
        assert!(a.vertices().skip(nv).all(|v| avs[v].deleted()));
        assert!(b.vertices().skip(nv).all(|v| bvs[v].deleted()));
        for (av, bv) in a.vertices().take(nv).zip(b.vertices().take(nv)) {
            print!("Checking vertex {}...", av);
            assert_eq!(a.topology().vertex(av), b.topology().vertex(bv));
            assert_eq!(avs[av], bvs[bv]);
            assert_eq!(apts[av], bpts[bv]);
            println!("✔ Passed");
        }
        // Compare edge status.
        let aes = a.edge_status_prop();
        let aes = aes.try_borrow().unwrap();
        let bes = b.edge_status_prop();
        let bes = bes.try_borrow().unwrap();
        let ne = usize::min(a.num_edges(), b.num_edges());
        assert!(a.edges().skip(ne).all(|e| aes[e].deleted()));
        assert!(b.edges().skip(ne).all(|e| bes[e].deleted()));
        for (ae, be) in a.edges().take(ne).zip(b.edges().take(ne)) {
            print!("Checking edge {}...", ae);
            assert_eq!(aes[ae], bes[be]);
            println!("✔ Passed");
        }
        // Compare halfedges.
        let ahs = a.halfedge_status_prop();
        let ahs = ahs.try_borrow().unwrap();
        let bhs = b.halfedge_status_prop();
        let bhs = bhs.try_borrow().unwrap();
        let nh = usize::min(a.num_halfedges(), b.num_halfedges());
        assert!(a.halfedges().skip(nh).all(|e| ahs[e].deleted()));
        assert!(b.halfedges().skip(nh).all(|e| bhs[e].deleted()));
        for (ah, bh) in a.halfedges().take(nh).zip(b.halfedges().take(nh)) {
            print!("Checking halfedge {}...", ah);
            assert_eq!(a.topology().halfedge(ah), b.topology().halfedge(bh));
            assert_eq!(ahs[ah], bhs[bh]);
            println!("✔ Passed");
        }
        // Compare faces.
        let afs = a.face_status_prop();
        let afs = afs.try_borrow().unwrap();
        let bfs = b.face_status_prop();
        let bfs = bfs.try_borrow().unwrap();
        let nf = usize::min(a.num_faces(), b.num_faces());
        assert!(a.faces().skip(nf).all(|f| afs[f].deleted()));
        assert!(b.faces().skip(nf).all(|f| bfs[f].deleted()));
        for (af, bf) in a.faces().take(nf).zip(b.faces().take(nf)) {
            print!("Checking face {}...", af);
            assert_eq!(a.topology().face(af), b.topology().face(bf));
            assert_eq!(afs[af], bfs[bf]);
            println!("✔ Passed");
        }
    }

    fn test_edge_collapse(mut mesh: PolyMeshF32, h: HH) {
        let mut history = TopolHistory::default();
        let area = mesh.try_calc_area().expect("Cannot compute area");
        let volume = mesh.try_calc_volume().expect("Cannot compute volume");
        let closed = mesh.is_closed();
        let copy = mesh.try_clone().expect("Cannot clone mesh");
        // Compute expected element counts after collapse.
        let is_triangle_1 = mesh.loop_ccw_iter(h).take(4).count() == 3;
        let is_triangle_2 = mesh.loop_ccw_iter(h.opposite()).take(4).count() == 3;
        let (nv, ne, nf) = (
            mesh.num_vertices() - 1,
            mesh.num_edges()
                - match (is_triangle_1, is_triangle_2) {
                    (true, true) => 3,
                    (true, false) | (false, true) => 2,
                    (false, false) => 1,
                },
            mesh.num_faces()
                - match (is_triangle_1, is_triangle_2) {
                    (true, true) => 2,
                    (true, false) | (false, true) => 1,
                    (false, false) => 0,
                },
        );
        {
            let checkpt = history.check_point();
            history
                .try_commit_edge_collapse(&mesh, h)
                .expect("Cannot commit edge collapse to history");
            let mut hcache = Vec::new();
            // Collapse and check.
            mesh.try_collapse_edge(h, &mut hcache)
                .expect("Cannot collapse edge");
            if closed {
                assert!(area > mesh.try_calc_area().expect("Cannot compute area"));
            }
            let newvol = mesh.try_calc_volume().expect("Cannot compute volume");
            if volume == 0.0 {
                assert_eq!(0.0, newvol);
            } else {
                assert!(volume > newvol);
            }
            {
                let (vstatus, estatus, fstatus) = (
                    mesh.vertex_status_prop(),
                    mesh.edge_status_prop(),
                    mesh.face_status_prop(),
                );
                let (vstatus, estatus, fstatus) = (
                    vstatus.try_borrow().expect("Cannot borrow property"),
                    estatus.try_borrow().expect("Cannot borrow property"),
                    fstatus.try_borrow().expect("Cannot borrow property"),
                );
                assert_eq!(
                    nv,
                    mesh.vertices().filter(|&v| !vstatus[v].deleted()).count()
                );
                assert_eq!(ne, mesh.edges().filter(|&e| !estatus[e].deleted()).count());
                assert_eq!(nf, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            }
            // Revert and check again.
            assert!(history
                .try_restore(checkpt, &mut mesh,)
                .expect("Cannot restore from history"));
        }
        // Nothing should be deleted.
        mesh.check_for_deleted()
            .expect("Unexpected deleted elements");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(area, mesh.try_calc_area().expect("Cannot compute area"));
        assert_eq!(
            volume,
            mesh.try_calc_volume().expect("Cannot compute volume")
        );
        compare_meshes(&mesh, &copy);
    }

    #[test]
    fn t_box_revert_edge_collapse() {
        let mesh = PolyMeshF32::unit_box().expect("Cannot make a box");
        let h = mesh
            .find_halfedge(0.into(), 1.into())
            .expect("Cannot find halfedge");
        test_edge_collapse(mesh, h);
    }

    #[test]
    fn t_octahedron_revert_edge_collapse() {
        let mesh = PolyMeshF32::octahedron(1.0).expect("Cannot make an icosahedron");
        let h = mesh
            .find_halfedge(0.into(), 1.into())
            .expect("Cannot find halfedge");
        test_edge_collapse(mesh, h);
    }

    #[test]
    fn t_icosahedron_revert_edge_collapse() {
        let mesh = PolyMeshF32::icosahedron(1.0).expect("Cannot make an icosahedron");
        let h = mesh
            .find_halfedge(0.into(), 1.into())
            .expect("Cannot find halfedge");
        test_edge_collapse(mesh, h);
    }

    #[test]
    fn t_triangles_revert_edge_collapse() {
        let mesh = {
            let mut mesh = PolyMeshF32::default();
            mesh.add_vertices(&[
                vec3(2.0, 2.0, 0.0),
                vec3(4.0, 2.0, 0.0),
                vec3(1.0, 0.0, 0.0),
                vec3(3.0, 0.0, 0.0),
                vec3(5.0, 0.0, 0.0),
                vec3(0.0, 2.0, 0.0),
                vec3(6.0, 2.0, 0.0),
                vec3(1.0, 4.0, 0.0),
                vec3(3.0, 4.0, 0.0),
                vec3(5.0, 4.0, 0.0),
            ])
            .expect("Cannot add vertices");
            for (a, b, c) in [
                (0u32, 5, 2),
                (0, 2, 3),
                (0, 3, 1),
                (1, 3, 4),
                (1, 4, 6),
                (0, 7, 5),
                (0, 8, 7),
                (0, 1, 8),
                (1, 9, 8),
                (1, 6, 9),
            ] {
                mesh.add_tri_face(a.into(), b.into(), c.into())
                    .expect("Cannot add face");
            }
            mesh
        };
        let h = mesh
            .find_halfedge(0.into(), 1.into())
            .expect("Cannot find halfedge");
        test_edge_collapse(mesh, h);
    }

    #[test]
    fn t_box_insert_edge_revert() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
        let copy = mesh.try_clone().expect("Cannot clone mesh");
        let area = mesh.try_calc_area().expect("Cannot compute area");
        let volume = mesh.try_calc_volume().expect("Cannot compute volume");
        let mut history = TopolHistory::default();
        let h = mesh
            .find_halfedge(4.into(), 5.into())
            .expect("Cannot find halfedge");
        {
            let ckpt = history.check_point();
            history
                .try_commit_edge_insertion(&mesh, h.next(&mesh), h)
                .expect("Cannot commit edge insertion to history");
            mesh.insert_edge(h.next(&mesh), h)
                .expect("Cannot insert edge");
            assert_eq!(8, mesh.num_vertices());
            assert_eq!(13, mesh.num_edges());
            assert_eq!(7, mesh.num_faces());
            assert_eq!(area, mesh.try_calc_area().expect("Cannot compute area"));
            assert_eq!(
                volume,
                mesh.try_calc_volume().expect("Cannot compute volume")
            );
            assert!(history
                .try_restore(ckpt, &mut mesh,)
                .expect("Cannot restore from history"));
        }
        mesh.check_topology().expect("Topological errors found");
        compare_meshes(&mesh, &copy);
        assert_eq!(area, mesh.try_calc_area().expect("Cannot compute area"));
        assert_eq!(
            volume,
            mesh.try_calc_volume().expect("Cannot compute volume")
        );
    }
}
