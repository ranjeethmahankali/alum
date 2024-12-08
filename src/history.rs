use crate::{
    element::{Face, Halfedge, Vertex},
    topol::Topology,
    EPropBuf, EditableTopology, FPropBuf, HPropBuf, Status, VPropBuf, EH, FH, HH, VH,
};

enum Element {
    Vertex(VH, Vertex, Status),
    Halfedge(HH, Halfedge, Status),
    EdgeStatus(EH, Status),
    Face(FH, Face, Status),
}

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
        self.edges.clear();
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
        self.halfedges.extend_from_slice(&[h1, o, o1]);
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
        self.halfedges.extend_from_slice(&[hn, hp, o, on, op]);
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

    pub fn check_point(&self) -> usize {
        self.ops.len()
    }

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
        self.ops.truncate(check_point);
        true
    }
}

#[cfg(test)]
mod test {
    use super::TopolHistory;
    use crate::{use_glam::PolyMeshF32, EditableTopology, HasIterators, HasTopology};

    #[test]
    fn t_box_revert_edge_collapse() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot make a box");
        let h = mesh
            .find_halfedge(0.into(), 1.into())
            .expect("Cannot find halfedge");
        let mut history = TopolHistory::default();
        let area = mesh.try_calc_area().expect("Cannot compute area");
        let volume = mesh.try_calc_volume().expect("Cannot compute volume");
        {
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
            let checkpt = history.check_point();
            history.commit_edge_collapse(&mesh, h, &vstatus, &hstatus, &estatus, &fstatus);
            let mut hcache = Vec::new();
            // Collapse and check.
            mesh.collapse_edge(
                h,
                &mut vstatus,
                &mut hstatus,
                &mut estatus,
                &mut fstatus,
                &mut hcache,
            );
            assert!(area > mesh.try_calc_area().expect("Cannot compute area"));
            assert!(volume > mesh.try_calc_volume().expect("Cannot compute volume"));
            assert_eq!(
                7,
                mesh.vertices().filter(|&v| !vstatus[v].deleted()).count()
            );
            assert_eq!(11, mesh.edges().filter(|&e| !estatus[e].deleted()).count());
            assert_eq!(6, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            // Revert and check again.
            assert!(history.restore(
                checkpt,
                &mut mesh,
                &mut vstatus,
                &mut hstatus,
                &mut estatus,
                &mut fstatus
            ));
        }
        // Nothing should be deleted.
        mesh.check_for_deleted()
            .expect("Unexpected deleted elements");
        assert_eq!(8, mesh.num_vertices());
        assert_eq!(12, mesh.num_edges());
        assert_eq!(6, mesh.num_faces());
        assert_eq!(area, mesh.try_calc_area().expect("Cannot compute area"));
        assert_eq!(
            volume,
            mesh.try_calc_volume().expect("Cannot compute volume")
        );
    }

    #[test]
    fn t_icosahedron_revert_edge_collapse() {
        let mut mesh = PolyMeshF32::icosahedron(1.0).expect("Cannot make an icosahedron");
        let h = mesh
            .find_halfedge(0.into(), 1.into())
            .expect("Cannot find halfedge");
        let mut history = TopolHistory::default();
        let area = mesh.try_calc_area().expect("Cannot compute area");
        let volume = mesh.try_calc_volume().expect("Cannot compute volume");
        {
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
            let checkpt = history.check_point();
            history.commit_edge_collapse(&mesh, h, &vstatus, &hstatus, &estatus, &fstatus);
            let mut hcache = Vec::new();
            // Collapse and check.
            mesh.collapse_edge(
                h,
                &mut vstatus,
                &mut hstatus,
                &mut estatus,
                &mut fstatus,
                &mut hcache,
            );
            assert!(area > (mesh.try_calc_area().expect("Cannot compute area")));
            assert!(volume > mesh.try_calc_volume().expect("Cannot compute volume"));
            assert_eq!(
                11,
                mesh.vertices().filter(|&v| !vstatus[v].deleted()).count()
            );
            assert_eq!(27, mesh.edges().filter(|&e| !estatus[e].deleted()).count());
            assert_eq!(18, mesh.faces().filter(|&f| !fstatus[f].deleted()).count());
            // Revert and check again.
            assert!(history.restore(
                checkpt,
                &mut mesh,
                &mut vstatus,
                &mut hstatus,
                &mut estatus,
                &mut fstatus
            ));
        }
        // Nothing should be deleted.
        mesh.check_for_deleted()
            .expect("Unexpected deleted elements");
        assert_eq!(12, mesh.num_vertices());
        assert_eq!(30, mesh.num_edges());
        assert_eq!(20, mesh.num_faces());
        assert_eq!(area, mesh.try_calc_area().expect("Cannot compute area"));
        assert_eq!(
            volume,
            mesh.try_calc_volume().expect("Cannot compute volume")
        );
    }
}
