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
    cache: Vec<Element>,
    vertices: Vec<VH>,
    halfedges: Vec<HH>,
    edges: Vec<EH>,
    faces: Vec<FH>,
}

impl TopolHistory {
    fn commit_vertices(&mut self, mesh: &Topology, vstatus: &VPropBuf<Status>) {
        self.vertices.sort();
        self.vertices.dedup();
        self.cache.extend(
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
        self.cache.extend(
            self.halfedges
                .drain(..)
                .map(|h| Element::Halfedge(h, *mesh.halfedge(h), hstatus[h])),
        );
        self.cache.extend(
            self.edges
                .drain(..)
                .map(|e| Element::EdgeStatus(e, estatus[e])),
        );
    }

    fn commit_faces(&mut self, mesh: &Topology, fstatus: &FPropBuf<Status>) {
        self.faces.sort();
        self.faces.dedup();
        self.cache.extend(
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

    pub fn record_edge_collapse(
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

    fn commit(
        &mut self,
        mesh: &Topology,
        vstatus: &VPropBuf<Status>,
        hstatus: &HPropBuf<Status>,
        estatus: &EPropBuf<Status>,
        fstatus: &FPropBuf<Status>,
    ) {
        self.commit_edges(mesh, hstatus, estatus);
        self.commit_faces(mesh, fstatus);
        self.commit_vertices(mesh, vstatus);
    }

    pub fn check_point(&self) -> usize {
        self.cache.len()
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
        if check_point >= self.cache.len() {
            return false;
        }
        let mesh = mesh.topology_mut();
        for op in self.cache[check_point..].iter().rev() {
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
        self.cache.truncate(check_point);
        true
    }
}
