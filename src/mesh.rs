use std::cell::{Ref, RefMut};

use crate::{
    element::{EH, FH, HH, VH},
    error::Error,
    iterator,
    property::{EProperty, FProperty, HProperty, TPropData, VProperty},
    status::Status,
    topol::{TopolCache, Topology},
};

pub struct PolyMeshT<VecT: TPropData> {
    topol: Topology,
    points: VProperty<VecT>,
    cache: TopolCache,
}

impl<VecT: TPropData> Default for PolyMeshT<VecT> {
    fn default() -> Self {
        Self::new()
    }
}

impl<VecT: TPropData> PolyMeshT<VecT> {
    pub fn new() -> Self {
        let mut topol = Topology::new();
        let points = topol.create_vertex_prop();
        PolyMeshT {
            topol,
            points,
            cache: TopolCache::default(),
        }
    }

    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        let mut topol = Topology::with_capacity(nverts, nedges, nfaces);
        let points = VProperty::with_capacity(nverts, &mut topol.vprops);
        PolyMeshT {
            topol,
            points,
            cache: TopolCache::default(),
        }
    }

    pub fn create_vertex_prop<T: TPropData>(&mut self) -> VProperty<T> {
        self.topol.create_vertex_prop()
    }

    pub fn create_halfedge_prop<T: TPropData>(&mut self) -> HProperty<T> {
        self.topol.create_halfedge_prop()
    }

    pub fn create_edge_prop<T: TPropData>(&mut self) -> EProperty<T> {
        self.topol.create_edge_prop()
    }

    pub fn create_face_prop<T: TPropData>(&mut self) -> FProperty<T> {
        self.topol.create_face_prop()
    }

    pub fn reserve(&mut self, nverts: usize, nedges: usize, nfaces: usize) -> Result<(), Error> {
        self.topol.reserve(nverts, nedges, nfaces)
    }

    pub fn clear(&mut self) -> Result<(), Error> {
        self.topol.clear()
    }

    pub fn num_vertices(&self) -> usize {
        self.topol.num_vertices()
    }

    pub fn num_edges(&self) -> usize {
        self.topol.num_edges()
    }

    pub fn num_halfedges(&self) -> usize {
        self.topol.num_halfedges()
    }

    pub fn num_faces(&self) -> usize {
        self.topol.num_faces()
    }

    pub fn vertices(&self) -> impl Iterator<Item = VH> {
        self.topol.vertices()
    }

    pub fn halfedges(&self) -> impl Iterator<Item = HH> {
        self.topol.halfedges()
    }

    pub fn edges(&self) -> impl Iterator<Item = EH> {
        self.topol.edges()
    }

    pub fn faces(&self) -> impl Iterator<Item = FH> {
        self.topol.faces()
    }

    pub fn is_valid_vertex(&self, v: VH) -> bool {
        self.topol.is_valid_vertex(v)
    }

    pub fn is_valid_halfedge(&self, h: HH) -> bool {
        self.topol.is_valid_halfedge(h)
    }

    pub fn is_valid_edge(&self, e: EH) -> bool {
        self.topol.is_valid_edge(e)
    }

    pub fn is_valid_face(&self, f: FH) -> bool {
        self.topol.is_valid_face(f)
    }

    pub fn is_manifold_vertex(&self, v: VH) -> bool {
        self.topol.is_manifold_vertex(v)
    }

    pub fn is_boundary_edge(&self, e: EH) -> bool {
        self.topol.is_boundary_edge(e)
    }

    pub fn vertex_valence(&self, v: VH) -> usize {
        self.topol.vertex_valence(v)
    }

    pub fn face_valence(&self, f: FH) -> usize {
        self.topol.face_valence(f)
    }

    pub fn point(&self, vi: VH) -> Result<VecT, Error> {
        Ok(*(self.points.get(vi)?))
    }

    pub fn set_point(&mut self, v: VH, pos: VecT) -> Result<(), Error> {
        self.points.set(v, pos)
    }

    pub fn vertex_status(&self, v: VH) -> Result<Ref<'_, Status>, Error> {
        self.topol.vertex_status(v)
    }

    pub fn vertex_status_mut(&mut self, v: VH) -> Result<RefMut<'_, Status>, Error> {
        self.topol.vertex_status_mut(v)
    }

    pub fn halfedge_status(&self, h: HH) -> Result<Ref<'_, Status>, Error> {
        self.topol.halfedge_status(h)
    }

    pub fn halfedge_status_mut(&mut self, h: HH) -> Result<RefMut<'_, Status>, Error> {
        self.topol.halfedge_status_mut(h)
    }

    pub fn edge_status(&self, e: EH) -> Result<Ref<'_, Status>, Error> {
        self.topol.edge_status(e)
    }

    pub fn edge_status_mut(&mut self, e: EH) -> Result<RefMut<'_, Status>, Error> {
        self.topol.edge_status_mut(e)
    }

    pub fn face_status(&self, f: FH) -> Result<Ref<'_, Status>, Error> {
        self.topol.face_status(f)
    }

    pub fn face_status_mut(&mut self, f: FH) -> Result<RefMut<'_, Status>, Error> {
        self.topol.face_status_mut(f)
    }

    pub fn from_vertex(&self, h: HH) -> VH {
        self.topol.from_vertex(h)
    }

    pub fn halfedge_face(&self, h: HH) -> Option<FH> {
        self.topol.halfedge_face(h)
    }

    pub fn face_halfedge(&self, f: FH) -> HH {
        self.topol.face_halfedge(f)
    }

    pub fn halfedge_edge(&self, h: HH) -> EH {
        self.topol.halfedge_edge(h)
    }

    pub fn edge_halfedge(&self, e: EH, flag: bool) -> HH {
        self.topol.edge_halfedge(e, flag)
    }

    pub fn cw_rotated_halfedge(&self, h: HH) -> HH {
        self.topol.cw_rotated_halfedge(h)
    }

    pub fn ccw_rotated_halfedge(&self, h: HH) -> HH {
        self.topol.ccw_rotated_halfedge(h)
    }

    pub fn voh_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::voh_ccw_iter(&self.topol, v)
    }

    pub fn voh_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::voh_cw_iter(&self.topol, v)
    }

    pub fn vih_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::vih_ccw_iter(&self.topol, v)
    }

    pub fn vih_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::vih_cw_iter(&self.topol, v)
    }

    pub fn vf_ccw_iter(&self, v: VH) -> impl Iterator<Item = FH> + use<'_, VecT> {
        iterator::vf_ccw_iter(&self.topol, v)
    }

    pub fn vf_cw_iter(&self, v: VH) -> impl Iterator<Item = FH> + use<'_, VecT> {
        iterator::vf_cw_iter(&self.topol, v)
    }

    pub fn vv_ccw_iter(&self, v: VH) -> impl Iterator<Item = VH> + use<'_, VecT> {
        iterator::vv_ccw_iter(&self.topol, v)
    }

    pub fn vv_cw_iter(&self, v: VH) -> impl Iterator<Item = VH> + use<'_, VecT> {
        iterator::vv_cw_iter(&self.topol, v)
    }

    pub fn ve_ccw_iter(&self, v: VH) -> impl Iterator<Item = EH> + use<'_, VecT> {
        iterator::ve_ccw_iter(&self.topol, v)
    }

    pub fn ve_cw_iter(&self, v: VH) -> impl Iterator<Item = EH> + use<'_, VecT> {
        iterator::ve_cw_iter(&self.topol, v)
    }

    pub fn ev_iter(&self, e: EH) -> impl Iterator<Item = VH> + use<'_, VecT> {
        iterator::ev_iter(&self.topol, e)
    }

    pub fn eh_iter(&self, e: EH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::eh_iter(&self.topol, e)
    }

    pub fn ef_iter(&self, e: EH) -> impl Iterator<Item = FH> + use<'_, VecT> {
        iterator::ef_iter(&self.topol, e)
    }

    pub fn fh_ccw_iter(&self, f: FH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::fh_ccw_iter(&self.topol, f)
    }

    pub fn fh_cw_iter(&self, f: FH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::fh_cw_iter(&self.topol, f)
    }

    pub fn fv_ccw_iter(&self, f: FH) -> impl Iterator<Item = VH> + use<'_, VecT> {
        iterator::fv_ccw_iter(&self.topol, f)
    }

    pub fn fv_cw_iter(&self, f: FH) -> impl Iterator<Item = VH> + use<'_, VecT> {
        iterator::fv_cw_iter(&self.topol, f)
    }

    pub fn fe_ccw_iter(&self, f: FH) -> impl Iterator<Item = EH> + use<'_, VecT> {
        iterator::fe_ccw_iter(&self.topol, f)
    }

    pub fn fe_cw_iter(&self, f: FH) -> impl Iterator<Item = EH> + use<'_, VecT> {
        iterator::fe_cw_iter(&self.topol, f)
    }

    pub fn ff_ccw_iter(&self, f: FH) -> impl Iterator<Item = FH> + use<'_, VecT> {
        iterator::ff_ccw_iter(&self.topol, f)
    }

    pub fn ff_cw_iter(&self, f: FH) -> impl Iterator<Item = FH> + use<'_, VecT> {
        iterator::ff_cw_iter(&self.topol, f)
    }

    pub fn add_vertex(&mut self, pos: VecT) -> Result<VH, Error> {
        let vi = self.topol.add_vertex()?;
        self.points.set(vi, pos)?;
        Ok(vi)
    }

    pub fn add_face(&mut self, verts: &[VH]) -> Result<FH, Error> {
        self.topol.add_face(verts, &mut self.cache)
    }

    pub fn add_tri_face(&mut self, v0: VH, v1: VH, v2: VH) -> Result<FH, Error> {
        self.add_face(&[v0, v1, v2])
    }

    pub fn add_quad_face(&mut self, v0: VH, v1: VH, v2: VH, v3: VH) -> Result<FH, Error> {
        self.add_face(&[v0, v1, v2, v3])
    }

    pub fn delete_vertex(&mut self, delete_isolated_vertices: bool, v: VH) -> Result<(), Error> {
        self.topol
            .delete_vertex(v, delete_isolated_vertices, &mut self.cache)
    }

    pub fn delete_face(&mut self, f: FH, delete_isolated_vertices: bool) -> Result<(), Error> {
        self.topol.delete_face(
            f,
            delete_isolated_vertices,
            &mut self.cache.halfedges,
            &mut self.cache.edges,
            &mut self.cache.vertices,
        )
    }

    pub fn delete_edge(&mut self, e: EH, delete_isolated_vertices: bool) -> Result<(), Error> {
        self.topol.delete_edge(
            e,
            delete_isolated_vertices,
            &mut self.cache.halfedges,
            &mut self.cache.edges,
            &mut self.cache.vertices,
        )
    }

    pub fn garbage_collection(&mut self) -> Result<(), Error> {
        self.topol.garbage_collection(&mut self.cache)
    }
}
