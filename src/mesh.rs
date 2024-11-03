use std::ops::{Add, Div, Mul, Sub};

use crate::{
    element::{Handle, EH, FH, HH, VH},
    error::Error,
    iterator,
    property::{EProperty, FProperty, HProperty, TPropData, VProperty},
    status::Status,
    topol::{TopolCache, Topology},
    vector::{TScalar, TVec3},
};

pub struct PolyMeshT<VecT: TVec3> {
    topol: Topology,
    cache: TopolCache,
    points: VProperty<VecT>,
    vnormals: Option<VProperty<VecT>>,
    fnormals: Option<FProperty<VecT>>,
}

pub type PolyMeshF32 = PolyMeshT<glam::Vec3>;

impl<VecT: TVec3> Default for PolyMeshT<VecT> {
    fn default() -> Self {
        Self::new()
    }
}

impl<VecT: TVec3> PolyMeshT<VecT> {
    pub fn new() -> Self {
        let mut topol = Topology::new();
        let points = topol.new_vprop();
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        let mut topol = Topology::with_capacity(nverts, nedges, nfaces);
        let points = topol.new_vprop_with_capacity(nverts);
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    pub fn create_vertex_prop<T: TPropData>(&mut self) -> VProperty<T> {
        self.topol.new_vprop()
    }

    pub fn create_halfedge_prop<T: TPropData>(&mut self) -> HProperty<T> {
        self.topol.new_hprop()
    }

    pub fn create_edge_prop<T: TPropData>(&mut self) -> EProperty<T> {
        self.topol.new_eprop()
    }

    pub fn create_face_prop<T: TPropData>(&mut self) -> FProperty<T> {
        self.topol.new_fprop()
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

    pub(crate) fn topology(&self) -> &Topology {
        &self.topol
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
        self.points.get(vi)
    }

    pub fn vertex_status(&self, v: VH) -> Result<Status, Error> {
        self.topol.vertex_status(v)
    }

    pub fn halfedge_status(&self, h: HH) -> Result<Status, Error> {
        self.topol.halfedge_status(h)
    }

    pub fn edge_status(&self, e: EH) -> Result<Status, Error> {
        self.topol.edge_status(e)
    }

    pub fn face_status(&self, f: FH) -> Result<Status, Error> {
        self.topol.face_status(f)
    }

    pub fn set_point(&mut self, v: VH, pos: VecT) -> Result<(), Error> {
        self.points.set(v, pos)
    }

    pub fn points(&self) -> VProperty<VecT> {
        self.points.clone()
    }

    pub fn has_vertex_normals(&self) -> bool {
        self.vnormals.is_some()
    }

    pub fn vertex_normals(&self) -> Option<VProperty<VecT>> {
        self.vnormals.clone()
    }

    pub fn request_vertex_normals(&mut self) -> VProperty<VecT> {
        self.vnormals
            .get_or_insert_with(|| self.topol.new_vprop())
            .clone()
    }

    pub fn has_face_normals(&self) -> bool {
        self.fnormals.is_some()
    }

    pub fn face_normals(&self) -> Option<FProperty<VecT>> {
        self.fnormals.clone()
    }

    pub fn request_face_normals(&mut self) -> FProperty<VecT> {
        self.fnormals
            .get_or_insert_with(|| self.topol.new_fprop())
            .clone()
    }

    pub fn update_face_normals(&mut self) -> Result<FProperty<VecT>, Error>
    where
        VecT::Scalar: TScalar
            + Mul<Output = VecT::Scalar>
            + Add<Output = VecT::Scalar>
            + TPropData
            + Div<Output = VecT::Scalar>
            + PartialOrd,
        VecT: TVec3 + Add<Output = VecT> + Sub<Output = VecT>,
    {
        let mut fprop = self.request_face_normals();
        {
            let mut fnormals = fprop.try_borrow_mut()?;
            let fnormals: &mut [VecT] = &mut fnormals;
            let points = self.points();
            let points = points.try_borrow()?;
            for f in self.faces() {
                fnormals[f.index() as usize] = self.calc_face_normal(f, &points);
            }
        }
        Ok(fprop)
    }

    pub fn update_vertex_normals_fast(&mut self) -> Result<VProperty<VecT>, Error>
    where
        VecT::Scalar: TScalar
            + Mul<Output = VecT::Scalar>
            + Add<Output = VecT::Scalar>
            + TPropData
            + Div<Output = VecT::Scalar>
            + PartialOrd,
        VecT: TVec3 + Add<Output = VecT> + Sub<Output = VecT>,
    {
        let mut vprop = self.request_vertex_normals();
        {
            let fprop = match self.face_normals() {
                Some(prop) => prop,
                None => self.update_face_normals()?,
            };
            let fnormals = fprop.try_borrow()?;
            let mut vnormals = vprop.try_borrow_mut()?;
            let vnormals: &mut [VecT] = &mut vnormals;
            for v in self.vertices() {
                vnormals[v.index() as usize] = self.calc_vertex_normal_fast(v, &fnormals);
            }
        }
        Ok(vprop)
    }

    pub fn to_vertex(&self, h: HH) -> VH {
        self.topol.to_vertex(h)
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

    pub fn voh_ccw_circulator(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::voh_ccw_circulator(&self.topol, h)
    }

    pub fn voh_cw_circulator(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::voh_cw_circulator(&self.topol, h)
    }

    pub fn fh_ccw_circulator(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::fh_ccw_circulator(&self.topol, h)
    }

    pub fn fh_cw_circulator(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT> {
        iterator::fh_cw_circulator(&self.topol, h)
    }

    pub fn triangulated_vertices(&self) -> impl Iterator<Item = [VH; 3]> + use<'_, VecT> {
        self.topol.triangulated_vertices()
    }

    pub fn add_vertex(&mut self, pos: VecT) -> Result<VH, Error> {
        let vi = self.topol.add_vertex()?;
        self.points.set(vi, pos)?;
        Ok(vi)
    }

    pub fn add_vertices(&mut self, pos: &[VecT], dst: &mut [VH]) -> Result<(), Error> {
        if pos.len() != dst.len() {
            return Err(Error::MismatchedArrayLengths(pos.len(), dst.len()));
        }
        self.topol.add_vertices(dst)?;
        let verts: &[VH] = dst; // Make immutable.
        {
            let mut points = self.points.try_borrow_mut()?;
            let points: &mut [VecT] = &mut points;
            for (i, v) in verts.iter().enumerate() {
                points[v.index() as usize] = pos[i];
            }
        }
        Ok(())
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
