use std::ops::Range;

use crate::{
    error::Error,
    iterator,
    property::Property,
    topol::{TopolCache, Topology},
};

pub struct Mesh {
    topol: Topology,
    cache: TopolCache,
    points: Property<glam::Vec3>,
}

impl Default for Mesh {
    fn default() -> Self {
        Self::new()
    }
}

impl Mesh {
    pub fn new() -> Self {
        let mut topol = Topology::new();
        let points = topol.create_vertex_prop::<glam::Vec3>();
        Mesh {
            topol,
            points,
            cache: TopolCache::default(),
        }
    }

    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        let mut topol = Topology::with_capacity(nverts, nedges, nfaces);
        let points = topol.create_vertex_prop::<glam::Vec3>();
        Mesh {
            topol,
            points,
            cache: TopolCache::default(),
        }
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

    pub fn vertex_iter(&self) -> Range<u32> {
        self.topol.vertex_iter()
    }

    pub fn halfedge_iter(&self) -> Range<u32> {
        self.topol.halfedge_iter()
    }

    pub fn edge_iter(&self) -> Range<u32> {
        self.topol.edge_iter()
    }

    pub fn face_iter(&self) -> Range<u32> {
        self.topol.face_iter()
    }

    pub fn is_boundary_edge(&self, e: u32) -> bool {
        self.topol.is_boundary_edge(e)
    }

    pub fn point(&self, vi: u32) -> Result<glam::Vec3, Error> {
        self.points.get(vi)
    }

    pub fn from_vertex(&self, h: u32) -> u32 {
        self.topol.from_vertex(h)
    }

    pub fn halfedge_face(&self, h: u32) -> Option<u32> {
        self.topol.halfedge_face(h)
    }

    pub fn face_halfedge(&self, f: u32) -> u32 {
        self.topol.face_halfedge(f)
    }

    pub const fn halfedge_edge(&self, h: u32) -> u32 {
        self.topol.halfedge_edge(h)
    }

    pub const fn edge_halfedge(&self, e: u32, flag: bool) -> u32 {
        self.topol.edge_halfedge(e, flag)
    }

    pub fn cw_rotated_halfedge(&self, h: u32) -> u32 {
        self.topol.cw_rotated_halfedge(h)
    }

    pub fn ccw_rotated_halfedge(&self, h: u32) -> u32 {
        self.topol.ccw_rotated_halfedge(h)
    }

    pub fn voh_ccw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::voh_ccw_iter(&self.topol, v)
    }

    pub fn voh_cw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::voh_cw_iter(&self.topol, v)
    }

    pub fn vih_ccw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::vih_ccw_iter(&self.topol, v)
    }

    pub fn vih_cw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::vih_cw_iter(&self.topol, v)
    }

    pub fn vf_ccw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::vf_ccw_iter(&self.topol, v)
    }

    pub fn vf_cw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::vf_cw_iter(&self.topol, v)
    }

    pub fn vv_ccw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::vv_ccw_iter(&self.topol, v)
    }

    pub fn vv_cw_iter<'a>(&'a self, v: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::vv_cw_iter(&self.topol, v)
    }

    pub fn fh_ccw_iter<'a>(&'a self, f: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::fh_ccw_iter(&self.topol, f)
    }

    pub fn fh_cw_iter<'a>(&'a self, f: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::fh_cw_iter(&self.topol, f)
    }

    pub fn fv_ccw_iter<'a>(&'a self, f: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::fv_ccw_iter(&self.topol, f)
    }

    pub fn fv_cw_iter<'a>(&'a self, f: u32) -> impl Iterator<Item = u32> + use<'a> {
        iterator::fv_cw_iter(&self.topol, f)
    }

    pub fn add_vertex(&mut self, pos: glam::Vec3) -> Result<u32, Error> {
        let vi = self.topol.add_vertex()?;
        self.points.set(vi, pos)?;
        Ok(vi)
    }

    pub fn add_face(&mut self, verts: &[u32]) -> Result<u32, Error> {
        self.topol.add_face(verts, &mut self.cache)
    }

    pub fn add_tri_face(&mut self, v0: u32, v1: u32, v2: u32) -> Result<u32, Error> {
        self.add_face(&[v0, v1, v2])
    }

    pub fn add_quad_face(&mut self, v0: u32, v1: u32, v2: u32, v3: u32) -> Result<u32, Error> {
        self.add_face(&[v0, v1, v2, v3])
    }
}
