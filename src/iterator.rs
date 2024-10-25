use crate::topol::Topology;

struct OutgoingHalfedgeIter<'a, const CCW: bool> {
    topol: &'a Topology,
    hstart: Option<u32>,
    hcurrent: Option<u32>,
}

impl<'a> Iterator for OutgoingHalfedgeIter<'a, true> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = self
                    .topol
                    .opposite_halfedge(self.topol.prev_halfedge(current));
                self.hcurrent = match self.hstart {
                    Some(start) if start != next => Some(next),
                    _ => None,
                };
                Some(current)
            }
            None => None,
        }
    }
}

impl<'a> Iterator for OutgoingHalfedgeIter<'a, false> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = self
                    .topol
                    .next_halfedge(self.topol.opposite_halfedge(current));
                self.hcurrent = match self.hstart {
                    Some(start) if start != next => Some(next),
                    _ => None,
                };
                Some(current)
            }
            None => None,
        }
    }
}

struct FaceHalfedgeIter<'a, const CCW: bool> {
    topol: &'a Topology,
    hstart: u32,
    hcurrent: Option<u32>,
}

impl<'a> Iterator for FaceHalfedgeIter<'a, true> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = self.topol.next_halfedge(current);
                self.hcurrent = if next == self.hstart {
                    None
                } else {
                    Some(next)
                };
                Some(current)
            }
            None => None,
        }
    }
}

impl<'a> Iterator for FaceHalfedgeIter<'a, false> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.hcurrent {
            Some(current) => {
                let next = self.topol.prev_halfedge(current);
                self.hcurrent = if next == self.hstart {
                    None
                } else {
                    Some(next)
                };
                Some(current)
            }
            None => None,
        }
    }
}

pub(crate) fn voh_ccw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<true> {
        topol,
        hstart: h,
        hcurrent: h,
    }
}

pub(crate) fn voh_cw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<false> {
        topol,
        hstart: h,
        hcurrent: h,
    }
}

pub(crate) fn vih_ccw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<true> {
        topol,
        hstart: h,
        hcurrent: h,
    }
    .map(|h| topol.opposite_halfedge(h))
}

pub(crate) fn vih_cw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<false> {
        topol,
        hstart: h,
        hcurrent: h,
    }
    .map(|h| topol.opposite_halfedge(h))
}

pub(crate) fn vf_ccw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    voh_ccw_iter(topol, v).filter_map(|h| topol.halfedge_face(h))
}

pub(crate) fn vf_cw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    voh_cw_iter(topol, v).filter_map(|h| topol.halfedge_face(h))
}

pub(crate) fn vv_ccw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    voh_ccw_iter(topol, v).map(|h| topol.to_vertex(h))
}

pub(crate) fn vv_cw_iter<'a>(topol: &'a Topology, v: u32) -> impl Iterator<Item = u32> + use<'a> {
    voh_cw_iter(topol, v).map(|h| topol.to_vertex(h))
}

pub(crate) fn fh_ccw_iter<'a>(topol: &'a Topology, f: u32) -> impl Iterator<Item = u32> + use<'a> {
    let h = topol.face_halfedge(f);
    FaceHalfedgeIter::<true> {
        topol,
        hstart: h,
        hcurrent: Some(h),
    }
}

pub(crate) fn fh_cw_iter<'a>(topol: &'a Topology, f: u32) -> impl Iterator<Item = u32> + use<'a> {
    let h = topol.face_halfedge(f);
    FaceHalfedgeIter::<false> {
        topol,
        hstart: h,
        hcurrent: Some(h),
    }
}

pub(crate) fn fv_ccw_iter<'a>(topol: &'a Topology, f: u32) -> impl Iterator<Item = u32> + use<'a> {
    fh_ccw_iter(topol, f).map(|h| topol.to_vertex(h))
}

pub(crate) fn fv_cw_iter<'a>(topol: &'a Topology, f: u32) -> impl Iterator<Item = u32> + use<'a> {
    fh_cw_iter(topol, f).map(|h| topol.to_vertex(h))
}

#[cfg(test)]
mod test {
    use crate::{
        iterator::{fv_ccw_iter, fv_cw_iter, vf_ccw_iter, vf_cw_iter, vv_ccw_iter, vv_cw_iter},
        topol::{TopolCache, Topology},
    };

    fn quad_box() -> Topology {
        let mut topol = Topology::with_capacity(8, 12, 6);
        let verts: Vec<_> = (0..8)
            .map(|_| topol.add_vertex().expect("Unable to add a vertex"))
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
                .add_face(indices, &mut cache)
                .expect("Unable to add a face")
        })
        .collect();
        assert_eq!(faces, (0..6).collect::<Vec<_>>());
        assert_eq!(topol.num_vertices(), 8);
        assert_eq!(topol.num_halfedges(), 24);
        assert_eq!(topol.num_edges(), 12);
        assert_eq!(topol.num_faces(), 6);
        topol
    }

    #[test]
    fn t_box_vf_ccw_iter() {
        let qbox = quad_box();
        assert_eq!(&vf_ccw_iter(&qbox, 0).collect::<Vec<_>>(), &[4, 0, 1]);
        assert_eq!(&vf_ccw_iter(&qbox, 1).collect::<Vec<_>>(), &[2, 1, 0]);
        assert_eq!(&vf_ccw_iter(&qbox, 2).collect::<Vec<_>>(), &[3, 2, 0]);
        assert_eq!(&vf_ccw_iter(&qbox, 3).collect::<Vec<_>>(), &[4, 3, 0]);
        assert_eq!(&vf_ccw_iter(&qbox, 4).collect::<Vec<_>>(), &[5, 4, 1]);
        assert_eq!(&vf_ccw_iter(&qbox, 5).collect::<Vec<_>>(), &[5, 1, 2]);
        assert_eq!(&vf_ccw_iter(&qbox, 6).collect::<Vec<_>>(), &[5, 2, 3]);
        assert_eq!(&vf_ccw_iter(&qbox, 7).collect::<Vec<_>>(), &[5, 3, 4]);
    }

    #[test]
    fn t_box_vf_cw_iter() {
        let qbox = quad_box();
        assert_eq!(&vf_cw_iter(&qbox, 0).collect::<Vec<_>>(), &[4, 1, 0]);
        assert_eq!(&vf_cw_iter(&qbox, 1).collect::<Vec<_>>(), &[2, 0, 1]);
        assert_eq!(&vf_cw_iter(&qbox, 2).collect::<Vec<_>>(), &[3, 0, 2]);
        assert_eq!(&vf_cw_iter(&qbox, 3).collect::<Vec<_>>(), &[4, 0, 3]);
        assert_eq!(&vf_cw_iter(&qbox, 4).collect::<Vec<_>>(), &[5, 1, 4]);
        assert_eq!(&vf_cw_iter(&qbox, 5).collect::<Vec<_>>(), &[5, 2, 1]);
        assert_eq!(&vf_cw_iter(&qbox, 6).collect::<Vec<_>>(), &[5, 3, 2]);
        assert_eq!(&vf_cw_iter(&qbox, 7).collect::<Vec<_>>(), &[5, 4, 3]);
    }

    #[test]
    fn t_box_vv_ccw_iter() {
        let qbox = quad_box();
        assert_eq!(&vv_ccw_iter(&qbox, 0).collect::<Vec<_>>(), &[4, 3, 1]);
        assert_eq!(&vv_ccw_iter(&qbox, 1).collect::<Vec<_>>(), &[2, 5, 0]);
        assert_eq!(&vv_ccw_iter(&qbox, 2).collect::<Vec<_>>(), &[3, 6, 1]);
        assert_eq!(&vv_ccw_iter(&qbox, 3).collect::<Vec<_>>(), &[0, 7, 2]);
        assert_eq!(&vv_ccw_iter(&qbox, 4).collect::<Vec<_>>(), &[5, 7, 0]);
        assert_eq!(&vv_ccw_iter(&qbox, 5).collect::<Vec<_>>(), &[6, 4, 1]);
        assert_eq!(&vv_ccw_iter(&qbox, 6).collect::<Vec<_>>(), &[7, 5, 2]);
        assert_eq!(&vv_ccw_iter(&qbox, 7).collect::<Vec<_>>(), &[4, 6, 3]);
    }

    #[test]
    fn t_box_vv_cw_iter() {
        let qbox = quad_box();
        assert_eq!(&vv_cw_iter(&qbox, 0).collect::<Vec<_>>(), &[4, 1, 3]);
        assert_eq!(&vv_cw_iter(&qbox, 1).collect::<Vec<_>>(), &[2, 0, 5]);
        assert_eq!(&vv_cw_iter(&qbox, 2).collect::<Vec<_>>(), &[3, 1, 6]);
        assert_eq!(&vv_cw_iter(&qbox, 3).collect::<Vec<_>>(), &[0, 2, 7]);
        assert_eq!(&vv_cw_iter(&qbox, 4).collect::<Vec<_>>(), &[5, 0, 7]);
        assert_eq!(&vv_cw_iter(&qbox, 5).collect::<Vec<_>>(), &[6, 1, 4]);
        assert_eq!(&vv_cw_iter(&qbox, 6).collect::<Vec<_>>(), &[7, 2, 5]);
        assert_eq!(&vv_cw_iter(&qbox, 7).collect::<Vec<_>>(), &[4, 3, 6]);
    }

    #[test]
    fn t_box_fv_ccw_iter() {
        let qbox = quad_box();
        assert_eq!(&fv_ccw_iter(&qbox, 0).collect::<Vec<_>>(), &[0, 3, 2, 1]);
        assert_eq!(&fv_ccw_iter(&qbox, 1).collect::<Vec<_>>(), &[0, 1, 5, 4]);
        assert_eq!(&fv_ccw_iter(&qbox, 2).collect::<Vec<_>>(), &[1, 2, 6, 5]);
        assert_eq!(&fv_ccw_iter(&qbox, 3).collect::<Vec<_>>(), &[2, 3, 7, 6]);
        assert_eq!(&fv_ccw_iter(&qbox, 4).collect::<Vec<_>>(), &[3, 0, 4, 7]);
        assert_eq!(&fv_ccw_iter(&qbox, 5).collect::<Vec<_>>(), &[4, 5, 6, 7]);
    }

    #[test]
    fn t_box_fv_cw_iter() {
        let qbox = quad_box();
        assert_eq!(&fv_cw_iter(&qbox, 0).collect::<Vec<_>>(), &[0, 1, 2, 3]);
        assert_eq!(&fv_cw_iter(&qbox, 1).collect::<Vec<_>>(), &[0, 4, 5, 1]);
        assert_eq!(&fv_cw_iter(&qbox, 2).collect::<Vec<_>>(), &[1, 5, 6, 2]);
        assert_eq!(&fv_cw_iter(&qbox, 3).collect::<Vec<_>>(), &[2, 6, 7, 3]);
        assert_eq!(&fv_cw_iter(&qbox, 4).collect::<Vec<_>>(), &[3, 7, 4, 0]);
        assert_eq!(&fv_cw_iter(&qbox, 5).collect::<Vec<_>>(), &[4, 7, 6, 5]);
    }
}
