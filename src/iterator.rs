use crate::topol::{Topology, EH, FH, HH, VH};

struct OutgoingHalfedgeIter<'a, const CCW: bool> {
    topol: &'a Topology,
    hstart: Option<HH>,
    hcurrent: Option<HH>,
}

impl<'a> Iterator for OutgoingHalfedgeIter<'a, true> {
    type Item = HH;

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
    type Item = HH;

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
    hstart: HH,
    hcurrent: Option<HH>,
}

impl<'a> Iterator for FaceHalfedgeIter<'a, true> {
    type Item = HH;

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
    type Item = HH;

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

pub(crate) fn vv_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = VH> + use<'_> {
    voh_ccw_iter(topol, v).map(|h| topol.to_vertex(h))
}

pub(crate) fn vv_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = VH> + use<'_> {
    voh_cw_iter(topol, v).map(|h| topol.to_vertex(h))
}

pub(crate) fn vih_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<true> {
        topol,
        hstart: h,
        hcurrent: h,
    }
    .map(|h| topol.opposite_halfedge(h))
}

pub(crate) fn vih_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<false> {
        topol,
        hstart: h,
        hcurrent: h,
    }
    .map(|h| topol.opposite_halfedge(h))
}

pub(crate) fn voh_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<true> {
        topol,
        hstart: h,
        hcurrent: h,
    }
}

pub(crate) fn voh_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = HH> + use<'_> {
    let h = topol.vertex_halfedge(v);
    OutgoingHalfedgeIter::<false> {
        topol,
        hstart: h,
        hcurrent: h,
    }
}

pub(crate) fn ve_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = EH> + use<'_> {
    voh_ccw_iter(topol, v).map(|h| topol.halfedge_edge(h))
}

pub(crate) fn ve_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = EH> + use<'_> {
    voh_cw_iter(topol, v).map(|h| topol.halfedge_edge(h))
}

pub(crate) fn vf_ccw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = FH> + use<'_> {
    voh_ccw_iter(topol, v).filter_map(|h| topol.halfedge_face(h))
}

pub(crate) fn vf_cw_iter(topol: &Topology, v: VH) -> impl Iterator<Item = FH> + use<'_> {
    voh_cw_iter(topol, v).filter_map(|h| topol.halfedge_face(h))
}

pub(crate) fn ev_iter(topol: &Topology, e: EH) -> impl Iterator<Item = VH> + use<'_> {
    eh_iter(topol, e).map(|h| topol.to_vertex(h))
}

pub(crate) fn eh_iter(topol: &Topology, e: EH) -> impl Iterator<Item = HH> + use<'_> {
    [false, true]
        .iter()
        .map(move |flag| topol.edge_halfedge(e, *flag))
}

pub(crate) fn ef_iter(topol: &Topology, e: EH) -> impl Iterator<Item = FH> + use<'_> {
    eh_iter(topol, e).filter_map(|h| topol.halfedge_face(h))
}

pub(crate) fn fv_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = VH> + use<'_> {
    fh_ccw_iter(topol, f).map(|h| topol.to_vertex(h))
}

pub(crate) fn fv_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = VH> + use<'_> {
    fh_cw_iter(topol, f).map(|h| topol.to_vertex(h))
}

pub(crate) fn fh_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = HH> + use<'_> {
    let h = topol.face_halfedge(f);
    FaceHalfedgeIter::<true> {
        topol,
        hstart: h,
        hcurrent: Some(h),
    }
}

pub(crate) fn fh_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = HH> + use<'_> {
    let h = topol.face_halfedge(f);
    FaceHalfedgeIter::<false> {
        topol,
        hstart: h,
        hcurrent: Some(h),
    }
}

pub(crate) fn fe_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = EH> + use<'_> {
    fh_ccw_iter(topol, f).map(|h| topol.halfedge_edge(h))
}

pub(crate) fn fe_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = EH> + use<'_> {
    fh_cw_iter(topol, f).map(|h| topol.halfedge_edge(h))
}

pub(crate) fn ff_ccw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = FH> + use<'_> {
    fh_ccw_iter(topol, f).filter_map(|h| topol.halfedge_face(topol.opposite_halfedge(h)))
}

pub(crate) fn ff_cw_iter(topol: &Topology, f: FH) -> impl Iterator<Item = FH> + use<'_> {
    fh_cw_iter(topol, f).filter_map(|h| topol.halfedge_face(topol.opposite_halfedge(h)))
}

#[cfg(test)]
mod test {
    use crate::{
        iterator::{
            ff_ccw_iter, ff_cw_iter, fv_ccw_iter, fv_cw_iter, vf_ccw_iter, vf_cw_iter,
            vih_ccw_iter, vih_cw_iter, voh_ccw_iter, voh_cw_iter, vv_ccw_iter, vv_cw_iter,
        },
        topol::{Handle, TopolCache, Topology},
    };

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
                .add_face(
                    &indices.iter().map(|idx| (*idx).into()).collect::<Vec<_>>(),
                    &mut cache,
                )
                .expect("Unable to add a face")
        })
        .collect();
        assert_eq!(faces, (0u32..6).map(|i| i.into()).collect::<Vec<_>>());
        assert_eq!(topol.num_vertices(), 8);
        assert_eq!(topol.num_halfedges(), 24);
        assert_eq!(topol.num_edges(), 12);
        assert_eq!(topol.num_faces(), 6);
        topol
    }

    #[test]
    fn t_box_vv_ccw_iter() {
        let qbox = quad_box();
        for (vi, vis) in [
            (0u32, [4u32, 3, 1]),
            (1u32, [2u32, 5, 0]),
            (2u32, [3u32, 6, 1]),
            (3u32, [0u32, 7, 2]),
            (4u32, [5u32, 7, 0]),
            (5u32, [6u32, 4, 1]),
            (6u32, [7u32, 5, 2]),
            (7u32, [4u32, 6, 3]),
        ] {
            assert!(vv_ccw_iter(&qbox, vi.into()).eq(vis.iter().map(|i| (*i).into())));
        }
    }

    #[test]
    fn t_box_vv_cw_iter() {
        let qbox = quad_box();
        for (vi, fis) in [
            (0u32, [4, 1, 3]),
            (1u32, [2, 0, 5]),
            (2u32, [3, 1, 6]),
            (3u32, [0, 2, 7]),
            (4u32, [5, 0, 7]),
            (5u32, [6, 1, 4]),
            (6u32, [7, 2, 5]),
            (7u32, [4, 3, 6]),
        ] {
            assert!(vv_cw_iter(&qbox, vi.into()).eq(fis.iter().map(|i| (*i).into())));
        }
    }

    #[test]
    fn t_box_vih_iter() {
        let qbox = quad_box();
        for v in qbox.vertex_iter() {
            assert!(
                vih_ccw_iter(&qbox, v).all(|h| qbox.to_vertex(h) == v && qbox.from_vertex(h) != v)
            );
            assert!(
                vih_cw_iter(&qbox, v).all(|h| qbox.to_vertex(h) == v && qbox.from_vertex(h) != v)
            );
        }
    }

    #[test]
    fn t_box_voh_iter() {
        let qbox = quad_box();
        for v in qbox.vertex_iter() {
            assert!(
                voh_ccw_iter(&qbox, v).all(|h| qbox.from_vertex(h) == v && qbox.to_vertex(h) != v)
            );
            assert!(
                voh_cw_iter(&qbox, v).all(|h| qbox.from_vertex(h) == v && qbox.to_vertex(h) != v)
            );
        }
    }

    #[test]
    fn t_box_vf_ccw_iter() {
        let qbox = quad_box();
        for (vi, fis) in [
            (0u32, [4u32, 0, 1]),
            (1u32, [2u32, 1, 0]),
            (2u32, [3u32, 2, 0]),
            (3u32, [4u32, 3, 0]),
            (4u32, [5u32, 4, 1]),
            (5u32, [5u32, 1, 2]),
            (6u32, [5u32, 2, 3]),
            (7u32, [5u32, 3, 4]),
        ] {
            assert!(vf_ccw_iter(&qbox, vi.into()).eq(fis.iter().map(|i| (*i).into())));
        }
    }

    #[test]
    fn t_box_vf_cw_iter() {
        let qbox = quad_box();
        for (vi, fis) in [
            (0u32, [4u32, 1, 0]),
            (1, [2, 0, 1]),
            (2, [3, 0, 2]),
            (3, [4, 0, 3]),
            (4, [5, 1, 4]),
            (5, [5, 2, 1]),
            (6, [5, 3, 2]),
            (7, [5, 4, 3]),
        ] {
            assert!(vf_cw_iter(&qbox, vi.into()).eq(fis.iter().map(|i| (*i).into())));
        }
    }

    #[test]
    fn t_box_fv_ccw_iter() {
        let qbox = quad_box();
        for (fi, vis) in [
            (0u32, [0, 3, 2, 1]),
            (1u32, [0, 1, 5, 4]),
            (2u32, [1, 2, 6, 5]),
            (3u32, [2, 3, 7, 6]),
            (4u32, [3, 0, 4, 7]),
            (5u32, [4, 5, 6, 7]),
        ] {
            assert!(fv_ccw_iter(&qbox, fi.into()).eq(vis.iter().map(|i| (*i).into())));
        }
    }

    #[test]
    fn t_box_fv_cw_iter() {
        let qbox = quad_box();
        for (fi, vis) in [
            (0u32, [0, 1, 2, 3]),
            (1u32, [0, 4, 5, 1]),
            (2u32, [1, 5, 6, 2]),
            (3u32, [2, 6, 7, 3]),
            (4u32, [3, 7, 4, 0]),
            (5u32, [4, 7, 6, 5]),
        ] {
            assert!(fv_cw_iter(&qbox, fi.into()).eq(vis.iter().map(|i| (*i).into())));
        }
    }

    #[test]
    fn t_box_ff_ccw_iter() {
        let qbox = quad_box();
        for (fi, fis) in [
            (0u32, [1, 4, 3, 2]),
            (1u32, [4, 0, 2, 5]),
            (2u32, [1, 0, 3, 5]),
            (3u32, [2, 0, 4, 5]),
            (4u32, [3, 0, 1, 5]),
            (5u32, [4, 1, 2, 3]),
        ] {
            assert_eq!(
                ff_ccw_iter(&qbox, fi.into())
                    .map(|f| f.index())
                    .collect::<Vec<_>>(),
                &fis
            );
        }
    }

    #[test]
    fn t_box_ff_cw_iter() {
        let qbox = quad_box();
        for (fi, fis) in [
            (0u32, [1, 2, 3, 4]),
            (1u32, [4, 5, 2, 0]),
            (2u32, [1, 5, 3, 0]),
            (3u32, [2, 5, 4, 0]),
            (4u32, [3, 5, 1, 0]),
            (5u32, [4, 3, 2, 1]),
        ] {
            assert_eq!(
                ff_cw_iter(&qbox, fi.into())
                    .map(|f| f.index())
                    .collect::<Vec<_>>(),
                &fis
            );
        }
    }
}
