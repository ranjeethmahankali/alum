use crate::{
    element::VH,
    error::Error,
    mesh::PolyMeshT,
    property::TPropData,
    vector::{FromFloat, TVec},
};

impl<VecT> PolyMeshT<VecT, 3>
where
    VecT: TVec<3>,
{
    /// Makes a box with the following topology, spanning from the min point to
    /// the max point.
    ///
    ///  ```text
    ///       7-----------6
    ///      /|          /|
    ///     / |         / |
    ///    4-----------5  |
    ///    |  |        |  |
    ///    |  3--------|--2
    ///    | /         | /
    ///    |/          |/
    ///    0-----------1
    ///  ```
    pub fn quad_box(min: VecT, max: VecT) -> Result<Self, Error> {
        const BOX_POS: [(bool, bool, bool); 8] = [
            (false, false, false),
            (true, false, false),
            (true, true, false),
            (false, true, false),
            (false, false, true),
            (true, false, true),
            (true, true, true),
            (false, true, true),
        ];
        const BOX_IDX: [(usize, usize, usize, usize); 6] = [
            (0, 3, 2, 1),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (3, 0, 4, 7),
            (4, 5, 6, 7),
        ];
        let mut qbox = Self::with_capacity(8, 12, 6);
        let verts = {
            let mut pos: [VecT; 8] = [VecT::zero(); 8];
            for (i, &(xf, yf, zf)) in BOX_POS.iter().enumerate() {
                pos[i] = VecT::new([
                    if xf { max.coord(0) } else { min.coord(0) },
                    if yf { max.coord(1) } else { min.coord(1) },
                    if zf { max.coord(2) } else { min.coord(2) },
                ]);
            }
            let mut verts: [VH; 8] = [0.into(); 8];
            qbox.add_vertices(&pos, &mut verts)?;
            verts
        };
        // TODO: Find a more efficient way of adding a large number of faces
        // that doesn't require repeatedly borrowing the property buffers.
        for (a, b, c, d) in BOX_IDX {
            qbox.add_quad_face(verts[a], verts[b], verts[c], verts[d])?;
        }
        Ok(qbox)
    }

    /// Create a mesh representing a box with quadrilateral faces, of size 1,
    /// spanning from the origin to (1, 1, 1).
    pub fn unit_box() -> Result<Self, Error>
    where
        VecT::Scalar: FromFloat + TPropData,
    {
        Self::quad_box(
            VecT::new([VecT::Scalar::from_f64(0.); 3]),
            VecT::new([VecT::Scalar::from_f64(1.); 3]),
        )
    }
}

#[cfg(test)]
mod test {
    use crate::mesh::PolyMeshF32;

    #[test]
    fn t_quad_box() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a quad box mesh");
        assert_eq!(qbox.num_vertices(), 8);
        assert_eq!(qbox.num_halfedges(), 24);
        assert_eq!(qbox.num_edges(), 12);
        assert_eq!(qbox.num_faces(), 6);
        for v in qbox.vertices() {
            assert_eq!(qbox.vf_ccw_iter(v).count(), 3);
        }
    }
}
