use std::ops::{Mul, Neg};

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

/// Platonic solids.
impl<VecT> PolyMeshT<VecT, 3>
where
    VecT: TVec<3>,
    VecT::Scalar: TPropData + FromFloat + Mul<Output = VecT::Scalar> + Neg<Output = VecT::Scalar>,
{
    pub fn tetrahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(4, 6, 4);
        let a = radius * VecT::Scalar::from_f64(1.0f64 / 3.0);
        let b = radius * VecT::Scalar::from_f64((8.0 / 9.0f64).sqrt());
        let c = radius * VecT::Scalar::from_f64((2.0 / 9.0f64).sqrt());
        let d = radius * VecT::Scalar::from_f64((2.0 / 3.0f64).sqrt());
        let verts = {
            let mut verts: [VH; 4] = [0.into(); 4];
            mesh.add_vertices(
                &[
                    VecT::new([
                        VecT::Scalar::from_f64(0.0),
                        VecT::Scalar::from_f64(0.0),
                        VecT::Scalar::from_f64(1.0),
                    ]),
                    VecT::new([-c, d, -a]),
                    VecT::new([-c, -d, -a]),
                    VecT::new([b, VecT::Scalar::from_f64(0.0), -a]),
                ],
                &mut verts,
            )?;
            verts
        };
        mesh.add_tri_face(verts[0], verts[1], verts[2])?;
        mesh.add_tri_face(verts[0], verts[2], verts[3])?;
        mesh.add_tri_face(verts[0], verts[3], verts[1])?;
        mesh.add_tri_face(verts[3], verts[2], verts[1])?;
        Ok(mesh)
    }

    pub fn hexahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let a = radius * VecT::Scalar::from_f64(1.0f64 / 3.0f64.sqrt());
        let mut mesh = Self::with_capacity(8, 12, 6);
        let verts = {
            let mut verts: [VH; 8] = [0.into(); 8];
            mesh.add_vertices(
                &[
                    VecT::new([-a, -a, -a]),
                    VecT::new([a, -a, -a]),
                    VecT::new([a, a, -a]),
                    VecT::new([-a, a, -a]),
                    VecT::new([-a, -a, a]),
                    VecT::new([a, -a, a]),
                    VecT::new([a, a, a]),
                    VecT::new([-a, a, a]),
                ],
                &mut verts,
            )?;
            verts
        };
        mesh.add_quad_face(verts[3], verts[2], verts[1], verts[0])?;
        mesh.add_quad_face(verts[2], verts[6], verts[5], verts[1])?;
        mesh.add_quad_face(verts[5], verts[6], verts[7], verts[4])?;
        mesh.add_quad_face(verts[0], verts[4], verts[7], verts[3])?;
        mesh.add_quad_face(verts[3], verts[7], verts[6], verts[2])?;
        mesh.add_quad_face(verts[1], verts[5], verts[4], verts[0])?;
        Ok(mesh)
    }

    pub fn octahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(6, 12, 8);
        let verts = {
            let mut verts: [VH; 6] = [0.into(); 6];
            let zero = VecT::Scalar::from_f64(0.);
            mesh.add_vertices(
                &[
                    VecT::new([radius, zero, zero]),
                    VecT::new([zero, radius, zero]),
                    VecT::new([-radius, zero, zero]),
                    VecT::new([zero, -radius, zero]),
                    VecT::new([zero, zero, radius]),
                    VecT::new([zero, zero, -radius]),
                ],
                &mut verts,
            )?;
            verts
        };
        mesh.add_tri_face(verts[0], verts[4], verts[3])?;
        mesh.add_tri_face(verts[1], verts[4], verts[0])?;
        mesh.add_tri_face(verts[2], verts[4], verts[1])?;
        mesh.add_tri_face(verts[3], verts[4], verts[2])?;
        mesh.add_tri_face(verts[3], verts[5], verts[0])?;
        mesh.add_tri_face(verts[0], verts[5], verts[1])?;
        mesh.add_tri_face(verts[1], verts[5], verts[2])?;
        mesh.add_tri_face(verts[2], verts[5], verts[3])?;
        Ok(mesh)
    }

    pub fn icosahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let a = radius * VecT::Scalar::from_f64(1.0f64);
        let b = radius * VecT::Scalar::from_f64(1.0 / ((1.0f64 + 5.0f64.sqrt()) * 0.5f64));
        let zero = VecT::Scalar::from_f64(0.0f64);
        let mut mesh = Self::with_capacity(12, 30, 20);
        let verts = {
            let mut verts: [VH; 12] = [0.into(); 12];
            mesh.add_vertices(
                &[
                    VecT::new([zero, b, -a]),
                    VecT::new([b, a, zero]),
                    VecT::new([-b, a, zero]),
                    VecT::new([zero, b, a]),
                    VecT::new([zero, -b, a]),
                    VecT::new([-a, zero, b]),
                    VecT::new([zero, -b, -a]),
                    VecT::new([a, zero, -b]),
                    VecT::new([a, zero, b]),
                    VecT::new([-a, zero, -b]),
                    VecT::new([b, -a, zero]),
                    VecT::new([-b, -a, zero]),
                ],
                &mut verts,
            )?;
            verts
        };
        mesh.add_tri_face(verts[2], verts[1], verts[0])?;
        mesh.add_tri_face(verts[1], verts[2], verts[3])?;
        mesh.add_tri_face(verts[5], verts[4], verts[3])?;
        mesh.add_tri_face(verts[4], verts[9], verts[3])?;
        mesh.add_tri_face(verts[7], verts[6], verts[0])?;
        mesh.add_tri_face(verts[6], verts[9], verts[0])?;
        mesh.add_tri_face(verts[11], verts[10], verts[4])?;
        mesh.add_tri_face(verts[10], verts[11], verts[6])?;
        mesh.add_tri_face(verts[9], verts[5], verts[2])?;
        mesh.add_tri_face(verts[5], verts[9], verts[11])?;
        mesh.add_tri_face(verts[9], verts[7], verts[1])?;
        mesh.add_tri_face(verts[7], verts[9], verts[10])?;
        mesh.add_tri_face(verts[2], verts[5], verts[3])?;
        mesh.add_tri_face(verts[9], verts[1], verts[3])?;
        mesh.add_tri_face(verts[9], verts[2], verts[0])?;
        mesh.add_tri_face(verts[1], verts[7], verts[0])?;
        mesh.add_tri_face(verts[11], verts[9], verts[6])?;
        mesh.add_tri_face(verts[7], verts[10], verts[6])?;
        mesh.add_tri_face(verts[5], verts[11], verts[4])?;
        mesh.add_tri_face(verts[10], verts[9], verts[4])?;
        Ok(mesh)
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
