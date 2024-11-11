use std::ops::{Add, Div, Mul, Neg};

use crate::{
    element::{Handle, VH},
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
        const BOX_IDX: [(u32, u32, u32, u32); 6] = [
            (0, 3, 2, 1),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (3, 0, 4, 7),
            (4, 5, 6, 7),
        ];
        let mut qbox = Self::with_capacity(8, 12, 6);
        let _verts = {
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
            qbox.add_quad_face(a.into(), b.into(), c.into(), d.into())?;
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

impl<VecT> PolyMeshT<VecT, 3>
where
    VecT: TVec<3> + Add<Output = VecT> + Div<VecT::Scalar, Output = VecT>,
    VecT::Scalar: FromFloat + Add<Output = VecT::Scalar>,
{
    pub fn dual_mesh(&self) -> Result<Self, Error> {
        let mut outmesh =
            Self::with_capacity(self.num_faces(), self.num_edges(), self.num_vertices());
        let points = self.points();
        let points = points.try_borrow()?;
        let fstatus = self.topol.fstatus.try_borrow()?;
        for f in self.faces() {
            if fstatus[f.index() as usize].deleted() {
                return Err(Error::GarbageCollectionRequired);
            }
            outmesh.add_vertex(self.calc_face_centroid(f, &points))?;
        }
        let mut fverts: Vec<VH> = Vec::new();
        for v in self.vertices().filter(|v| !self.is_boundary_vertex(*v)) {
            fverts.clear();
            fverts.extend(self.vf_ccw_iter(v).map(|f| -> VH { f.index().into() }));
            outmesh.add_face(&fverts)?;
        }
        Ok(outmesh)
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
        let _verts = {
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
        mesh.add_tri_face(0.into(), 1.into(), 2.into())?;
        mesh.add_tri_face(0.into(), 2.into(), 3.into())?;
        mesh.add_tri_face(0.into(), 3.into(), 1.into())?;
        mesh.add_tri_face(3.into(), 2.into(), 1.into())?;
        Ok(mesh)
    }

    pub fn hexahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let a = radius * VecT::Scalar::from_f64(1.0f64 / 3.0f64.sqrt());
        let mut mesh = Self::with_capacity(8, 12, 6);
        let _verts = {
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
        mesh.add_quad_face(3.into(), 2.into(), 1.into(), 0.into())?;
        mesh.add_quad_face(2.into(), 6.into(), 5.into(), 1.into())?;
        mesh.add_quad_face(5.into(), 6.into(), 7.into(), 4.into())?;
        mesh.add_quad_face(0.into(), 4.into(), 7.into(), 3.into())?;
        mesh.add_quad_face(3.into(), 7.into(), 6.into(), 2.into())?;
        mesh.add_quad_face(1.into(), 5.into(), 4.into(), 0.into())?;
        Ok(mesh)
    }

    pub fn octahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(6, 12, 8);
        let _verts = {
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
        mesh.add_tri_face(0.into(), 4.into(), 3.into())?;
        mesh.add_tri_face(1.into(), 4.into(), 0.into())?;
        mesh.add_tri_face(2.into(), 4.into(), 1.into())?;
        mesh.add_tri_face(3.into(), 4.into(), 2.into())?;
        mesh.add_tri_face(3.into(), 5.into(), 0.into())?;
        mesh.add_tri_face(0.into(), 5.into(), 1.into())?;
        mesh.add_tri_face(1.into(), 5.into(), 2.into())?;
        mesh.add_tri_face(2.into(), 5.into(), 3.into())?;
        Ok(mesh)
    }

    pub fn icosahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(12, 30, 20);
        let _verts = {
            let mut verts: [VH; 12] = [0.into(); 12];
            mesh.add_vertices(
                &[
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(0.5257311121191336),
                        radius * VecT::Scalar::from_f64(-0.8506508083520399),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.5257311121191336),
                        radius * VecT::Scalar::from_f64(0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.5257311121191336),
                        radius * VecT::Scalar::from_f64(0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(0.5257311121191336),
                        radius * VecT::Scalar::from_f64(0.8506508083520399),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(-0.5257311121191336),
                        radius * VecT::Scalar::from_f64(0.8506508083520399),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(0.5257311121191336),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(-0.5257311121191336),
                        radius * VecT::Scalar::from_f64(-0.8506508083520399),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(-0.5257311121191336),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(0.5257311121191336),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                        radius * VecT::Scalar::from_f64(-0.5257311121191336),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.5257311121191336),
                        radius * VecT::Scalar::from_f64(-0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.5257311121191336),
                        radius * VecT::Scalar::from_f64(-0.8506508083520399),
                        radius * VecT::Scalar::from_f64(0.0),
                    ]),
                ],
                &mut verts,
            )?;
            verts
        };
        mesh.add_tri_face(2.into(), 1.into(), 0.into())?;
        mesh.add_tri_face(1.into(), 2.into(), 3.into())?;
        mesh.add_tri_face(5.into(), 4.into(), 3.into())?;
        mesh.add_tri_face(4.into(), 8.into(), 3.into())?;
        mesh.add_tri_face(7.into(), 6.into(), 0.into())?;
        mesh.add_tri_face(6.into(), 9.into(), 0.into())?;
        mesh.add_tri_face(11.into(), 10.into(), 4.into())?;
        mesh.add_tri_face(10.into(), 11.into(), 6.into())?;
        mesh.add_tri_face(9.into(), 5.into(), 2.into())?;
        mesh.add_tri_face(5.into(), 9.into(), 11.into())?;
        mesh.add_tri_face(8.into(), 7.into(), 1.into())?;
        mesh.add_tri_face(7.into(), 8.into(), 10.into())?;
        mesh.add_tri_face(2.into(), 5.into(), 3.into())?;
        mesh.add_tri_face(8.into(), 1.into(), 3.into())?;
        mesh.add_tri_face(9.into(), 2.into(), 0.into())?;
        mesh.add_tri_face(1.into(), 7.into(), 0.into())?;
        mesh.add_tri_face(11.into(), 9.into(), 6.into())?;
        mesh.add_tri_face(7.into(), 10.into(), 6.into())?;
        mesh.add_tri_face(5.into(), 11.into(), 4.into())?;
        mesh.add_tri_face(10.into(), 8.into(), 4.into())?;
        Ok(mesh)
    }

    pub fn dodecahedron(radius: VecT::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(20, 30, 12);
        let _verts = {
            let mut verts: [VH; 20] = [0.into(); 20];
            mesh.add_vertices(
                &[
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(0.9341723589627157),
                        radius * VecT::Scalar::from_f64(-0.35682208977308993),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(0.9341723589627157),
                        radius * VecT::Scalar::from_f64(0.35682208977308993),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(0.9341723589627157),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(0.9341723589627157),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(-0.9341723589627157),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(-0.9341723589627157),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(-0.9341723589627157),
                        radius * VecT::Scalar::from_f64(0.35682208977308993),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.),
                        radius * VecT::Scalar::from_f64(-0.9341723589627157),
                        radius * VecT::Scalar::from_f64(-0.35682208977308993),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.9341723589627157),
                        radius * VecT::Scalar::from_f64(0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.9341723589627157),
                        radius * VecT::Scalar::from_f64(-0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.9341723589627157),
                        radius * VecT::Scalar::from_f64(0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.9341723589627157),
                        radius * VecT::Scalar::from_f64(-0.35682208977308993),
                        radius * VecT::Scalar::from_f64(0.),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                    ]),
                    VecT::new([
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                        radius * VecT::Scalar::from_f64(-0.5773502691896257),
                        radius * VecT::Scalar::from_f64(0.5773502691896257),
                    ]),
                ],
                &mut verts,
            )?;
        };
        mesh.add_face(&[15.into(), 4.into(), 5.into(), 14.into(), 0.into()])?;
        mesh.add_face(&[15.into(), 0.into(), 1.into(), 13.into(), 10.into()])?;
        mesh.add_face(&[14.into(), 8.into(), 12.into(), 1.into(), 0.into()])?;
        mesh.add_face(&[13.into(), 1.into(), 12.into(), 2.into(), 3.into()])?;
        mesh.add_face(&[19.into(), 3.into(), 2.into(), 18.into(), 6.into()])?;
        mesh.add_face(&[18.into(), 2.into(), 12.into(), 8.into(), 9.into()])?;
        mesh.add_face(&[17.into(), 7.into(), 16.into(), 5.into(), 4.into()])?;
        mesh.add_face(&[17.into(), 4.into(), 15.into(), 10.into(), 11.into()])?;
        mesh.add_face(&[19.into(), 11.into(), 10.into(), 13.into(), 3.into()])?;
        mesh.add_face(&[16.into(), 9.into(), 8.into(), 14.into(), 5.into()])?;
        mesh.add_face(&[19.into(), 6.into(), 7.into(), 17.into(), 11.into()])?;
        mesh.add_face(&[18.into(), 9.into(), 16.into(), 7.into(), 6.into()])?;
        Ok(mesh)
    }
}

#[cfg(test)]
mod test {
    use core::f32;

    use crate::{macros::assert_f32_eq, mesh::PolyMeshF32};

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
        assert_eq!(1., qbox.try_calc_volume().expect("Cannot compute volume"));
        assert_eq!(6., qbox.try_calc_area().expect("Cannot compute area"));
    }

    #[test]
    fn t_tetrahedron() {
        let tet = PolyMeshF32::tetrahedron(1.0).expect("Cannot create a tetrahedron");
        assert_eq!(4, tet.num_vertices());
        assert_eq!(12, tet.num_halfedges());
        assert_eq!(6, tet.num_edges());
        assert_eq!(4, tet.num_faces());
        assert_eq!(
            8.0 / 3.0f32.sqrt(),
            tet.try_calc_area().expect("Cannot compute area")
        );
        assert_f32_eq!(
            8.0 / (9.0 * 3.0f32.sqrt()),
            tet.try_calc_volume().expect("Cannot compute volume")
        );
    }

    #[test]
    fn t_hexahedron() {
        let hex = PolyMeshF32::hexahedron(1.0).expect("Cannot create hexahedron");
        assert_eq!(hex.num_vertices(), 8);
        assert_eq!(hex.num_halfedges(), 24);
        assert_eq!(hex.num_edges(), 12);
        assert_eq!(hex.num_faces(), 6);
        assert_f32_eq!(8.0, hex.try_calc_area().expect("Cannot compute area"), 1e-6);
        assert_f32_eq!(
            8.0 / (3.0 * 3.0f32.sqrt()),
            hex.try_calc_volume().expect("Cannot compute volume")
        );
    }

    #[test]
    fn t_octahedron() {
        let oct = PolyMeshF32::octahedron(1.0).expect("Cannot create octahedron");
        assert_eq!(oct.num_vertices(), 6);
        assert_eq!(oct.num_halfedges(), 24);
        assert_eq!(oct.num_edges(), 12);
        assert_eq!(oct.num_faces(), 8);
        assert_eq!(
            4.0 * 3.0f32.sqrt(),
            oct.try_calc_area().expect("Cannot compute area")
        );
        assert_f32_eq!(
            4.0 / 3.0,
            oct.try_calc_volume().expect("Cannot compute volume")
        );
    }

    #[test]
    fn t_icosahedron() {
        let ico = PolyMeshF32::icosahedron(1.0).expect("Cannot create icosahedron");
        assert_eq!(12, ico.num_vertices());
        assert_eq!(60, ico.num_halfedges());
        assert_eq!(30, ico.num_edges());
        assert_eq!(20, ico.num_faces());
        assert_f32_eq!(
            {
                let phi = (1.0 + 5.0f32.sqrt()) / 2.0;
                20.0 * 3.0f32.sqrt() / (phi * phi + 1.0)
            },
            ico.try_calc_area().expect("Cannot compute area"),
            1e-6
        );
        assert_f32_eq!(
            {
                let phi = (1.0 + 5.0f32.sqrt()) / 2.0;
                20.0 * phi * phi / (3.0 * (phi * phi + 1.0) * (phi * phi + 1.0).sqrt())
            },
            ico.try_calc_volume().expect("Cannot compute volume"),
            1e-6
        );
    }

    #[test]
    fn t_dodecahedron() {
        let dod = PolyMeshF32::dodecahedron(1.0).expect("Cannot create dodecahedron");
        assert_eq!(20, dod.num_vertices());
        assert_eq!(60, dod.num_halfedges());
        assert_eq!(30, dod.num_edges());
        assert_eq!(12, dod.num_faces());
        assert_f32_eq!(
            {
                let phi = (1.0 + 5.0f32.sqrt()) / 2.0;
                20.0f32 / (phi * (3.0f32 - phi).sqrt())
            },
            dod.try_calc_area().expect("Cannot compute area"),
            1e-6
        );
        assert_f32_eq!(
            {
                let phi = (1.0 + 5.0f32.sqrt()) / 2.0;
                40.0 / (3.0 * 3.0f32.sqrt() * (6.0 - 2.0 * phi))
            },
            dod.try_calc_volume().expect("Cannot compute volume")
        );
    }
}
