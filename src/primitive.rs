use std::ops::{Add, Div, Mul, Neg};

use crate::{
    element::{Handle, VH},
    error::Error,
    mesh::{Adaptor, FloatScalarAdaptor, PolyMeshT},
};

impl<A> PolyMeshT<3, A>
where
    A: Adaptor<3>,
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
    pub fn quad_box(min: A::Vector, max: A::Vector) -> Result<Self, Error> {
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
        let mut pos: [A::Vector; 8] = [A::zero_vector(); 8];
        for (i, &(xf, yf, zf)) in BOX_POS.iter().enumerate() {
            pos[i] = A::vector([
                A::vector_coord(if xf { &max } else { &min }, 0),
                A::vector_coord(if yf { &max } else { &min }, 1),
                A::vector_coord(if zf { &max } else { &min }, 2),
            ]);
        }
        let verts = qbox.add_vertices(&pos)?;
        assert_eq!(
            verts,
            0..8,
            "Vertices are expected to be in one contiguous range"
        );
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
        A: FloatScalarAdaptor<3>,
    {
        Self::quad_box(
            A::vector([A::scalarf64(0.0); 3]),
            A::vector([A::scalarf64(1.0); 3]),
        )
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM> + FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    /// Create a mesh that is the dual of this mesh.
    ///
    /// <https://en.wikipedia.org/wiki/Dual_polyhedron>.
    ///
    /// The dual mesh will contain vertices at the centroids of the input mesh,
    /// and an edge connecting the vertices for every pair of adjacent faces in
    /// the input mesh. The output mesh will contain as many vertices as the
    /// faces of the input mesh, as many edges as the interior edges of the
    /// input mesh, and as many faces as the interior vertices of the input
    /// mesh.
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
impl<A> PolyMeshT<3, A>
where
    A: Adaptor<3> + FloatScalarAdaptor<3>,
    A::Scalar: Mul<Output = A::Scalar> + Neg<Output = A::Scalar>,
{
    /// Create a tetrahedron centered at the origin, with a unit
    /// circumradius. The vertices of the mesh will lie on the unit sphere.
    pub fn tetrahedron(radius: A::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(4, 6, 4);
        let a = radius * A::scalarf64(1.0f64 / 3.0);
        let b = radius * A::scalarf64((8.0 / 9.0f64).sqrt());
        let c = radius * A::scalarf64((2.0 / 9.0f64).sqrt());
        let d = radius * A::scalarf64((2.0 / 3.0f64).sqrt());
        let verts = mesh.add_vertices(&[
            A::vector([A::scalarf64(0.0), A::scalarf64(0.0), A::scalarf64(1.0)]),
            A::vector([-c, d, -a]),
            A::vector([-c, -d, -a]),
            A::vector([b, A::scalarf64(0.0), -a]),
        ])?;
        assert_eq!(
            verts,
            0..4,
            "Vertices are expected to be in one contiguous range"
        );
        mesh.add_tri_face(0.into(), 1.into(), 2.into())?;
        mesh.add_tri_face(0.into(), 2.into(), 3.into())?;
        mesh.add_tri_face(0.into(), 3.into(), 1.into())?;
        mesh.add_tri_face(3.into(), 2.into(), 1.into())?;
        Ok(mesh)
    }

    /// Create a hexehedron centered at the origin, with a unit
    /// circumradius. The vertices of the mesh will lie on the unit sphere.
    pub fn hexahedron(radius: A::Scalar) -> Result<Self, Error> {
        let a = radius * A::scalarf64(1.0f64 / 3.0f64.sqrt());
        let mut mesh = Self::with_capacity(8, 12, 6);
        let verts = mesh.add_vertices(&[
            A::vector([-a, -a, -a]),
            A::vector([a, -a, -a]),
            A::vector([a, a, -a]),
            A::vector([-a, a, -a]),
            A::vector([-a, -a, a]),
            A::vector([a, -a, a]),
            A::vector([a, a, a]),
            A::vector([-a, a, a]),
        ])?;
        assert_eq!(
            verts,
            0..8,
            "Vertices are expected to be in one contiguous range"
        );
        mesh.add_quad_face(3.into(), 2.into(), 1.into(), 0.into())?;
        mesh.add_quad_face(2.into(), 6.into(), 5.into(), 1.into())?;
        mesh.add_quad_face(5.into(), 6.into(), 7.into(), 4.into())?;
        mesh.add_quad_face(0.into(), 4.into(), 7.into(), 3.into())?;
        mesh.add_quad_face(3.into(), 7.into(), 6.into(), 2.into())?;
        mesh.add_quad_face(1.into(), 5.into(), 4.into(), 0.into())?;
        Ok(mesh)
    }

    /// Create an octahedron centered at the origin, with a unit
    /// circumradius. The vertices of the mesh will lie on the unit sphere.
    pub fn octahedron(radius: A::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(6, 12, 8);
        let zero = A::scalarf64(0.0);
        let verts = mesh.add_vertices(&[
            A::vector([radius, zero, zero]),
            A::vector([zero, radius, zero]),
            A::vector([-radius, zero, zero]),
            A::vector([zero, -radius, zero]),
            A::vector([zero, zero, radius]),
            A::vector([zero, zero, -radius]),
        ])?;
        assert_eq!(
            verts,
            0..6,
            "Vertices are expected to be in one contiguous range"
        );
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

    /// Create an icosahedron centered at the origin, with a unit
    /// circumradius. The vertices of the mesh will lie on the unit sphere.
    pub fn icosahedron(radius: A::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(12, 30, 20);
        let verts = mesh.add_vertices(&[
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.5257311121191336),
                radius * A::scalarf64(-0.8506508083520399),
            ]),
            A::vector([
                radius * A::scalarf64(0.5257311121191336),
                radius * A::scalarf64(0.8506508083520399),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(-0.5257311121191336),
                radius * A::scalarf64(0.8506508083520399),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.5257311121191336),
                radius * A::scalarf64(0.8506508083520399),
            ]),
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.5257311121191336),
                radius * A::scalarf64(0.8506508083520399),
            ]),
            A::vector([
                radius * A::scalarf64(-0.8506508083520399),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.5257311121191336),
            ]),
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.5257311121191336),
                radius * A::scalarf64(-0.8506508083520399),
            ]),
            A::vector([
                radius * A::scalarf64(0.8506508083520399),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.5257311121191336),
            ]),
            A::vector([
                radius * A::scalarf64(0.8506508083520399),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.5257311121191336),
            ]),
            A::vector([
                radius * A::scalarf64(-0.8506508083520399),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.5257311121191336),
            ]),
            A::vector([
                radius * A::scalarf64(0.5257311121191336),
                radius * A::scalarf64(-0.8506508083520399),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(-0.5257311121191336),
                radius * A::scalarf64(-0.8506508083520399),
                radius * A::scalarf64(0.0),
            ]),
        ])?;
        assert_eq!(
            verts,
            0..12,
            "Vertices are expected to be in one contiguous range"
        );
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

    /// Create a dodecahedron centered at the origin, with a unit
    /// circumradius. The vertices of the mesh will lie on the unit sphere.
    pub fn dodecahedron(radius: A::Scalar) -> Result<Self, Error> {
        let mut mesh = Self::with_capacity(20, 30, 12);
        let verts = mesh.add_vertices(&[
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.9341723589627157),
                radius * A::scalarf64(-0.35682208977308993),
            ]),
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.9341723589627157),
                radius * A::scalarf64(0.35682208977308993),
            ]),
            A::vector([
                radius * A::scalarf64(-0.35682208977308993),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.9341723589627157),
            ]),
            A::vector([
                radius * A::scalarf64(0.35682208977308993),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(0.9341723589627157),
            ]),
            A::vector([
                radius * A::scalarf64(0.35682208977308993),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.9341723589627157),
            ]),
            A::vector([
                radius * A::scalarf64(-0.35682208977308993),
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.9341723589627157),
            ]),
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.9341723589627157),
                radius * A::scalarf64(0.35682208977308993),
            ]),
            A::vector([
                radius * A::scalarf64(0.0),
                radius * A::scalarf64(-0.9341723589627157),
                radius * A::scalarf64(-0.35682208977308993),
            ]),
            A::vector([
                radius * A::scalarf64(-0.9341723589627157),
                radius * A::scalarf64(0.35682208977308993),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(-0.9341723589627157),
                radius * A::scalarf64(-0.35682208977308993),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(0.9341723589627157),
                radius * A::scalarf64(0.35682208977308993),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(0.9341723589627157),
                radius * A::scalarf64(-0.35682208977308993),
                radius * A::scalarf64(0.0),
            ]),
            A::vector([
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
            ]),
            A::vector([
                radius * A::scalarf64(0.5773502691896257),
                radius * A::scalarf64(-0.5773502691896257),
                radius * A::scalarf64(0.5773502691896257),
            ]),
        ])?;
        assert_eq!(
            verts,
            0..20,
            "Vertices are expected to be in one contiguous range"
        );
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

#[cfg(all(test, feature = "use_glam"))]
mod test {
    use core::f32;

    use crate::{alum_glam::PolyMeshF32, macros::assert_f32_eq};

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
