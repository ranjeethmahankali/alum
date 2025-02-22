use crate::{
    EH, EditableTopology, Error, FH, FloatScalarAdaptor, HH, Handle, HasIterators, HasTopology,
    PolyMeshT, topol::Topology,
};
use std::{
    marker::PhantomData,
    ops::{Add, Div, Mul},
};

/// This struct doesn't conatin any data. This purpose of this struct is to
/// provide a convenient scope inside its `impl` where I can impose all the
/// trait bounds once and write several free functions within that scope, with
/// those trait bounds applied.
struct CatmullClark<const DIM: usize, A>(PhantomData<A>);

impl<const DIM: usize, A> CatmullClark<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    /// Reserve the space for the final mesh after `niter` iterations of
    /// subdivision. Also reserve space for buffers used to store the points
    /// corresponding to vertices, edges and faces. These buffers are used to
    /// store the points corresponding to the mesh before subdivision, so the
    /// capacity of these buffers will correspond to one less iteration of
    /// subdivision (hence less memory required) than that of the mesh
    /// itself. This reserve takes care of that.
    fn reserve(
        niter: usize,
        mesh: &mut Topology,
        mut vertex_points: Option<&mut Vec<A::Vector>>,
        edge_points: &mut Vec<A::Vector>,
        face_points: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        if let Some(vpts) = vertex_points.as_mut() {
            vpts.clear();
        }
        edge_points.clear();
        face_points.clear();
        debug_assert!(niter > 0);
        let mut nv = mesh.num_vertices();
        let mut ne = mesh.num_edges();
        let mut nf = mesh.num_faces();
        for i in 0..niter {
            if i == niter - 1 {
                // Reserve the pos arrays.
                if let Some(vpts) = vertex_points.as_mut() {
                    vpts.reserve(nv);
                }
                edge_points.reserve(ne);
                face_points.reserve(nf);
            }
            let v = nv + ne + nf;
            let f = if i == 0 {
                mesh.faces().map(|f| f.valence(mesh)).sum::<usize>()
            } else {
                nf * 4
            };
            let e = 2 * ne + f;
            (nv, ne, nf) = (v, e, f);
        }
        mesh.reserve(nv, ne, nf)
    }

    /// Compute the locations of the points to split the faces and edges.
    fn calc_face_edge_points(
        mesh: &PolyMeshT<DIM, A>,
        update_points: bool,
        face_points: &mut Vec<A::Vector>,
        edge_points: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        let points = mesh.points();
        let points = points.try_borrow()?;
        // Compute face points.
        face_points.clear();
        face_points.extend(mesh.faces().map(|f| mesh.calc_face_centroid(f, &points)));
        // Compute edge points.
        edge_points.clear();
        edge_points.extend(mesh.edges().map(|e| {
            let (h, oh) = e.halfedges();
            match (h.face(mesh), oh.face(mesh)) {
                (Some(fa), Some(fb)) if update_points => {
                    (points[h.head(mesh)]
                        + points[oh.head(mesh)]
                        + face_points[fa.index() as usize]
                        + face_points[fb.index() as usize])
                        * A::scalarf64(0.25)
                }
                _ => (points[h.head(mesh)] + points[oh.head(mesh)]) * A::scalarf64(0.5),
            }
        }));
        Ok(())
    }

    /// Compute the new locations of the vertices, and move them there. This
    /// function makes use of the preallocated buffer `vertex_points` to avoid
    /// repeated allocations.
    fn update_vertex_positions(
        mesh: &mut PolyMeshT<DIM, A>,
        face_points: &[A::Vector],
        edge_points: &[A::Vector],
        vertex_points: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        // Compute vertex positions.
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut()?;
        vertex_points.clear();
        {
            vertex_points.extend(mesh.vertices().map(|v| {
                if v.is_boundary(mesh) {
                    let (count, sum) = mesh
                        .ve_ccw_iter(v)
                        .filter(|e| e.is_boundary(mesh))
                        .fold((1usize, points[v]), |(count, total), e| {
                            (count + 1, total + edge_points[e.index() as usize])
                        });
                    sum / A::scalarf64(count as f64)
                } else {
                    let valence = v.valence(mesh) as f64;
                    (((mesh.vf_ccw_iter(v).fold(A::zero_vector(), |total, f| {
                        total + face_points[f.index() as usize]
                    }) + mesh
                        .vv_ccw_iter(v)
                        .fold(A::zero_vector(), |total, v| total + points[v]))
                        / A::scalarf64(valence))
                        + (points[v] * A::scalarf64(valence - 2.0)))
                        / A::scalarf64(valence)
                }
            }));
        }
        // Update the vertex positions.
        points.copy_from_slice(vertex_points);
        Ok(())
    }

    /// Split a face according to the catmull clark scheme. This must be called
    /// after all the edges of the face are already split. `hloop`, `spliths`,
    /// and `subfaces` are preallocated buffers used by this function to avoid
    /// repeated allocations inside the call.
    fn split_face(
        mesh: &mut PolyMeshT<DIM, A>,
        f: FH,
        num_old_verts: u32,
        face_points: &[A::Vector],
        hloop: &mut Vec<HH>,
        spliths: &mut Vec<HH>,
        subfaces: &mut Vec<FH>,
    ) -> Result<(), Error> {
        // Find a halfedge that points to an old vertex, and collect the
        // loop of halfedges starting from there.
        let hstart = mesh
            .fh_ccw_iter(f)
            .find(|&h| h.head(mesh).index() < num_old_verts)
            .ok_or(Error::CannotSplitFace(f))?;
        hloop.clear();
        hloop.extend(mesh.loop_ccw_iter(hstart));
        let hloop: &[HH] = hloop; // Immutable.
        debug_assert!(hloop.len() % 2 == 0);
        let valence = mesh.loop_ccw_iter(hstart).count() / 2;
        debug_assert_eq!(valence * 2, hloop.len());
        let ne = mesh.num_edges();
        // New vertex in the middle.
        let fv = mesh.add_vertex(face_points[f.index() as usize])?;
        // Create new edges and faces, with some math to get the indices
        // of edges we haven't added yet. This is a bit sketchy, but
        // should be safe with ample testing.
        spliths.clear();
        subfaces.clear();
        for (lei, hpair) in hloop.chunks_exact(2).enumerate() {
            let h1 = hpair[0];
            let pei = (ne + ((lei + valence - 1) % valence)) as u32;
            let nei = (ne + ((lei + 1) % valence)) as u32;
            let enew = mesh.topol.new_edge(
                h1.tail(mesh),
                fv,
                h1.prev(mesh),
                (2 * pei + 1).into(),
                (2 * nei).into(),
                h1,
            )?;
            debug_assert_eq!(enew.index(), (lei + ne) as u32);
            spliths.push(enew.halfedge(false));
            let flocal = if h1 == hstart {
                f
            } else {
                mesh.topol.new_face(h1)?
            };
            subfaces.push(flocal);
        }
        for (i, hpair) in hloop.chunks_exact(2).enumerate() {
            let rh = spliths[i];
            let orh = rh.opposite();
            let flocal = subfaces[i];
            let pflocal = subfaces[(i + valence - 1) % valence];
            let h1 = hpair[0];
            let h2 = hpair[1];
            // Link halfedges and faces.
            mesh.topol.face_mut(flocal).halfedge = h1;
            mesh.topol.halfedge_mut(rh).face = Some(pflocal);
            mesh.topol.halfedge_mut(orh).face = Some(flocal);
            mesh.topol.halfedge_mut(h1).face = Some(flocal);
            mesh.topol.halfedge_mut(h2).face = Some(flocal);
            // Link halfedges.
            mesh.topol.link_halfedges(rh.prev(mesh), rh);
            mesh.topol.link_halfedges(orh, h1);
            mesh.topol.link_halfedges(h1, h2);
        }
        mesh.topol.vertex_mut(fv).halfedge = Some(spliths[0].opposite());
        Ok(())
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    /// Subdivide the mesh according to the [Catmull-Clark
    /// scheme](https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface).
    ///
    /// Subdivisions are carried out for the given number of
    /// `iterations`. `update_points` determines whether the vertices of the
    /// mesh are moved. If this is `false`, only the topology of the mesh is
    /// subdivided, leaving the shape of the mesh unmodified. If this is `true`,
    /// then the shape is smoothed according to the Catmull-Clark scheme.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
    /// assert_eq!((8, 12, 6), (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces()));
    /// mesh.subdivide_catmull_clark(1, true).expect("Subdivision failed");
    /// // The mesh now has more faces.
    /// assert_eq!((26, 48, 24), (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces()));
    /// mesh.check_topology().expect("Topological errors found");
    /// ```
    pub fn subdivide_catmull_clark(
        &mut self,
        iterations: usize,
        update_points: bool,
    ) -> Result<(), Error> {
        if iterations == 0 {
            return Ok(());
        }
        self.check_for_deleted()?;
        let mut fpos = Vec::new();
        let mut epos = Vec::new();
        let mut vpos = if update_points {
            Some(Vec::new())
        } else {
            None
        };
        // Use vectors instead of properties because we don't want these to
        // change when we add new topology.
        CatmullClark::<DIM, A>::reserve(
            iterations,
            &mut self.topol,
            vpos.as_mut(),
            &mut epos,
            &mut fpos,
        )?;
        // Temporary storage to use inside the loop.
        let mut hloop = Vec::new();
        let mut spliths = Vec::new();
        let mut subfaces = Vec::new();
        for _ in 0..iterations {
            CatmullClark::<DIM, A>::calc_face_edge_points(
                self,
                update_points,
                &mut fpos,
                &mut epos,
            )?;
            // Make them immutable from here.
            let fpos: &[A::Vector] = &fpos;
            let epos: &[A::Vector] = &epos;
            if let Some(vpos) = &mut vpos {
                CatmullClark::<DIM, A>::update_vertex_positions(self, fpos, epos, vpos)?;
            }
            let num_old_verts = self.num_vertices() as u32;
            for (ei, pos) in epos.iter().enumerate() {
                let ev = self.add_vertex(*pos)?;
                let e: EH = (ei as u32).into();
                self.split_edge(e.halfedge(false), ev, true)?;
            }
            for f in self.faces() {
                CatmullClark::<DIM, A>::split_face(
                    self,
                    f,
                    num_old_verts,
                    fpos,
                    &mut hloop,
                    &mut spliths,
                    &mut subfaces,
                )?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::{HasTopology, obj::test::bunny_mesh, use_glam::PolyMeshF32};

    #[test]
    fn t_box_catmull_clark() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
        mesh.subdivide_catmull_clark(1, true)
            .expect("Subdivision failed");
        assert_eq!(26, mesh.num_vertices());
        assert_eq!(48, mesh.num_edges());
        assert_eq!(24, mesh.num_faces());
        mesh.check_topology().expect("Topological errors found");
    }

    #[test]
    fn t_bunny_subdiv_catmull_clark() {
        let mut mesh = bunny_mesh();
        mesh.subdivide_catmull_clark(3, true)
            .expect("Cannot subivide");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(mesh.try_calc_area().expect("Cannot compute area"), 5.566642);
        assert_eq!(238464, mesh.num_faces());
        assert_eq!(477096, mesh.num_edges());
        assert_eq!(238630, mesh.num_vertices());
    }
}
