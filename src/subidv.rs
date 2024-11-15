use crate::{iterator, topol::Topology, Adaptor, Error, FloatScalarAdaptor, Handle, PolyMeshT, HH};
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
    A: Adaptor<DIM> + FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    fn count_topol(mesh: &Topology, niter: usize) -> (usize, usize, usize) {
        debug_assert!(niter > 0);
        let mut nv = mesh.num_vertices();
        let mut ne = mesh.num_edges();
        let mut nf = mesh.num_faces();
        for i in 1..niter {
            let v = nv + ne + nf;
            let f = if i == 1 {
                mesh.faces().map(|f| mesh.face_valence(f)).sum::<usize>()
            } else {
                nf * 4
            };
            let e = 2 * ne + f;
            (nv, ne, nf) = (v, e, f);
        }
        (nv, ne, nf)
    }

    fn calc_face_edge_points(
        mesh: &mut PolyMeshT<DIM, A>,
        update_points: bool,
        fpos: &mut Vec<A::Vector>,
        epos: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        let points = mesh.points();
        let points = points.try_borrow()?;
        // Compute face points.
        fpos.clear();
        fpos.extend(mesh.faces().map(|f| mesh.calc_face_centroid(f, &points)));
        // Compute edge points.
        epos.clear();
        epos.extend(mesh.edges().map(|e| {
            let (h, oh) = mesh.halfedge_pair(e);
            match (mesh.halfedge_face(h), mesh.halfedge_face(oh)) {
                (Some(fa), Some(fb)) if update_points => {
                    (points[mesh.head_vertex(h).index() as usize]
                        + points[mesh.head_vertex(oh).index() as usize]
                        + fpos[fa.index() as usize]
                        + fpos[fb.index() as usize])
                        * A::scalarf64(0.25)
                }
                _ => {
                    (points[mesh.head_vertex(h).index() as usize]
                        + points[mesh.head_vertex(oh).index() as usize])
                        * A::scalarf64(0.5)
                }
            }
        }));
        Ok(())
    }

    fn update_vertex_positions(
        mesh: &mut PolyMeshT<DIM, A>,
        fpos: &[A::Vector],
        epos: &[A::Vector],
        vpos: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        // Compute vertex positions.
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut()?;
        vpos.clear();
        {
            let points: &[A::Vector] = &points;
            vpos.extend(mesh.vertices().map(|v| {
                if mesh.is_boundary_vertex(v) {
                    let (count, sum) = mesh
                        .ve_ccw_iter(v)
                        .filter(|e| mesh.is_boundary_edge(*e))
                        .fold((1usize, points[v.index() as usize]), |(count, total), e| {
                            (count + 1, total + epos[e.index() as usize])
                        });
                    sum / A::scalarf64(count as f64)
                } else {
                    let valence = mesh.vertex_valence(v) as f64;
                    (((mesh.vf_ccw_iter(v).fold(A::zero_vector(), |total, f| {
                        total + fpos[f.index() as usize]
                    }) + mesh.vv_ccw_iter(v).fold(A::zero_vector(), |total, v| {
                        total + points[v.index() as usize]
                    })) / A::scalarf64(valence))
                        + (points[v.index() as usize] * A::scalarf64(valence - 2.0)))
                        / A::scalarf64(valence)
                }
            }));
        }
        // Update the vertex positions.
        points.copy_from_slice(vpos);
        Ok(())
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM> + FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    pub fn subidivide_catmull_clark(
        &mut self,
        iters: usize,
        update_points: bool,
    ) -> Result<(), Error> {
        if iters == 0 {
            return Ok(());
        }
        {
            // Ensure no deleted topology.
            let fstatus = self.topol.fstatus.try_borrow()?;
            if let Some(f) = self.faces().find(|f| fstatus[f.index() as usize].deleted()) {
                return Err(Error::DeletedFace(f));
            }
            let estatus = self.topol.estatus.try_borrow()?;
            if let Some(e) = self.edges().find(|e| estatus[e.index() as usize].deleted()) {
                return Err(Error::DeletedEdge(e));
            }
            let vstatus = self.topol.hstatus.try_borrow()?;
            if let Some(v) = self
                .vertices()
                .find(|v| vstatus[v.index() as usize].deleted())
            {
                return Err(Error::DeletedVertex(v));
            }
        }
        // Use vectors instead of properties because we don't want these to
        // change when we add new topology.
        let (nverts, nedges, nfaces) = CatmullClark::<DIM, A>::count_topol(&self.topol, iters);
        self.reserve(nverts, nedges, nfaces)?;
        let mut fpos = Vec::with_capacity(nfaces);
        let mut epos = Vec::with_capacity(nedges);
        let mut vpos = if update_points {
            Some(Vec::with_capacity(nverts))
        } else {
            None
        };
        // Temporary storage to use inside the loop.
        let mut fhs = Vec::new();
        let mut hhs = Vec::new();
        let mut ffs = Vec::new();
        for _ in 0..iters {
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
                self.split_edge((ei as u32).into(), *pos, true)?;
            }
            for f in self.faces() {
                // Find a halfedge that points to an old vertex, and collect the
                // loop of halfedges starting from there.
                let hstart = self
                    .fh_ccw_iter(f)
                    .find(|&h| self.head_vertex(h).index() < num_old_verts)
                    .ok_or(Error::CannotSplitFace(f))?;
                fhs.clear();
                fhs.extend(iterator::loop_ccw_iter(&self.topol, hstart));
                let fhs: &[HH] = &fhs; // Immutable.
                debug_assert!(fhs.len() % 2 == 0);
                let valence = iterator::loop_ccw_iter(&self.topol, hstart).count() / 2;
                debug_assert_eq!(valence * 2, fhs.len());
                let ne = self.num_edges();
                // New vertex in the middle.
                let fv = self.add_vertex(fpos[f.index() as usize])?;
                // Create new edges and faces, with some math to get the indices
                // of edges we haven't added yet. This is a bit sketchy, but
                // should be safe with ample testing.
                hhs.clear();
                ffs.clear();
                for (lei, hpair) in fhs.chunks_exact(2).enumerate() {
                    let h1 = hpair[0];
                    let pei = (ne + ((lei + valence - 1) % valence)) as u32;
                    let nei = (ne + ((lei + 1) % valence)) as u32;
                    let enew = self.topol.new_edge(
                        self.tail_vertex(h1),
                        fv,
                        self.prev_halfedge(h1),
                        (2 * pei + 1).into(),
                        (2 * nei).into(),
                        h1,
                    )?;
                    debug_assert_eq!(enew.index(), (lei + ne) as u32);
                    hhs.push(self.edge_halfedge(enew, false));
                    let flocal = if h1 == hstart {
                        f
                    } else {
                        self.topol.new_face(h1)?
                    };
                    ffs.push(flocal);
                }
                for (i, hpair) in fhs.chunks_exact(2).enumerate() {
                    let rh = hhs[i];
                    let orh = self.opposite_halfedge(rh);
                    let flocal = ffs[i];
                    let pflocal = ffs[(i + valence - 1) % valence];
                    let h1 = hpair[0];
                    let h2 = hpair[1];
                    // Link halfedges and faces.
                    self.topol.face_mut(flocal).halfedge = h1;
                    self.topol.halfedge_mut(rh).face = Some(pflocal);
                    self.topol.halfedge_mut(orh).face = Some(flocal);
                    self.topol.halfedge_mut(h1).face = Some(flocal);
                    self.topol.halfedge_mut(h2).face = Some(flocal);
                    // Link halfedges.
                    self.topol.link_halfedges(self.prev_halfedge(rh), rh);
                    self.topol.link_halfedges(orh, h1);
                    self.topol.link_halfedges(h1, h2);
                }
                self.topol.vertex_mut(fv).halfedge = Some(self.opposite_halfedge(hhs[0]));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::{alum_glam::PolyMeshF32, obj::test::bunny_mesh};

    #[test]
    fn t_box_catmull_clark() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
        mesh.subidivide_catmull_clark(1, true)
            .expect("Subdivision failed");
        assert_eq!(26, mesh.num_vertices());
        assert_eq!(48, mesh.num_edges());
        assert_eq!(24, mesh.num_faces());
        mesh.check_topology().expect("Topological errors found");
    }

    #[test]
    fn t_bunny_subdiv() {
        let mut mesh = bunny_mesh();
        mesh.subidivide_catmull_clark(3, true)
            .expect("Cannot subivide");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(mesh.try_calc_area().expect("Cannot compute area"), 5.566642);
        assert_eq!(238464, mesh.num_faces());
        assert_eq!(477096, mesh.num_edges());
        assert_eq!(238630, mesh.num_vertices());
    }
}
