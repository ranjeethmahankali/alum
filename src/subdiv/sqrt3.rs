use core::f64;
use std::{
    marker::PhantomData,
    ops::{Add, Div, Mul},
};

use crate::{
    topol::Topology, EditableTopology, Error, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, PolyMeshT, VPropRef, EH,
};

const WEIGHTS: [(f64, f64); 64] = [
    (f64::MAX, f64::MAX), // Should never happen.
    (0.7777777777777778, 0.2222222222222222),
    (0.33333333333333337, 0.3333333333333333),
    (0.4444444444444444, 0.1851851851851852),
    (0.5555555555555556, 0.1111111111111111),
    (0.6242259987499883, 0.07515480025000235),
    (0.6666666666666667, 0.05555555555555555),
    (0.6941088448574964, 0.04369873644892909),
    (0.7126903958192328, 0.035913700522595904),
    (0.7257876540264395, 0.03046803844150672),
    (0.7353371098610995, 0.026466289013890053),
    (0.7425007850735958, 0.023409019538764012),
    (0.7480056452854309, 0.02099952955954743),
    (0.7523235612562689, 0.01905203374951778),
    (0.7557708595338709, 0.017444938604723507),
    (0.7585656572539112, 0.016095622849739248),
    (0.7608621183358415, 0.014946117604009906),
    (0.7627716065343013, 0.01395461138033522),
    (0.7643761379524241, 0.013090214558198662),
    (0.7657371648223632, 0.012329622904086146),
    (0.7669014480655897, 0.011654927596720516),
    (0.7679050679524757, 0.011052139621310682),
    (0.7687762163587772, 0.010510171983691944),
    (0.7695371749661777, 0.010020122827557494),
    (0.7702057391753485, 0.009574760867693813),
    (0.7707962580285846, 0.009168149678856614),
    (0.771320403872456, 0.008795369081828614),
    (0.7717877490177386, 0.008452305591935607),
    (0.7722062027070719, 0.008135492760461719),
    (0.7725823457133526, 0.00784198807884991),
    (0.7729216890519568, 0.0075692770316014395),
    (0.7732288758338877, 0.00731519755374556),
    (0.7735078400896067, 0.0070778799971997885),
    (0.7737619327250459, 0.006855699008331943),
    (0.7739940221519782, 0.006647234642588876),
    (0.7742065752441399, 0.006451240707310287),
    (0.7744017228916018, 0.006266618808566616),
    (0.7745813134106019, 0.006092396934848598),
    (0.7747469563117161, 0.005927711676007472),
    (0.7749000583639806, 0.005771793375282549),
    (0.7750418534655862, 0.005623953663360345),
    (0.7751734275067441, 0.005483574938859899),
    (0.7752957391611397, 0.005350101448544294),
    (0.7754096373500245, 0.005223031689534314),
    (0.7755158759735407, 0.005101911909692259),
    (0.7756151263870157, 0.004986330524732986),
    (0.7757079880080735, 0.004875913304172316),
    (0.775794997367729, 0.004770319204941938),
    (0.7758766358608468, 0.0046692367528990255),
    (0.7759533364051658, 0.004572380889690494),
    (0.7760254891809951, 0.004479490216380098),
    (0.7760934465937636, 0.004390324576592872),
    (0.7761575275773454, 0.004304662931204898),
    (0.7762180213363371, 0.004222301484220055),
    (0.7762751906093207, 0.004143052025753321),
    (0.7763292745219, 0.004066740463238181),
    (0.7763804910873873, 0.003993205516296656),
    (0.7764290394039888, 0.003922297554315986),
    (0.7764751015898577, 0.003853877558795557),
    (0.7765188444911416, 0.003787816195065397),
    (0.7765604211929497, 0.003723992980117506),
    (0.7765999723587975, 0.003662295535101679),
    (0.7766376274204211, 0.003602618912573853),
    (0.7766735056367559, 0.0035448649898927637),
];

fn compute_weight(valence: usize) -> (f64, f64) {
    let valence = valence as f64;
    let alpha = (4.0 - 2.0 * f64::cos(f64::consts::TAU / valence)) / 9.0;
    (1.0 - alpha, alpha / valence)
}

fn get_weight(valence: usize) -> (f64, f64) {
    match WEIGHTS.get(valence) {
        Some((a, b)) => (*a, *b),
        None => compute_weight(valence),
    }
}

struct Sqrt3Scheme<const DIM: usize, A>(PhantomData<A>);

impl<const DIM: usize, A> Sqrt3Scheme<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>,
{
    fn compute_edge_points(
        mesh: &Topology,
        e: EH,
        points: &VPropRef<A::Vector>,
    ) -> (A::Vector, A::Vector) {
        let h = {
            let mut h = e.halfedge(false);
            if !h.is_boundary(mesh) {
                h = e.halfedge(true);
            }
            h
        };
        let (p1, p2, p3, p4) = (
            points[h.next(mesh).head(mesh)],
            points[h.head(mesh)],
            points[h.tail(mesh)],
            points[h.prev(mesh).tail(mesh)],
        );
        (
            (p1 + p2 * A::scalarf64(16.0) + p3 * A::scalarf64(10.0)) / A::scalarf64(27.0),
            (p2 * A::scalarf64(10.0) + p3 * A::scalarf64(16.0) + p4) / A::scalarf64(27.0),
        )
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>,
{
    /// Subdivide the mesh using the [sqrt-3 subdivision
    /// scheme](https://www.graphics.rwth-aachen.de/publication/03138/).
    ///
    /// This subdivision scheme works with triangles, so if the given mesh is
    /// not a strictly triangle mesh, it will be triangulated before any
    /// subdivisions are performed. Subdivision is performed for the given
    /// number of `iterations`.
    ///
    /// Consecutive iterations of the this subdivision scheme are not
    /// identical. The `phase` flag tells this implementation how to subdivide
    /// the mesh. This is treated as a suggestion, and can be ignored in order
    /// to preserve the topology of the mesh. Essentially, `phase` is meant to
    /// be used to keep track of the internal state of the subdivision when
    /// calling this function multiple times. When calling this function for the
    /// first time, set the `phase` to `false`. Store the returned boolean value
    /// and pass it in as the phase for the subsequent call to this
    /// function. While it is not necessary to keep track of the phase, but
    /// doing so will result in better quality subdivision.
    pub fn subdivide_sqrt3(&mut self, iterations: usize, mut phase: bool) -> Result<bool, Error> {
        if iterations == 0 {
            return Ok(phase);
        }
        self.check_for_deleted()?;
        self.triangulate()?;
        // Phase cannot be true if any face has more than one boundary edge.
        if self.faces().any(|f| {
            self.fh_ccw_iter(f)
                .filter(|h| h.opposite().is_boundary(self))
                .count()
                > 1
        }) {
            phase = false;
        }
        let mut epoints: Vec<(A::Vector, A::Vector)> = Vec::new();
        let mut vpoints: Vec<A::Vector> = Vec::new();
        let mut fpoints: Vec<Option<A::Vector>> = Vec::new();
        // TODO: Reserve all memory required at once.
        let mut points = self.points();
        for _ in 0..iterations {
            {
                let points = points.try_borrow()?;
                // Compute edge points only if we're in phase.
                if phase {
                    epoints.resize(self.num_edges(), (A::zero_vector(), A::zero_vector()));
                    for e in self.edges().filter(|e| e.is_boundary(self)) {
                        epoints[e.index() as usize] =
                            Sqrt3Scheme::<DIM, A>::compute_edge_points(&self.topol, e, &points);
                    }
                }
                // Compute relaxed vertex positions.
                vpoints.clear();
                vpoints.extend(self.vertices().map(|v| {
                    if let Some(h) = v.halfedge(self) {
                        if h.is_boundary(self) {
                            if phase {
                                let ph = h.prev(self);
                                debug_assert!(ph.is_boundary(self));
                                ((points[h.head(self)] + points[ph.tail(self)]) * A::scalarf64(4.0)
                                    + points[v] * A::scalarf64(19.0))
                                    / A::scalarf64(27.0)
                            } else {
                                points[v]
                            }
                        } else {
                            let (valence, sum) = self
                                .vv_ccw_iter(v)
                                .fold((0usize, A::zero_vector()), |(valence, sum), nv| {
                                    (valence + 1, sum + points[nv])
                                });
                            let (a, b) = get_weight(valence);
                            sum * A::scalarf64(b) + points[v] * A::scalarf64(a)
                        }
                    } else {
                        // Isolated vertices don't move.
                        points[v]
                    }
                }));
                // Compute centroids of faces not on the boundary.
                fpoints.clear();
                fpoints.extend(self.faces().map(|f| {
                    if !phase || !f.is_boundary(self, false) {
                        Some(self.calc_face_centroid(f, &points))
                    } else {
                        None
                    }
                }));
            }
            // Make then immutable.
            let vpoints: &[A::Vector] = &vpoints;
            let epoints: &[(A::Vector, A::Vector)] = &epoints;
            let fpoints: &[Option<A::Vector>] = &fpoints;
            {
                // Update vertex positions.
                let mut points = points.try_borrow_mut()?;
                points.copy_from_slice(vpoints);
            }
            let num_old_edges = self.num_edges();
            // Split faces.
            for f in self.faces() {
                if let Some(fpoint) = fpoints[f.index() as usize] {
                    let fv = self.add_vertex(fpoint)?;
                    self.split_face(f, fv, true)?;
                } else {
                    debug_assert!(phase && !epoints.is_empty());
                    let hright = self
                        .fh_ccw_iter(f)
                        .find(|h| h.opposite().is_boundary(self))
                        .ok_or(Error::CannotSplitFace(f))?;
                    let (pa, pb) = epoints[hright.edge().index() as usize];
                    let (va, vb) = (self.add_vertex(pa)?, self.add_vertex(pb)?);
                    let hleft = self.split_edge(hright, va, true)?;
                    self.split_edge(hright, vb, true)?;
                    self.insert_edge(hleft, hleft.prev(self))?;
                    self.insert_edge(hright.next(self), hright)?;
                }
            }
            // Flip old edges.
            for e in self.edges().take(num_old_edges) {
                if e.is_boundary(self) {
                    continue;
                }
                self.swap_edge_ccw(e)?;
            }
            phase = !phase;
        }
        Ok(phase)
    }
}

#[cfg(test)]
mod test {
    use crate::{use_glam::PolyMeshF32, HasTopology};

    #[test]
    fn t_box_sqrt3_subdivision() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot make a box");
        assert_eq!(
            (8, 12, 6),
            (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces())
        );
        let mut phase = false;
        phase = mesh.subdivide_sqrt3(1, phase).expect("Cannot subdivide");
        assert_eq!(
            (20, 54, 36),
            (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces())
        );
        phase = mesh.subdivide_sqrt3(1, phase).expect("Cannot subdivide");
        assert_eq!(
            (56, 162, 108),
            (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces())
        );
        mesh.subdivide_sqrt3(1, phase).expect("Cannot subdivide");
        assert_eq!(
            (164, 486, 324),
            (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces())
        );
        assert!(mesh.check_topology().is_ok());
    }
}
