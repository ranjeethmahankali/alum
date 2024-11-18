use crate::{
    topol::Topology, Adaptor, EditableTopology, Error, Handle, HasIterators, HasTopology,
    PolyMeshT, Status, HH, VH,
};
use heap::Heap;

mod heap;

pub trait DecimationModule {
    fn collapse_cost(&self, mesh: &Topology, h: HH) -> Option<f64>;

    fn before_collapse(&self, mesh: &Topology, h: HH);

    fn after_collapse(&self, mesh: &Topology, h: HH);
}

fn queue_vertex_collapse(
    mesh: &Topology,
    v: VH,
    module: &impl DecimationModule,
    heap: &mut Heap<(f64, VH, HH)>,
) {
    mesh.voh_ccw_iter(v)
        .fold((None, f64::MAX), |(hopt, best_cost), h| {
            match module.collapse_cost(mesh, h) {
                Some(cost) if cost < best_cost => (Some(h), cost),
                _ => (hopt, best_cost),
            }
        });
    todo!()
}

fn is_collapse_legal(mesh: &Topology, h: HH, estatus: &[Status], vstatus: &mut [Status]) -> bool {
    // Check for:
    // - Vertex and edge are not features.
    // - At least two incident vertices at the tail vertex.
    // - No topological errors.
    let v = h.tail(mesh);
    !vstatus[v.index() as usize].feature()
        && mesh.vf_ccw_iter(v).take(2).count() == 2
        && mesh.check_edge_collapse(h, estatus, vstatus)
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    fn decimate_while<M, F>(
        &mut self,
        module: &impl DecimationModule,
        pred: F,
    ) -> Result<usize, Error>
    where
        F: Fn(usize, usize, usize) -> bool,
        M: DecimationModule,
    {
        // Initialize heap with vertex collapses.
        let mut heap = Heap::<(f64, VH, HH)>::new();
        {
            let status = self.topol.vstatus.try_borrow()?;
            for v in self
                .vertices()
                .filter(|v| !status[v.index() as usize].deleted())
            {
                queue_vertex_collapse(&self.topol, v, module, &mut heap);
            }
        }
        let mut vstatus = self.topol.vstatus.clone();
        let mut hstatus = self.topol.hstatus.clone();
        let mut estatus = self.topol.estatus.clone();
        let mut fstatus = self.topol.fstatus.clone();
        let (mut nc, mut nv, mut nf) = (0usize, self.num_vertices(), self.num_faces());
        let mut one_ring = Vec::new();
        while let Some((_, v0, h)) = heap.pop() {
            if !pred(nc, nv, nf) {
                break;
            }
            {
                let estatus = estatus.try_borrow()?;
                let mut vstatus = vstatus.try_borrow_mut()?;
                if !is_collapse_legal(&self.topol, h, &estatus, &mut vstatus) {
                    continue;
                }
            }
            one_ring.clear(); // Collect the one ring for later.
            one_ring.extend(self.vv_ccw_iter(v0));
            // Update stats.
            (nc, nv, nf) = (
                nc + 1,
                nv - 1,
                nf - if h.edge().is_boundary(self) { 1 } else { 2 },
            );
            module.before_collapse(&self.topol, h);
            {
                let mut vstatus = vstatus.try_borrow_mut()?;
                let mut hstatus = hstatus.try_borrow_mut()?;
                let mut estatus = estatus.try_borrow_mut()?;
                let mut fstatus = fstatus.try_borrow_mut()?;
                self.topol.collapse_edge(
                    h,
                    &mut vstatus,
                    &mut hstatus,
                    &mut estatus,
                    &mut fstatus,
                    &mut self.cache,
                );
            }
            module.after_collapse(&self.topol, h);
            {
                // Update heap.
                for &v in one_ring.iter() {
                    queue_vertex_collapse(&self.topol, v, module, &mut heap);
                }
            }
        }
        heap.clear();
        Ok(nc)
    }

    pub fn decimate<M>(
        &mut self,
        module: &impl DecimationModule,
        num_collapses: usize,
    ) -> Result<usize, Error>
    where
        M: DecimationModule,
    {
        self.decimate_while::<M, _>(module, |n, _v, _f| n < num_collapses)
    }

    pub fn decimate_to_vertex_count<M>(
        &mut self,
        module: &impl DecimationModule,
        vert_target: usize,
    ) -> Result<usize, Error>
    where
        M: DecimationModule,
    {
        let vtarget = usize::min(self.num_vertices(), vert_target);
        self.decimate_while::<M, _>(module, |_n, v, _f| v > vtarget)
    }

    pub fn decimate_to_face_count<M>(
        &mut self,
        module: &impl DecimationModule,
        face_target: usize,
    ) -> Result<usize, Error>
    where
        M: DecimationModule,
    {
        let ftarget = usize::min(face_target, self.num_faces());
        self.decimate_while::<M, _>(module, |_n, _v, f| f > ftarget)
    }
}
