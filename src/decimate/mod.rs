pub mod edge_length;
mod heap;
pub mod quadric;

use crate::{topol::Topology, EditableTopology, Error, Handle, HasIterators, Status, HH, VH};
use heap::Heap;
use std::cmp::Ordering;

pub trait Decimater<MeshT>
where
    MeshT: HasIterators,
    Self::Cost: PartialOrd,
{
    type Cost;

    fn collapse_cost(&self, mesh: &MeshT, h: HH) -> Option<Self::Cost>;

    fn before_collapse(&mut self, mesh: &MeshT, h: HH) -> Result<(), Error>;

    fn after_collapse(&mut self, mesh: &mut MeshT, v: VH) -> Result<(), Error>;
}

fn queue_vertex_collapse<MeshT, DecT>(
    mesh: &MeshT,
    v: VH,
    module: &DecT,
    heap_pos: &mut [Option<usize>],
    heap: &mut Heap<(DecT::Cost, VH, HH)>,
) where
    MeshT: HasIterators,
    DecT: Decimater<MeshT>,
{
    if let Some((h, cost)) = mesh // Find the best collapsible edge around this vertex.
        .voh_ccw_iter(v)
        .fold(None, |best, h| {
            match (best, module.collapse_cost(mesh, h)) {
                (None, None) => None,                                   // Found nothing.
                (None, Some(cost)) => Some((h, cost)),                  // Current wins by default.
                (Some((hbest, lowest)), None) => Some((hbest, lowest)), // Previous wins by default.
                (Some((hbest, lowest)), Some(cost)) => match cost.partial_cmp(&lowest) {
                    Some(Ordering::Less) => Some((h, cost)), // Current wins by comparison
                    _ => Some((hbest, lowest)),              // Previous wins.
                },
            }
        })
    {
        // Update or insert the found edge with the new cost.
        heap_pos[v.index() as usize] = Some(match heap_pos[v.index() as usize] {
            Some(pos) => heap.update(pos, (cost, v, h)),
            None => heap.push((cost, v, h)),
        });
    } else if let Some(pos) = std::mem::replace(&mut heap_pos[v.index() as usize], None) {
        // Remove the vertex from it's previous position in the heap.
        heap.remove(pos);
    }
}

fn is_collapse_legal(mesh: &Topology, h: HH, estatus: &[Status], vstatus: &mut [Status]) -> bool {
    let v = h.tail(mesh);
    !vstatus[v.index() as usize].feature() // Vertex not a feature.
        && !estatus[h.edge().index() as usize].feature() // Edge not a feature.
        && mesh.vf_ccw_iter(v).take(2).count() == 2 // Has at leaset two incident faces.
        && mesh.check_edge_collapse(h, estatus, vstatus) // No topological errors.
}

pub trait HasDecimation: EditableTopology {
    fn decimate_while<F, DecT>(&mut self, module: &mut DecT, pred: F) -> Result<usize, Error>
    where
        F: Fn(usize, usize, usize) -> bool,
        DecT: Decimater<Self>,
    {
        let mut heap_pos = self.create_vertex_prop(None);
        let mut heap_pos = heap_pos.try_borrow_mut()?;
        let mut cache = Vec::new();
        let mut vstatus = self.vertex_status_prop();
        let mut hstatus = self.halfedge_status_prop();
        let mut estatus = self.edge_status_prop();
        let mut fstatus = self.face_status_prop();
        // Initialize heap with vertex collapses.
        let mut heap = Heap::<(DecT::Cost, VH, HH)>::new();
        {
            let vstatus = vstatus.try_borrow()?;
            for v in self
                .vertices()
                .filter(|v| !vstatus[v.index() as usize].deleted())
            {
                queue_vertex_collapse(self, v, module, &mut heap_pos, &mut heap);
            }
        }
        let (mut nc, mut nv, mut nf) = (0usize, self.num_vertices(), self.num_faces());
        let mut one_ring = Vec::new();
        while let Some((_, v0, h)) = heap.pop() {
            heap_pos[v0.index() as usize] = None;
            if !pred(nc, nv, nf) {
                break;
            }
            {
                let estatus = estatus.try_borrow()?;
                let mut vstatus = vstatus.try_borrow_mut()?;
                if !is_collapse_legal(self.topology(), h, &estatus, &mut vstatus) {
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
            let v1 = h.head(self);
            module.before_collapse(self, h)?;
            {
                let mut vstatus = vstatus.try_borrow_mut()?;
                let mut hstatus = hstatus.try_borrow_mut()?;
                let mut estatus = estatus.try_borrow_mut()?;
                let mut fstatus = fstatus.try_borrow_mut()?;
                self.collapse_edge(
                    h,
                    &mut vstatus,
                    &mut hstatus,
                    &mut estatus,
                    &mut fstatus,
                    &mut cache,
                );
            }
            module.after_collapse(self, v1)?;
            {
                // Update heap.
                for &v in one_ring.iter() {
                    queue_vertex_collapse(self, v, module, &mut heap_pos, &mut heap);
                }
            }
        }
        heap.clear();
        Ok(nc)
    }

    fn decimate(
        &mut self,
        module: &mut impl Decimater<Self>,
        num_collapses: usize,
    ) -> Result<usize, Error> {
        self.decimate_while(module, |n, _v, _f| n < num_collapses)
    }

    fn decimate_to_vertex_count(
        &mut self,
        module: &mut impl Decimater<Self>,
        vert_target: usize,
    ) -> Result<usize, Error> {
        let vtarget = usize::min(self.num_vertices(), vert_target);
        self.decimate_while(module, |_n, v, _f| v > vtarget)
    }

    fn decimate_to_face_count(
        &mut self,
        module: &mut impl Decimater<Self>,
        face_target: usize,
    ) -> Result<usize, Error> {
        let ftarget = usize::min(face_target, self.num_faces());
        self.decimate_while(module, |_n, _v, f| f > ftarget)
    }
}
