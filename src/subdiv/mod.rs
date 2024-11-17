mod catmull;
mod loop_subd;

use crate::{topol::Topology, Error, Handle, HasTopology};

fn check_for_deleted(topol: &Topology) -> Result<(), Error> {
    // Ensure no deleted topology.
    let fstatus = topol.fstatus.try_borrow()?;
    if let Some(f) = topol
        .faces()
        .find(|f| fstatus[f.index() as usize].deleted())
    {
        return Err(Error::DeletedFace(f));
    }
    let estatus = topol.estatus.try_borrow()?;
    if let Some(e) = topol
        .edges()
        .find(|e| estatus[e.index() as usize].deleted())
    {
        return Err(Error::DeletedEdge(e));
    }
    let vstatus = topol.hstatus.try_borrow()?;
    if let Some(v) = topol
        .vertices()
        .find(|v| vstatus[v.index() as usize].deleted())
    {
        return Err(Error::DeletedVertex(v));
    }
    Ok(())
}
