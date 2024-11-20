mod catmull;
mod loop_subd;
mod sqrt3;

use crate::{topol::Topology, Error, HasTopology};

fn check_for_deleted(topol: &Topology) -> Result<(), Error> {
    // Ensure no deleted topology.
    let fstatus = topol.fstatus.try_borrow()?;
    if let Some(f) = topol.faces().find(|f| fstatus[*f].deleted()) {
        return Err(Error::DeletedFace(f));
    }
    let estatus = topol.estatus.try_borrow()?;
    if let Some(e) = topol.edges().find(|e| estatus[*e].deleted()) {
        return Err(Error::DeletedEdge(e));
    }
    let vstatus = topol.vstatus.try_borrow()?;
    if let Some(v) = topol.vertices().find(|v| vstatus[*v].deleted()) {
        return Err(Error::DeletedVertex(v));
    }
    Ok(())
}
