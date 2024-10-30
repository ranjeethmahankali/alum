use std::path::Path;

use crate::{
    element::VH,
    error::Error,
    mesh::PolyMeshT,
    vector::{TScalar, TVec3},
};

impl<VecT: TVec3> PolyMeshT<VecT>
where
    VecT::Scalar: TScalar,
{
    pub fn load_obj(path: &Path) -> Result<Self, Error> {
        let options = tobj::LoadOptions::default();
        let (models, _) =
            tobj::load_obj(path, &options).map_err(|e| Error::ObjLoadFailed(format!("{}", e)))?;
        let (nverts, nfaces) = models
            .iter()
            .fold((0usize, 0usize), |(nverts, nfaces), model| {
                let msh = &model.mesh;
                (
                    nverts + (msh.positions.len() / 3),
                    nfaces + msh.face_arities.len(),
                )
            });
        let nedges = nfaces * 3 / 2; // Estimate.
        let mut outmesh = PolyMeshT::<VecT>::with_capacity(nverts, nedges, nfaces);
        let mut voffset = 0u32;
        let mut positions = Vec::new();
        let mut vertices = Vec::new();
        for model in models {
            let mesh = model.mesh;
            if mesh.positions.len() % 3 != 0 {
                return Err(Error::IncorrectNumberOfCoordinates(mesh.positions.len()));
            }
            let nverts = (mesh.positions.len() / 3) as u32;
            positions.clear();
            positions.extend(mesh.positions.chunks(3).map(|triplet| {
                VecT::new(
                    VecT::Scalar::from_f32(triplet[0]),
                    VecT::Scalar::from_f32(triplet[1]),
                    VecT::Scalar::from_f32(triplet[2]),
                )
            }));
            vertices.resize(positions.len(), 0u32.into());
            outmesh.add_vertices(&positions, &mut vertices)?;
            // TODO: Load all faces at once for faster performance.
            let mut start = 0usize;
            let mut fvs = Vec::new();
            for size in mesh.face_arities {
                let size = size as usize;
                let indices = &mesh.indices[start..(start + size)];
                start += size;
                fvs.clear();
                fvs.extend(indices.iter().map(|i| {
                    let i = i + voffset;
                    let v: VH = i.into();
                    v
                }));
                outmesh.add_face(&fvs)?;
            }
            voffset += nverts;
        }
        Ok(outmesh)
    }
}

#[cfg(test)]
mod test {
    use std::{path::PathBuf, str::FromStr, time::Instant};

    use crate::mesh::PolyMeshF32;

    #[test]
    fn t_temp() {
        let before = Instant::now();
        PolyMeshF32::load_obj(
            &PathBuf::from_str("/home/rnjth94/dev/galproject/assets/bunny.obj").unwrap(),
        )
        .unwrap();
        let duration = Instant::now() - before;
        println!("{}ms", duration.as_millis());
        todo!()
    }
}
