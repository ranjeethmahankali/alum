use std::path::Path;

use crate::{
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
        let mut positions = Vec::new();
        let mut vertices = Vec::new();
        let mut fvs = Vec::new();
        for model in models {
            let mesh = &model.mesh;
            if mesh.positions.len() % 3 != 0 {
                return Err(Error::IncorrectNumberOfCoordinates(mesh.positions.len()));
            }
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
            // Faces.
            let mut start = 0usize;
            if mesh.face_arities.is_empty() {
                // Load triangles.
                for indices in mesh.indices.chunks(3) {
                    fvs.clear();
                    fvs.extend(indices.iter().map(|i| vertices[*i as usize]));
                    outmesh.add_face(&fvs)?;
                }
            } else {
                for size in &mesh.face_arities {
                    let size = *size as usize;
                    let indices = &mesh.indices[start..(start + size)];
                    start += size;
                    fvs.clear();
                    fvs.extend(indices.iter().map(|i| vertices[*i as usize]));
                    outmesh.add_face(&fvs)?;
                }
            }
        }
        Ok(outmesh)
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::mesh::PolyMeshF32;

    fn bunny_mesh() -> PolyMeshF32 {
        let path = {
            let mut dirpath = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
            dirpath.push("assets");
            dirpath.push("bunny.obj");
            dirpath
        };
        dbg!(&path);
        PolyMeshF32::load_obj(&path).expect("Cannot load mesh")
    }

    #[test]
    fn t_bunny_topol() {
        let mesh = bunny_mesh();
        assert_eq!(2503, mesh.num_vertices());
        assert_eq!(7473, mesh.num_edges());
        assert_eq!(4968, mesh.num_faces());
        assert_eq!(
            42,
            mesh.edges().filter(|e| mesh.is_boundary_edge(*e)).count()
        );
    }

    #[test]
    fn t_bunny_area() {
        let mesh = bunny_mesh();
        assert_eq!(
            mesh.try_calc_area().expect("Unable to compute the area"),
            0.05646857
        );
    }
}
