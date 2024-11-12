use std::{
    cell::{Ref, RefCell, RefMut},
    marker::PhantomData,
    rc::{Rc, Weak},
};

use crate::{
    element::{Handle, EH, FH, HH, VH},
    error::Error,
};

pub(crate) struct PropertyContainer<H>
where
    H: Handle,
{
    props: Vec<Box<dyn GenericProperty<H>>>,
    length: usize,
    _phantom: PhantomData<H>,
}

impl<H> Default for PropertyContainer<H>
where
    H: Handle,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<H> PropertyContainer<H>
where
    H: Handle,
{
    pub fn new() -> Self {
        PropertyContainer {
            props: Vec::new(),
            length: 0,
            _phantom: PhantomData,
        }
    }

    fn push_property(&mut self, prop: Box<dyn GenericProperty<H>>) {
        self.props.push(prop);
    }

    pub fn reserve(&mut self, n: usize) -> Result<(), Error> {
        for prop in self.props.iter_mut() {
            prop.reserve(n)?;
        }
        Ok(())
    }

    pub fn resize(&mut self, n: usize) -> Result<(), Error> {
        for prop in self.props.iter_mut() {
            prop.resize(n)?;
        }
        self.length = n;
        Ok(())
    }

    pub fn clear(&mut self) -> Result<(), Error> {
        for prop in self.props.iter_mut() {
            prop.clear()?;
        }
        self.length = 0;
        Ok(())
    }

    pub fn push_value(&mut self) -> Result<(), Error> {
        let (count, err) = self
            .props
            .iter_mut()
            .fold((0usize, Ok(())), |(count, err), prop| match err {
                Ok(()) => match prop.push() {
                    Ok(()) => (count + 1, Ok(())),
                    Err(e) => (count, Err(e)),
                },
                Err(e) => (count, Err(e)),
            });
        // If something went wrong, go back to how things were.
        if err.is_err() {
            for prop in self.props.iter_mut().take(count) {
                prop.resize(self.length)?;
            }
            return err;
        }
        self.length += 1;
        Ok(())
    }

    pub fn push_values(&mut self, num: usize) -> Result<(), Error> {
        let (count, err) = self
            .props
            .iter_mut()
            .fold((0usize, Ok(())), |(count, err), prop| match err {
                Ok(()) => match prop.push_many(num) {
                    Ok(()) => (count + 1, Ok(())),
                    Err(e) => (count, Err(e)),
                },
                Err(e) => (count, Err(e)),
            });
        // If something went wrong, go back to how things were.
        if err.is_err() {
            for prop in self.props.iter_mut().take(count) {
                prop.resize(self.length)?;
            }
            return err;
        }
        self.length += num;
        Ok(())
    }

    pub fn swap(&mut self, i: usize, j: usize) -> Result<(), Error> {
        for prop in self.props.iter_mut() {
            prop.swap(i, j)?;
        }
        Ok(())
    }

    pub fn copy(&mut self, src: H, dst: H) -> Result<(), Error> {
        for prop in self.props.iter_mut() {
            prop.copy(src.index() as usize, dst.index() as usize)?;
        }
        Ok(())
    }

    pub fn copy_many(&mut self, src: &[H], dst: &[H]) -> Result<(), Error> {
        for prop in self.props.iter_mut() {
            prop.copy_many(src, dst)?;
        }
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn garbage_collection(&mut self) {
        self.props.retain(|prop| prop.is_valid())
    }
}

trait GenericProperty<H>
where
    H: Handle,
{
    fn reserve(&mut self, n: usize) -> Result<(), Error>;

    fn resize(&mut self, n: usize) -> Result<(), Error>;

    fn clear(&mut self) -> Result<(), Error>;

    fn push(&mut self) -> Result<(), Error>;

    fn push_many(&mut self, num: usize) -> Result<(), Error>;

    fn swap(&mut self, i: usize, j: usize) -> Result<(), Error>;

    fn copy(&mut self, src: usize, dst: usize) -> Result<(), Error>;

    fn copy_many(&mut self, src: &[H], dst: &[H]) -> Result<(), Error>;

    fn is_valid(&self) -> bool;
}

/// This represents a property defined on the elements of the mesh. `T` is the
/// type of data associated with each element of the mesh, whose handle type is
/// `H`.
///
/// Why use properties instead of simple [`Vec<T>`] to associate values with
/// elements of a mesh? Say you use a simple [`Vec<T>`] to keep track of
/// properties. If you modify the mesh, by either adding new elements (vertices,
/// faces, etc.), or by deleting elements and garbage collecting. The [`Vec<T>`]
/// will go out of sync with the mesh. Instead using a [`Property<H, T>`]
/// guarantees the properties are always synchronized with the mesh, and that
/// every element of the mesh of type `H`, even the newly added ones will have a
/// value associated with it.
#[derive(Clone)]
pub struct Property<H, T>
where
    H: Handle,
    T: Clone + Copy,
{
    data: Rc<RefCell<Vec<T>>>,
    default: T,
    _phantom: PhantomData<H>,
}

impl<H, T> Property<H, T>
where
    H: Handle,
    T: Clone + Copy + 'static,
{
    pub(crate) fn new(container: &mut PropertyContainer<H>, default: T) -> Self {
        let prop = Property {
            data: Rc::new(RefCell::new(vec![default; container.len()])),
            default,
            _phantom: PhantomData,
        };
        container.push_property(prop.generic_ref());
        prop
    }

    pub(crate) fn with_capacity(
        n: usize,
        container: &mut PropertyContainer<H>,
        default: T,
    ) -> Self {
        let mut buf = Vec::with_capacity(n);
        buf.resize(container.len(), default);
        let prop = Property {
            data: Rc::new(RefCell::new(buf)),
            default,
            _phantom: PhantomData,
        };
        container.push_property(prop.generic_ref());
        prop
    }

    fn generic_ref(&self) -> Box<dyn GenericProperty<H>> {
        Box::new(PropertyRef {
            data: Rc::downgrade(&self.data),
            default: self.default,
        })
    }

    /// Try to borrow the property with read-only access.
    ///
    /// Properties use interior mutability pattern using a [`RefCell<T>`] to
    /// enforce runtime borrow checking rules. If borrowing fails,
    /// [`Error::BorrowedPropertyAccess`] is returned, otherwise a reference to
    /// the property is returned.
    pub fn try_borrow(&self) -> Result<Ref<[T]>, Error> {
        Ok(Ref::map(
            self.data
                .try_borrow()
                .map_err(|_| Error::BorrowedPropertyAccess)?,
            |p| -> &[T] { p },
        ))
    }

    /// Try to borrow the property with mutable access.
    ///
    /// Properties use interior mutability pattern using a [`RefCell<T>`] to
    /// enforce runtime borrow checking rules. If borrowing fails,
    /// [`Error::BorrowedPropertyAccess`] is returned, otherwise a mutable
    /// reference to the property is returned.
    pub fn try_borrow_mut(&mut self) -> Result<RefMut<[T]>, Error> {
        Ok(RefMut::map(
            self.data
                .try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?,
            |p| -> &mut [T] { p },
        ))
    }

    /// Read the property value of a mesh element.
    ///
    /// This function internally tries to borrow the property and returns an
    /// error if borrowing fails.
    pub fn get(&self, h: H) -> Result<T, Error> {
        Ok(*self
            .try_borrow()?
            .get(h.index() as usize)
            .ok_or(Error::OutOfBoundsAccess)?)
    }

    /// Get a mutable reference to the property value of a mesh element.
    ///
    /// This function internally tries to mutably borrow the property and
    /// returns an error if borrowing fails.
    pub fn get_mut(&mut self, h: H) -> Result<RefMut<T>, Error> {
        Ok(RefMut::map(
            self.data
                .try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?,
            |v| &mut v[h.index() as usize],
        ))
    }

    /// Set the property value of a mesh element.
    ///
    /// This function internally tries to mutably borrow the property and
    /// returns an error if borrowing fails.
    pub fn set(&mut self, h: H, val: T) -> Result<(), Error> {
        (*self.get_mut(h)?) = val;
        Ok(())
    }
}

/// Vertex property. A value of type `T` is defined on each vertex of the
/// mesh.
///
/// See the documentation of [`Property<H, T>`] for more context on how
/// properties work.
///
/// ```rust
/// use alum::alum_glam::PolyMeshF32;
///
/// let mut mesh = PolyMeshF32::icosahedron(1.0).expect("Cannot create an icosahedron");
/// // Crate a vertex property of type u32, with a default value of 42.
/// let vprop = mesh.create_vertex_prop(42u32);
/// let v = 2.into(); // Vertex indexed 2.
/// assert_eq!(42, vprop.get(v).expect("Cannot read vertex property"));
/// ```
pub type VProperty<T> = Property<VH, T>;

/// Halfedge property. A value of type `T` is defined on each halfedge of the
/// mesh.
///
/// See the documentation of [`Property<H, T>`] for more context on how
/// properties work.
///
/// ```rust
/// use alum::alum_glam::PolyMeshF32;
///
/// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot create an icosahedron");
/// // Crate a halfedge property of type u32, with a default value of 42.
/// let hprop = mesh.create_halfedge_prop(42u32);
/// let h = 2.into(); // Halfedge indexed 2.
/// assert_eq!(42, hprop.get(h).expect("Cannot read halfedge property"));
/// ```
pub type HProperty<T> = Property<HH, T>;

/// Edge property. A value of type `T` is defined on each edge of the
/// mesh.
///
/// See the documentation of [`Property<H, T>`] for more context on how
/// properties work.
///
/// ```rust
/// use alum::alum_glam::PolyMeshF32;
///
/// let mut mesh = PolyMeshF32::octahedron(1.0).expect("Cannot create an icosahedron");
/// // Crate a edge property of type u32, with a default value of 42.
/// let eprop = mesh.create_edge_prop(42u32);
/// let e = 2.into(); // Edge indexed 2.
/// assert_eq!(42, eprop.get(e).expect("Cannot read edge property"));
/// ```
pub type EProperty<T> = Property<EH, T>;

/// Face property. A value of type `T` is defined on each face of the
/// mesh.
///
/// See the documentation of [`Property<H, T>`] for more context on how
/// properties work.
///
/// ```rust
/// use alum::alum_glam::PolyMeshF32;
///
/// let mut mesh = PolyMeshF32::dodecahedron(1.0).expect("Cannot create an icosahedron");
/// // Crate a face property of type u32, with a default value of 42.
/// let fprop = mesh.create_face_prop(42u32);
/// let f = 2.into(); // Face indexed 2.
/// assert_eq!(42, fprop.get(f).expect("Cannot read face property"));
/// ```
pub type FProperty<T> = Property<FH, T>;

struct PropertyRef<T>
where
    T: Clone + Copy,
{
    data: Weak<RefCell<Vec<T>>>,
    default: T,
}

impl<H, T> GenericProperty<H> for PropertyRef<T>
where
    T: Clone + Copy,
    H: Handle,
{
    fn reserve(&mut self, n: usize) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            prop.try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)? // reserve memory.
                .reserve(n);
        }
        Ok(())
    }

    fn resize(&mut self, n: usize) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            prop.try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?
                .resize(n, self.default);
        }
        Ok(())
    }

    fn clear(&mut self) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            prop.try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?
                .clear();
        }
        Ok(())
    }

    fn push(&mut self) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            prop.try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?
                .push(self.default);
        }
        Ok(())
    }

    fn push_many(&mut self, num: usize) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            let mut prop = prop
                .try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?;
            let prop: &mut Vec<T> = &mut prop;
            prop.resize(prop.len() + num, self.default);
        }
        Ok(())
    }

    fn swap(&mut self, i: usize, j: usize) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            prop.try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?
                .swap(i, j);
        }
        Ok(())
    }

    fn copy(&mut self, src: usize, dst: usize) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            let mut buf = prop
                .try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?;
            let buf: &mut [T] = &mut buf;
            buf[dst] = buf[src];
        }
        Ok(())
    }

    fn copy_many(&mut self, src: &[H], dst: &[H]) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            if src.len() != dst.len() {
                return Err(Error::MismatchedArrayLengths(src.len(), dst.len()));
            }
            let mut buf = prop
                .try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?;
            let buf: &mut [T] = &mut buf;
            for (src, dst) in src
                .iter()
                .map(|h| h.index() as usize)
                .zip(dst.iter().map(|h| h.index() as usize))
            {
                buf[dst] = buf[src];
            }
        }
        Ok(())
    }

    fn is_valid(&self) -> bool {
        self.data.upgrade().is_some()
    }
}

#[cfg(test)]
mod test {
    use crate::topol::{test::quad_box, TopolCache};

    use super::{PropertyContainer, VProperty};

    #[test]
    fn t_garbage_collection() {
        let mut container = PropertyContainer::new();
        assert_eq!(container.props.len(), 0);
        {
            let _prop0 = VProperty::<u32>::new(&mut container, 0);
            assert_eq!(container.props.len(), 1);
            {
                let _prop1 = VProperty::<u16>::new(&mut container, 0);
                assert_eq!(container.props.len(), 2);
            }
            assert_eq!(container.props.len(), 2);
            assert_eq!(
                container
                    .props
                    .iter()
                    .filter(|prop| prop.is_valid())
                    .count(),
                1
            );
            container.garbage_collection();
            assert_eq!(container.props.len(), 1);
        }
        assert_eq!(container.props.len(), 1);
        assert_eq!(
            container
                .props
                .iter()
                .filter(|prop| prop.is_valid())
                .count(),
            0
        );
        container.garbage_collection();
        assert_eq!(container.props.len(), 0);
        let mut _prop = Some(VProperty::<u8>::new(&mut container, 0));
        assert_eq!(container.props.len(), 1);
        _prop = None;
        assert_eq!(container.props.len(), 1);
        assert!(!container.props[0].is_valid());
        container.garbage_collection();
        assert_eq!(container.props.len(), 0);
    }

    #[test]
    fn t_deleted_elements() {
        let mut qbox = quad_box();
        let (fis, vis) = {
            let mut vis = qbox.new_vprop(0u32);
            let mut fis = qbox.new_fprop(0u32);
            {
                let mut vis = vis.try_borrow_mut().expect("Cannot borrow property");
                for (i, v) in vis.iter_mut().enumerate() {
                    *v = i as u32;
                }
                let mut fis = fis.try_borrow_mut().expect("Cannot borrow property");
                for (i, v) in fis.iter_mut().enumerate() {
                    *v = i as u32;
                }
            }
            (fis, vis)
        };
        let mut cache = TopolCache::default();
        for fi in [1, 2, 5] {
            qbox.delete_face(
                fi.into(),
                true,
                &mut cache.halfedges,
                &mut cache.edges,
                &mut cache.vertices,
            )
            .expect("Cannot delete a face");
        }
        qbox.garbage_collection(&mut cache)
            .expect("Cannot garbage collect");
        let mut fis: Vec<u32> = fis
            .try_borrow()
            .expect("Cannot borrow face indices")
            .to_vec();
        fis.sort();
        let mut vis: Vec<u32> = vis
            .try_borrow()
            .expect("Cannot borrow vertex indices")
            .to_vec();
        vis.sort();
        // The properties should be preserved.
        assert_eq!(vis, &[0, 1, 2, 3, 4, 6, 7]);
        assert_eq!(fis, &[0, 3, 4]);
    }
}
