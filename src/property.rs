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

// 'static lifetime enforces the data stored inside properties is fully owned
// and doesn't contain any weird references.
pub trait TPropData: Default + Clone + Copy + 'static {}

// Unsigned integers.
impl TPropData for u8 {}
impl TPropData for u16 {}
impl TPropData for u32 {}
impl TPropData for u64 {}
impl TPropData for u128 {}
impl TPropData for usize {}
// Signed integers.
impl TPropData for i8 {}
impl TPropData for i16 {}
impl TPropData for i32 {}
impl TPropData for i64 {}
impl TPropData for i128 {}
impl TPropData for isize {}
// Floating point types.
impl TPropData for f32 {}
impl TPropData for f64 {}
// Other types.
impl TPropData for bool {}
impl TPropData for char {}
impl TPropData for glam::Vec3 {}

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

#[derive(Clone)]
pub struct Property<H: Handle, T: TPropData> {
    data: Rc<RefCell<Vec<T>>>,
    _phantom: PhantomData<H>,
}

impl<H: Handle, T: TPropData> Property<H, T> {
    pub(crate) fn new(container: &mut PropertyContainer<H>) -> Self {
        let prop = Property {
            data: Rc::new(RefCell::new(vec![T::default(); container.len()])),
            _phantom: PhantomData,
        };
        container.push_property(prop.generic_ref());
        prop
    }

    pub(crate) fn with_capacity(n: usize, container: &mut PropertyContainer<H>) -> Self {
        let mut buf = Vec::with_capacity(n);
        buf.resize(container.len(), T::default());
        let prop = Property {
            data: Rc::new(RefCell::new(buf)),
            _phantom: PhantomData,
        };
        container.push_property(prop.generic_ref());
        prop
    }

    fn generic_ref(&self) -> Box<dyn GenericProperty<H>> {
        Box::new(PropertyRef {
            data: Rc::downgrade(&self.data),
        })
    }

    pub fn try_borrow(&self) -> Result<Ref<Vec<T>>, Error> {
        self.data
            .try_borrow()
            .map_err(|_| Error::BorrowedPropertyAccess)
    }

    pub fn try_borrow_mut(&mut self) -> Result<RefMut<Vec<T>>, Error> {
        self.data
            .try_borrow_mut()
            .map_err(|_| Error::BorrowedPropertyAccess)
    }

    pub fn get(&self, i: H) -> Result<T, Error> {
        Ok(*self
            .try_borrow()?
            .get(i.index() as usize)
            .ok_or(Error::OutOfBoundsAccess)?)
    }

    pub fn get_mut(&mut self, i: H) -> Result<RefMut<T>, Error> {
        Ok(RefMut::map(self.try_borrow_mut()?, |v| {
            &mut v[i.index() as usize]
        }))
    }

    pub fn set(&mut self, i: H, val: T) -> Result<(), Error> {
        (*self.get_mut(i)?) = val;
        Ok(())
    }
}

impl<H: Handle, T: TPropData> Default for Property<H, T> {
    fn default() -> Self {
        Self {
            data: Default::default(),
            _phantom: PhantomData,
        }
    }
}

pub type VProperty<T> = Property<VH, T>;
pub type HProperty<T> = Property<HH, T>;
pub type EProperty<T> = Property<EH, T>;
pub type FProperty<T> = Property<FH, T>;

struct PropertyRef<T: TPropData> {
    data: Weak<RefCell<Vec<T>>>,
}

impl<T: TPropData, H: Handle> GenericProperty<H> for PropertyRef<T> {
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
                .resize(n, T::default());
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
                .push(T::default());
        }
        Ok(())
    }

    fn push_many(&mut self, num: usize) -> Result<(), Error> {
        if let Some(prop) = self.data.upgrade() {
            let mut prop = prop
                .try_borrow_mut()
                .map_err(|_| Error::BorrowedPropertyAccess)?;
            let prop: &mut Vec<T> = &mut prop;
            prop.resize(prop.len() + num, T::default());
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
    use super::{PropertyContainer, VProperty};

    #[test]
    fn t_garbage_collection() {
        let mut container = PropertyContainer::new();
        assert_eq!(container.props.len(), 0);
        {
            let _prop0 = VProperty::<u32>::new(&mut container);
            assert_eq!(container.props.len(), 1);
            {
                let _prop1 = VProperty::<u16>::new(&mut container);
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
        let mut _prop = Some(VProperty::<u8>::new(&mut container));
        assert_eq!(container.props.len(), 1);
        _prop = None;
        assert_eq!(container.props.len(), 1);
        assert!(!container.props[0].is_valid());
        container.garbage_collection();
        assert_eq!(container.props.len(), 0);
    }
}
