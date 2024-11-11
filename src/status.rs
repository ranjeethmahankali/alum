const DELETED: u8 = 1 << 0;
const LOCKED: u8 = 1 << 1;
const SELECTED: u8 = 1 << 2;
const HIDDEN: u8 = 1 << 3;
const FEATURE: u8 = 1 << 4;
const TAGGED: u8 = 1 << 5;
const TAGGED2: u8 = 1 << 6;
const UNUSED: u8 = 1 << 7;

/// Status of a mesh element.
///
/// This can be used to keep track of whether an element is deleted, marked as a
/// feature, tagged, marked hidden etc. This is used internally to keep track of
/// deletion, and garbage collection.
#[derive(Clone, Copy)]
pub struct Status {
    flags: u8,
}

impl Default for Status {
    /// New status with none of the flags set.
    fn default() -> Self {
        Self { flags: 0 }
    }
}

impl Status {
    fn check(&self, i: u8) -> bool {
        self.flags & i > 0
    }

    fn set(&mut self, i: u8, flag: bool) {
        if flag {
            self.flags |= i;
        } else {
            self.flags &= !i;
        }
    }

    /// Check the deleted flag.
    pub fn deleted(&self) -> bool {
        self.check(DELETED)
    }

    /// Set the deleted flag.
    pub fn set_deleted(&mut self, flag: bool) {
        self.set(DELETED, flag);
    }

    /// Check the locked flag.
    pub fn locked(&self) -> bool {
        self.check(LOCKED)
    }

    /// Set the locked flag.
    pub fn set_locked(&mut self, flag: bool) {
        self.set(LOCKED, flag)
    }

    /// Check the selected flag.
    pub fn selected(&self) -> bool {
        self.check(SELECTED)
    }

    /// Set the selected flag.
    pub fn set_selected(&mut self, flag: bool) {
        self.set(SELECTED, flag)
    }

    /// Check the hidden flag.
    pub fn hidden(&self) -> bool {
        self.check(HIDDEN)
    }

    /// Set the hidden flag.
    pub fn set_hidden(&mut self, flag: bool) {
        self.set(HIDDEN, flag)
    }

    /// Check the feature flag.
    pub fn feature(&self) -> bool {
        self.check(FEATURE)
    }

    /// Set the feature flag.
    pub fn set_feature(&mut self, flag: bool) {
        self.set(FEATURE, flag)
    }

    /// Check the tagged flag.
    pub fn tagged(&self) -> bool {
        self.check(TAGGED)
    }

    /// Set the tagged flag.
    pub fn set_tagged(&mut self, flag: bool) {
        self.set(TAGGED, flag)
    }

    /// Check the tagged2 flag.
    pub fn tagged2(&self) -> bool {
        self.check(TAGGED2)
    }

    /// Set the tagged2 flag.
    pub fn set_tagged2(&mut self, flag: bool) {
        self.set(TAGGED2, flag)
    }

    /// Check the unused flag.
    pub fn unused(&self) -> bool {
        self.check(UNUSED)
    }

    /// Set the unused flag.
    pub fn set_unused(&mut self, flag: bool) {
        self.set(UNUSED, flag)
    }
}
