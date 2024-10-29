use crate::property::TPropData;

const DELETED: u8 = 1 << 0;
const LOCKED: u8 = 1 << 1;
const SELECTED: u8 = 1 << 2;
const HIDDEN: u8 = 1 << 3;
const FEATURE: u8 = 1 << 4;
const TAGGED: u8 = 1 << 5;
const TAGGED2: u8 = 1 << 6;
const UNUSED: u8 = 1 << 7;

#[derive(Clone, Copy, Default)]
pub struct Status {
    flags: u8,
}

impl TPropData for Status {}

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

    pub fn deleted(&self) -> bool {
        self.check(DELETED)
    }

    pub fn set_deleted(&mut self, flag: bool) {
        self.set(DELETED, flag);
    }

    pub fn locked(&self) -> bool {
        self.check(LOCKED)
    }

    pub fn set_locked(&mut self, flag: bool) {
        self.set(LOCKED, flag)
    }

    pub fn selected(&self) -> bool {
        self.check(SELECTED)
    }

    pub fn set_selected(&mut self, flag: bool) {
        self.set(SELECTED, flag)
    }

    pub fn hidden(&self) -> bool {
        self.check(HIDDEN)
    }

    pub fn set_hidden(&mut self, flag: bool) {
        self.set(HIDDEN, flag)
    }

    pub fn feature(&self) -> bool {
        self.check(FEATURE)
    }

    pub fn set_feature(&mut self, flag: bool) {
        self.set(FEATURE, flag)
    }

    pub fn tagged(&self) -> bool {
        self.check(TAGGED)
    }

    pub fn set_tagged(&mut self, flag: bool) {
        self.set(TAGGED, flag)
    }

    pub fn tagged2(&self) -> bool {
        self.check(TAGGED2)
    }

    pub fn set_tagged2(&mut self, flag: bool) {
        self.set(TAGGED2, flag)
    }

    pub fn unused(&self) -> bool {
        self.check(UNUSED)
    }

    pub fn set_unused(&mut self, flag: bool) {
        self.set(UNUSED, flag)
    }
}
