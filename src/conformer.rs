use crate::{frame::Frame, structure::Structure};

#[derive(Clone, Copy)]
pub struct Conformer<'a> {
    structure: &'a Structure,
    frame: &'a Frame,
}

#[derive(Debug)]
pub struct FrameMismatch {
    pub atoms: usize,
    pub positions: usize,
    pub occupancies: usize,
    pub b_factors: usize,
}

impl<'a> Conformer<'a> {
    pub fn new(structure: &'a Structure, frame: &'a Frame) -> Result<Self, FrameMismatch> {
        let atoms = structure.atoms.len();
        if frame.positions.len() != atoms || frame.occupancies.len() != atoms || frame.b_factors.len() != atoms {
            return Err(FrameMismatch {
                atoms,
                positions: frame.positions.len(),
                occupancies: frame.occupancies.len(),
                b_factors: frame.b_factors.len(),
            });
        }
        Ok(Self { structure, frame })
    }

    pub fn structure(&self) -> &'a Structure {
        self.structure
    }

    pub fn frame(&self) -> &'a Frame {
        self.frame
    }
}
