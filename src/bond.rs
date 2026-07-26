use crate::atom::Atom;

#[derive(Debug)]
pub struct Bond<'a> {
    atom1: &'a Atom,
    atom2: &'a Atom,
}
