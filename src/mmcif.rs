use std::collections::HashMap;
use std::fmt;
use std::ops::Range;

use crate::atom::{Atom, Element};
use crate::conformer::Conformer;
use crate::molecule::Molecule;
use crate::position::Position;
use crate::residue::Residue;
use crate::topology::Topology;
use crate::{Cif, Row, Value};

#[derive(Debug)]
pub enum MmCifError {
    NoDataBlock,
    MissingField(&'static str),
    WrongType(&'static str),
}

impl fmt::Display for MmCifError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MmCifError::NoDataBlock => write!(f, "CIF file has no data blocks"),
            MmCifError::MissingField(tag) => write!(f, "missing field {tag}"),
            MmCifError::WrongType(tag) => write!(f, "field {tag} has an unexpected value type"),
        }
    }
}

impl std::error::Error for MmCifError {}

fn field<'a, T>(row: &Row<'a>, tag: &'static str, cast: fn(&'a Value) -> Option<T>) -> Result<T, MmCifError> {
    match row.get(tag) {
        None => Err(MmCifError::MissingField(tag)),
        Some(v) => cast(v).ok_or(MmCifError::WrongType(tag)),
    }
}

/// Like `field`, but a CIF `.`/`?` (inapplicable/unknown) is `Ok(None)` rather than
/// an error — e.g. `_struct_conn` partner seq_id, which is inapplicable for
/// non-polymer partners (a bound ion has no sequence position).
fn optional_field<'a, T>(
    row: &Row<'a>,
    tag: &'static str,
    cast: fn(&'a Value) -> Option<T>,
) -> Result<Option<T>, MmCifError> {
    match row.get(tag) {
        None => Err(MmCifError::MissingField(tag)),
        Some(Value::Inapplicable) | Some(Value::Unknown) => Ok(None),
        Some(v) => cast(v).map(Some).ok_or(MmCifError::WrongType(tag)),
    }
}

#[derive(Debug)]
pub enum EntityKind {
    Polymer,
    NonPolymer,
    Branched,
    Water,
}

fn entity_kind(v: &Value) -> Option<EntityKind> {
    match v.as_str()? {
        "polymer" => Some(EntityKind::Polymer),
        "non-polymer" => Some(EntityKind::NonPolymer),
        "branched" => Some(EntityKind::Branched),
        "water" => Some(EntityKind::Water),
        _ => None,
    }
}

#[derive(Debug)]
pub struct Entity {
    pub id: u8,
    pub kind: EntityKind,
    pub name: String,
}

#[derive(Debug)]
pub struct StructAsym {
    pub id: String,
    pub entity_id: u8,
}

fn element(v: &Value) -> Option<Element> {
    // `_atom_site.type_symbol` is conventionally all-uppercase ("FE", "ZN"),
    // not `Element::from_symbol`'s canonical mixed-case form ("Fe", "Zn").
    let symbol = v.as_str()?;
    let mut chars = symbol.chars();
    let normalized = match chars.next() {
        Some(first) => first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase(),
        None => return None,
    };
    Element::from_symbol(&normalized)
}

fn is_hetero(v: &Value) -> Option<bool> {
    match v.as_str()? {
        "ATOM" => Some(false),
        "HETATM" => Some(true),
        _ => None,
    }
}

#[derive(Debug)]
pub struct AtomSite {
    pub id: u32,
    pub label_entity_id: u8,
    pub label_asym_id: String,
    pub label_atom_id: String,
    // We use auth_seq_id because for non-polymers the label is `.`, which will end up merging all the atoms into one
    // residue
    pub auth_seq_id: i64,
    pub label_comp_id: String,
    pub type_symbol: Element,
    pub is_hetero: bool,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f32,
    pub b_factor: f32,
}

#[derive(Debug)]
pub struct ChemCompBond {
    pub comp_id: String,
    pub atom_id_1: String,
    pub atom_id_2: String,
}

#[derive(Debug)]
pub struct StructConn {
    pub ptnr1_asym_id: String,
    pub ptnr1_seq_id: Option<i64>,
    pub ptnr1_atom_id: String,
    pub ptnr2_asym_id: String,
    pub ptnr2_seq_id: Option<i64>,
    pub ptnr2_atom_id: String,
}

#[derive(Debug)]
pub struct MmCIF {
    pub entities: HashMap<u8, Entity>,
    pub struct_asyms: HashMap<String, StructAsym>,
    pub atom_sites: HashMap<String, Vec<AtomSite>>,
    pub chem_comp_bonds: HashMap<String, Vec<ChemCompBond>>,
    pub struct_conns: Vec<StructConn>,
}

impl MmCIF {
    pub fn new(cif: Cif) -> Result<Self, MmCifError> {
        let block = cif.blocks.into_iter().next().ok_or(MmCifError::NoDataBlock)?;

        let mut entities = HashMap::new();
        if let Some(loop_) = block.loop_by_category("_entity") {
            for row in loop_.iter() {
                let id = field(&row, "_entity.id", Value::as_int)? as u8;
                let kind = field(&row, "_entity.type", entity_kind)?;
                let name = field(&row, "_entity.pdbx_description", Value::as_str)?.to_string();
                entities.insert(id, Entity { id, kind, name });
            }
        }

        let mut struct_asyms = HashMap::new();
        if let Some(loop_) = block.loop_by_category("_struct_asym") {
            for row in loop_.iter() {
                let id = field(&row, "_struct_asym.id", Value::as_str)?.to_string();
                let entity_id = field(&row, "_struct_asym.entity_id", Value::as_int)? as u8;
                struct_asyms.insert(id.clone(), StructAsym { id, entity_id });
            }
        }

        let mut atom_sites: HashMap<String, Vec<AtomSite>> = HashMap::new();
        if let Some(loop_) = block.loop_by_category("_atom_site") {
            for row in loop_.iter() {
                let id = field(&row, "_atom_site.id", Value::as_int)? as u32;
                let label_entity_id = field(&row, "_atom_site.label_entity_id", Value::as_int)? as u8;
                let label_asym_id = field(&row, "_atom_site.label_asym_id", Value::as_str)?.to_string();
                let label_atom_id = field(&row, "_atom_site.label_atom_id", Value::as_str)?.to_string();
                let auth_seq_id = field(&row, "_atom_site.auth_seq_id", Value::as_int)?;
                let label_comp_id = field(&row, "_atom_site.label_comp_id", Value::as_str)?.to_string();
                let type_symbol = field(&row, "_atom_site.type_symbol", element)?;
                let is_hetero_atom = field(&row, "_atom_site.group_PDB", is_hetero)?;
                let x = field(&row, "_atom_site.Cartn_x", Value::as_float)?;
                let y = field(&row, "_atom_site.Cartn_y", Value::as_float)?;
                let z = field(&row, "_atom_site.Cartn_z", Value::as_float)?;
                let occupancy = field(&row, "_atom_site.occupancy", Value::as_float)? as f32;
                let b_factor = field(&row, "_atom_site.B_iso_or_equiv", Value::as_float)? as f32;
                atom_sites.entry(label_asym_id.clone()).or_default().push(AtomSite {
                    id,
                    label_entity_id,
                    label_asym_id,
                    label_atom_id,
                    auth_seq_id,
                    label_comp_id,
                    type_symbol,
                    is_hetero: is_hetero_atom,
                    x,
                    y,
                    z,
                    occupancy,
                    b_factor,
                });
            }
        }

        let mut chem_comp_bonds: HashMap<String, Vec<ChemCompBond>> = HashMap::new();
        if let Some(loop_) = block.loop_by_category("_chem_comp_bond") {
            for row in loop_.iter() {
                let comp_id = field(&row, "_chem_comp_bond.comp_id", Value::as_str)?.to_string();
                let atom_id_1 = field(&row, "_chem_comp_bond.atom_id_1", Value::as_str)?.to_string();
                let atom_id_2 = field(&row, "_chem_comp_bond.atom_id_2", Value::as_str)?.to_string();
                chem_comp_bonds.entry(comp_id.clone()).or_default().push(ChemCompBond {
                    comp_id,
                    atom_id_1,
                    atom_id_2,
                });
            }
        }

        let mut struct_conns = Vec::new();
        if let Some(loop_) = block.loop_by_category("_struct_conn") {
            for row in loop_.iter() {
                let ptnr1_symmetry = field(&row, "_struct_conn.ptnr1_symmetry", Value::as_str)?;
                let ptnr2_symmetry = field(&row, "_struct_conn.ptnr2_symmetry", Value::as_str)?;
                if ptnr1_symmetry != "1_555" || ptnr2_symmetry != "1_555" {
                    // Partner lives in a symmetry-generated copy, not this asymmetric
                    // unit's own atoms — out of scope for now.
                    continue;
                }

                let ptnr1_asym_id = field(&row, "_struct_conn.ptnr1_label_asym_id", Value::as_str)?.to_string();
                let ptnr1_seq_id = optional_field(&row, "_struct_conn.ptnr1_label_seq_id", Value::as_int)?;
                let ptnr1_atom_id = field(&row, "_struct_conn.ptnr1_label_atom_id", Value::as_str)?.to_string();
                let ptnr2_asym_id = field(&row, "_struct_conn.ptnr2_label_asym_id", Value::as_str)?.to_string();
                let ptnr2_seq_id = optional_field(&row, "_struct_conn.ptnr2_label_seq_id", Value::as_int)?;
                let ptnr2_atom_id = field(&row, "_struct_conn.ptnr2_label_atom_id", Value::as_str)?.to_string();

                struct_conns.push(StructConn {
                    ptnr1_asym_id,
                    ptnr1_seq_id,
                    ptnr1_atom_id,
                    ptnr2_asym_id,
                    ptnr2_seq_id,
                    ptnr2_atom_id,
                });
            }
        }

        Ok(Self {
            entities,
            struct_asyms,
            atom_sites,
            chem_comp_bonds,
            struct_conns,
        })
    }

    // TODO: return Result<Vec<Molecule>, MmCifError> instead of unwrapping —
    // a struct_asym with no matching entity_id or no atom_site rows is a
    // malformed-file problem, not an invariant of this program, same as the
    // field-level errors in `new`.
    pub fn build_molecules(&self) -> Vec<Molecule> {
        let mut molecules: Vec<Molecule> = vec![];
        for struct_asym in self.struct_asym() {
            let entity = self.entity(struct_asym.entity_id).unwrap();
            let mut residues: Vec<Residue> = vec![];
            let mut atoms: Vec<Atom> = vec![];
            let mut positions: Vec<Position> = vec![];
            let mut occupancies: Vec<f32> = vec![];
            let mut b_factors: Vec<f32> = vec![];
            let atom_sites = self.atom_sites(&struct_asym.id).unwrap();
            let mut residue_start = 0;
            for (i, atom_site) in atom_sites.iter().enumerate() {
                atoms.push(Atom::new(atom_site.type_symbol, atom_site.label_atom_id.clone()));
                positions.push(Position::new(atom_site.x, atom_site.y, atom_site.z));
                occupancies.push(atom_site.occupancy);
                b_factors.push(atom_site.b_factor);

                let residue_end = i + 1 == atom_sites.len() || atom_sites[i + 1].auth_seq_id != atom_site.auth_seq_id;
                if residue_end {
                    residues.push(Residue::new(
                        atom_site.label_comp_id.clone(),
                        Range {
                            start: residue_start,
                            end: i + 1, // end is exclusive
                        },
                        atom_sites[i].is_hetero,
                    ));
                    residue_start = i + 1;
                }
            }
            molecules.push(Molecule::new(
                Topology::new(entity.name.clone(), atoms, residues),
                Conformer::new(positions, occupancies, b_factors),
            ));
        }
        molecules
    }

    pub fn entities(&self) -> impl Iterator<Item = &Entity> {
        self.entities.values()
    }

    pub fn entity(&self, id: u8) -> Option<&Entity> {
        self.entities.get(&id)
    }

    pub fn atom_sites(&self, struct_asym_id: &str) -> Option<&Vec<AtomSite>> {
        self.atom_sites.get(struct_asym_id)
    }

    pub fn struct_asym(&self) -> impl Iterator<Item = &StructAsym> {
        self.struct_asyms.values()
    }
}
