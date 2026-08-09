use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::Range;

use crate::atom::{Atom, Element};
use crate::bond::Bond;
use crate::bond_order::{BondOrder, BondOrderSource};
use crate::component::Component;
use crate::frame::Frame;
use crate::geometry::Point3;
use crate::molecule::Molecule;
use crate::protein::Protein;
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

fn bond_order(v: &Value) -> Option<BondOrder> {
    match v.as_str()? {
        "sing" => Some(BondOrder::Single),
        "doub" => Some(BondOrder::Double),
        "trip" => Some(BondOrder::Triple),
        "arom" => Some(BondOrder::Aromatic),
        "pi" => Some(BondOrder::Pi),
        _ => None,
    }
}

fn flag(v: &Value) -> Option<bool> {
    match v.as_str()? {
        "N" => Some(false),
        "Y" => Some(true),
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

#[derive(Debug)]
pub struct AtomSite {
    pub id: u32,
    pub label_entity_id: u8,
    pub label_asym_id: String,
    pub label_atom_id: String,
    pub label_alt_id: Option<String>,
    // We use auth_seq_id because for non-polymers the label is `.`, which will end up merging all the atoms into one
    // residue
    pub auth_seq_id: i64,
    // Disambiguates residues that share an auth_seq_id (e.g. an inserted residue kept under a reference numbering scheme).
    // For example: 3PTB chain A residue 184 (TYR) vs 184A (GLY).
    pub pdbx_ins_code: Option<String>,
    pub label_comp_id: String,
    pub type_symbol: Element,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f32,
    pub b_factor: f32,
}

// (alt_id -> (positions, occupancies, b_factors)) accumulated per altloc while walking one
// struct_asym's atom_site rows; collapsed into `Vec<Frame>` once the walk is done.
type FrameBuildData = HashMap<Option<String>, (Vec<Point3>, Vec<f32>, Vec<f32>)>;

#[derive(Debug)]
pub struct ChemCompBond {
    pub value_order: BondOrder,
    pub aromatic_flag: bool,
}

#[derive(Debug)]
pub struct MmCIF {
    pub entities: HashMap<u8, Entity>,
    pub struct_asyms: HashMap<String, StructAsym>,
    pub atom_sites: HashMap<String, Vec<AtomSite>>,
    pub chem_comp_bonds: HashMap<(String, String, String), ChemCompBond>,
}

impl BondOrderSource for MmCIF {
    fn bond_order(&self, bond: &Bond) -> Option<BondOrder> {
        let topology = bond.atom_1.conformer.topology;
        let component_1 = topology.component_index_for(bond.atom_1.index);
        let component_2 = topology.component_index_for(bond.atom_2.index);
        if component_1 != component_2 {
            return None;
        }
        let component = &topology.components[component_1];
        let atom_1 = &topology.atoms[bond.atom_1.index];
        let atom_2 = &topology.atoms[bond.atom_2.index];
        // TODO: this ignores `chem_comp_bond.aromatic_flag` — a ring bond currently resolves to
        // whatever Kekulé value_order (Single/Double) the dictionary happens to record, never
        // `BondOrder::Aromatic`. If `aromatic_flag` is true, use `BondOrder::Aromatic` instead of
        // `value_order`.
        self.chem_comp_bonds_lookup(&component.name, &atom_1.name, &atom_2.name)
            .map(|chem_comp_bond| chem_comp_bond.value_order)
    }
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
                let label_alt_id = optional_field(&row, "_atom_site.label_alt_id", Value::as_str)?.map(str::to_string);
                let auth_seq_id = field(&row, "_atom_site.auth_seq_id", Value::as_int)?;
                let pdbx_ins_code =
                    optional_field(&row, "_atom_site.pdbx_PDB_ins_code", Value::as_str)?.map(str::to_string);
                let label_comp_id = field(&row, "_atom_site.label_comp_id", Value::as_str)?.to_string();
                let type_symbol = field(&row, "_atom_site.type_symbol", element)?;
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
                    label_alt_id,
                    auth_seq_id,
                    pdbx_ins_code,
                    label_comp_id,
                    type_symbol,
                    x,
                    y,
                    z,
                    occupancy,
                    b_factor,
                });
            }
        }

        let mut chem_comp_bonds: HashMap<(String, String, String), ChemCompBond> = HashMap::new();
        if let Some(loop_) = block.loop_by_category("_chem_comp_bond") {
            for row in loop_.iter() {
                let comp_id = field(&row, "_chem_comp_bond.comp_id", Value::as_str)?.to_string();
                let atom_id_1 = field(&row, "_chem_comp_bond.atom_id_1", Value::as_str)?.to_string();
                let atom_id_2 = field(&row, "_chem_comp_bond.atom_id_2", Value::as_str)?.to_string();
                let value_order = field(&row, "_chem_comp_bond.value_order", bond_order)?;
                let aromatic_flag = field(&row, "_chem_comp_bond.pdbx_aromatic_flag", flag)?;
                chem_comp_bonds.insert(
                    (comp_id, atom_id_1, atom_id_2),
                    ChemCompBond {
                        value_order,
                        aromatic_flag,
                    },
                );
            }
        }

        Ok(Self {
            entities,
            struct_asyms,
            atom_sites,
            chem_comp_bonds,
        })
    }

    // TODO: return Result<Vec<Protein>, MmCifError> / Result<Vec<Molecule>, MmCifError> instead of
    // unwrapping — a struct_asym with no matching entity_id or no atom_site rows is a
    // malformed-file problem, not an invariant of this program, same as the field-level errors in
    // `new`.
    pub fn build_proteins(&self) -> Vec<Protein> {
        let mut proteins = vec![];
        for struct_asym in self.struct_asym() {
            let entity = self.entity(struct_asym.entity_id).unwrap();
            if !matches!(entity.kind, EntityKind::Polymer) {
                continue;
            }
            let (name, topology, frames, residue_numbers) =
                self.build_topology_frames_components(struct_asym, entity, proteins.len());
            proteins.push(Protein::new(name, topology, frames, residue_numbers));
        }
        proteins
    }

    pub fn build_molecules(&self) -> Vec<Molecule> {
        let mut molecules = vec![];
        for struct_asym in self.struct_asym() {
            let entity = self.entity(struct_asym.entity_id).unwrap();
            if matches!(entity.kind, EntityKind::Polymer) {
                continue;
            }
            let (name, topology, frames, _residue_numbers) =
                self.build_topology_frames_components(struct_asym, entity, molecules.len());
            molecules.push(Molecule::new(molecules.len(), name, topology, frames));
        }
        molecules
    }

    /// Shared by `build_proteins`/`build_molecules`: turns one `struct_asym`'s atom_site rows into
    /// a `Topology` (atoms + `Component`s — one `Component` per distinct `(auth_seq_id,
    /// pdbx_ins_code)`, always at least one, even for a single-component ligand), one `Frame` per
    /// distinct `label_alt_id` (or a single frame if there are no altlocs), and the residue/
    /// nucleotide numbers parallel to `Topology.components` (unused by `build_molecules` today,
    /// but cheap to compute once here rather than duplicating this whole loop per caller).
    fn build_topology_frames_components(
        &self,
        struct_asym: &StructAsym,
        entity: &Entity,
        topology_id: usize,
    ) -> (String, Topology, Vec<Frame>, Vec<i64>) {
        let mut frame_data: FrameBuildData = HashMap::new();
        frame_data.insert(None, (vec![], vec![], vec![]));
        let mut components: Vec<Component> = vec![];
        let mut numbers: Vec<i64> = vec![];
        let mut atoms: Vec<Atom> = vec![];
        let atom_sites = self.atom_sites(&struct_asym.id).unwrap();
        let mut component_start = 0;
        // (auth_seq_id, pdbx_ins_code, label_atom_id) is the real atom identity. auth_seq_id
        // alone collides for inserted residues (e.g. 3PTB chain A: TYR 184 and GLY 184A both
        // have an atom named "N"), and label_alt_id must be ignored here since alt records
        // are the *same* atom, not different ones.
        let mut seen_atoms: HashSet<(i64, Option<String>, String)> = HashSet::new();
        for (i, atom_site) in atom_sites.iter().enumerate() {
            let identity = (
                atom_site.auth_seq_id,
                atom_site.pdbx_ins_code.clone(),
                atom_site.label_atom_id.clone(),
            );
            if seen_atoms.insert(identity) {
                atoms.push(Atom::new(atom_site.type_symbol, atom_site.label_atom_id.clone()));
            }
            match &atom_site.label_alt_id {
                None => {
                    for (_, frame) in frame_data.iter_mut() {
                        frame.0.push(Point3::new(atom_site.x, atom_site.y, atom_site.z));
                        frame.1.push(atom_site.occupancy);
                        frame.2.push(atom_site.b_factor);
                    }
                }
                alt_id => match frame_data.get_mut(alt_id) {
                    None => {
                        let template = frame_data.get(&None).unwrap();
                        let mut new_positions = template.0.clone();
                        new_positions.push(Point3::new(atom_site.x, atom_site.y, atom_site.z));
                        let mut new_occupancies = template.1.clone();
                        new_occupancies.push(atom_site.occupancy);
                        let mut new_b_factors = template.2.clone();
                        new_b_factors.push(atom_site.b_factor);
                        frame_data.insert(alt_id.clone(), (new_positions, new_occupancies, new_b_factors));
                    }
                    Some(frame) => {
                        frame.0.push(Point3::new(atom_site.x, atom_site.y, atom_site.z));
                        frame.1.push(atom_site.occupancy);
                        frame.2.push(atom_site.b_factor);
                    }
                },
            }

            let build_component = i + 1 == atom_sites.len()
                || atom_sites[i + 1].auth_seq_id != atom_site.auth_seq_id
                || atom_sites[i + 1].pdbx_ins_code != atom_site.pdbx_ins_code;
            if build_component {
                components.push(Component::new(
                    atom_site.label_comp_id.clone(),
                    Range {
                        start: component_start,
                        end: atoms.len(), // end is exclusive
                    },
                ));
                numbers.push(atom_site.auth_seq_id);
                component_start = atoms.len();
            }
        }
        let keep_none = frame_data.len() == 1;
        let frames: Vec<Frame> = frame_data
            .into_iter()
            .filter(|(alt_id, _)| alt_id.is_some() || keep_none)
            .enumerate()
            .map(|(frame_id, (_, (positions, occupancies, b_factors)))| {
                assert!(positions.len() == atoms.len());
                Frame::new(frame_id, positions, occupancies, b_factors)
            })
            .collect();
        let topology = Topology::new(topology_id, atoms, components);
        (entity.name.clone(), topology, frames, numbers)
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

    fn chem_comp_bonds_lookup(&self, comp_id: &str, atom_1: &str, atom_2: &str) -> Option<&ChemCompBond> {
        self.chem_comp_bonds
            .get(&(comp_id.to_string(), atom_1.to_string(), atom_2.to_string()))
            .or_else(|| {
                self.chem_comp_bonds
                    .get(&(comp_id.to_string(), atom_2.to_string(), atom_1.to_string()))
            })
    }
}
