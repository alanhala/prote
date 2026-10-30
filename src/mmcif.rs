use std::collections::{BTreeSet, HashMap};
use std::hash::Hash;
use std::vec;

use crate::bond::{Bond, ResidueTemplateBond};
use crate::bond_order::BondOrder;
use crate::element::Element;
use crate::frame::Frame;
use crate::geometry::Point3;
use crate::structure::{Atom, Entity, EntityInstance, EntityKind, PolymerKind, Residue, Structure};
use crate::{Cif, Loop, Row, Value};

// Entries in first-insertion order. Fields are private because `indices` must keep pointing into
// `entries`.
#[derive(Debug)]
pub struct OrderedMap<K, V> {
    indices: HashMap<K, usize>,
    entries: Vec<(K, V)>,
}

impl<K: Eq + Hash + Clone, V> OrderedMap<K, V> {
    pub fn new() -> Self {
        Self { indices: HashMap::new(), entries: vec![] }
    }

    pub fn get_or_insert_with(&mut self, key: K, value: impl FnOnce() -> V) -> &mut V {
        let index = *self.indices.entry(key).or_insert_with_key(|key| {
            self.entries.push((key.clone(), value()));
            self.entries.len() - 1
        });
        &mut self.entries[index].1
    }

    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.entries.iter().map(|(key, value)| (key, value))
    }

    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.entries.iter().map(|(_, value)| value)
    }

    pub fn index_of(&self, key: &K) -> Option<usize> {
        self.indices.get(key).copied()
    }
}

#[derive(Debug)]
pub enum MmCifError {
    NoDataBlock,
    MissingField(&'static str),
    WrongType(&'static str),
    UnknownEntity(u8),
    MissingAtom { atom: usize, model: i64 },
}

// A tag's position in its category, looked up once before reading the rows: searching every row's
// tags by name was most of the time spent reading `_atom_site`. A missing tag stays an error per
// row, so an absent category (no tags, no rows) still reads as empty.
#[derive(Clone, Copy)]
struct Column {
    tag: &'static str,
    index: Option<usize>,
}

impl Column {
    fn new(category: &Loop, tag: &'static str) -> Self {
        Self { tag, index: category.tag_index(tag) }
    }

    fn get<'a>(&self, row: &Row<'a>) -> Option<&'a Value> {
        row.at(self.index?)
    }
}

fn field<'a, T>(row: &Row<'a>, column: Column, cast: fn(&'a Value) -> Option<T>) -> Result<T, MmCifError> {
    match column.get(row) {
        None => Err(MmCifError::MissingField(column.tag)),
        Some(v) => cast(v).ok_or(MmCifError::WrongType(column.tag)),
    }
}

/// Like `field`, but a CIF `.`/`?` (inapplicable/unknown) is `Ok(None)` rather than
/// an error — e.g. `_struct_conn` partner seq_id, which is inapplicable for
/// non-polymer partners (a bound ion has no sequence position).
fn optional_field<'a, T>(row: &Row<'a>, column: Column, cast: fn(&'a Value) -> Option<T>) -> Result<Option<T>, MmCifError> {
    match column.get(row) {
        None => Err(MmCifError::MissingField(column.tag)),
        Some(Value::Inapplicable) | Some(Value::Unknown) => Ok(None),
        Some(v) => cast(v).map(Some).ok_or(MmCifError::WrongType(column.tag)),
    }
}

fn string(v: &Value) -> Option<String> {
    v.as_str().map(str::to_string)
}

// Identifiers like `auth_asym_id` are often digits ("1"), which the lexer reads as integers rather
// than strings.
fn code(v: &Value) -> Option<String> {
    string(v).or_else(|| v.as_int().map(|n| n.to_string()))
}

// A polymer's kind lives in a separate category (`_entity_poly`), so it's passed in.
fn entity_kind(v: &Value, polymer: Option<PolymerKind>) -> Option<EntityKind> {
    match v.as_str()? {
        "polymer" => Some(EntityKind::Polymer(polymer?)),
        "non-polymer" => Some(EntityKind::NonPolymer),
        "branched" => Some(EntityKind::Branched),
        "water" => Some(EntityKind::Water),
        _ => None,
    }
}

fn polymer_kind(v: &Value) -> Option<PolymerKind> {
    match v.as_str()? {
        "polypeptide(L)" => Some(PolymerKind::PolypeptideL),
        "polypeptide(D)" => Some(PolymerKind::PolypeptideD),
        "cyclic-pseudo-peptide" => Some(PolymerKind::CyclicPseudoPeptide),
        "polydeoxyribonucleotide" => Some(PolymerKind::Dna),
        "polyribonucleotide" => Some(PolymerKind::Rna),
        "polydeoxyribonucleotide/polyribonucleotide hybrid" => Some(PolymerKind::DnaRnaHybrid),
        "peptide nucleic acid" => Some(PolymerKind::PeptideNucleicAcid),
        "other" => Some(PolymerKind::Other),
        _ => None,
    }
}

// RCSB writes these in lowercase ("doub"), PDBe in uppercase ("DOUB").
fn bond_order(v: &Value) -> Option<BondOrder> {
    match v.as_str()?.to_ascii_lowercase().as_str() {
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

fn element(v: &Value) -> Option<Element> {
    let symbol = v.as_str()?;
    let mut chars = symbol.chars();
    let normalized = match chars.next() {
        Some(first) => first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase(),
        None => return None,
    };
    Element::from_symbol(&normalized)
}

#[derive(Debug)]
pub struct ChemCompBond {
    pub atom_1: String,
    pub atom_2: String,
    pub value_order: BondOrder,
    pub aromatic_flag: bool,
}

#[derive(Debug)]
pub struct AtomSite {
    pub type_symbol: Element,
    pub label_alt_id: Option<String>,
    pub position: Point3,
    pub occupancy: f32,
    pub b_factor: f32,
    pub model: i64,
}

#[derive(Debug)]
pub struct InstanceRows {
    pub entity_id: u8,
    pub chain: String,
    // Keyed by (auth_seq_id, pdbx_PDB_ins_code). We use auth_seq_id because for non-polymers the
    // label is `.`, which would merge all their atoms into one residue; the insertion code
    // separates residues that share it (e.g. 3PTB chain A: TYR 184 and GLY 184A).
    pub residues: OrderedMap<(i64, Option<String>), ResidueRows>,
}

#[derive(Debug)]
pub struct ResidueRows {
    pub name: String,
    pub seq_id: Option<i64>,
    pub atoms: OrderedMap<String, Vec<AtomSite>>,
}

#[derive(Debug)]
pub struct MmCIF {
    pub entities: OrderedMap<u8, Entity>,
    pub instances: OrderedMap<String, InstanceRows>,
    pub chem_comp_bonds: HashMap<String, Vec<ChemCompBond>>,
}

impl MmCIF {
    pub fn new(cif: Cif) -> Result<Self, MmCifError> {
        let block = cif.blocks.into_iter().next().ok_or(MmCifError::NoDataBlock)?;

        let mut polymer_kinds = HashMap::new();
        let rows = block.category("_entity_poly");
        let [entity_id, kind] = ["_entity_poly.entity_id", "_entity_poly.type"].map(|tag| Column::new(&rows, tag));
        for row in rows.iter() {
            let id = field(&row, entity_id, Value::as_int)? as u8;
            polymer_kinds.insert(id, field(&row, kind, polymer_kind)?);
        }

        let mut entities = OrderedMap::new();
        let rows = block.category("_entity");
        let [id, kind, description] = ["_entity.id", "_entity.type", "_entity.pdbx_description"].map(|tag| Column::new(&rows, tag));
        for row in rows.iter() {
            let id = field(&row, id, Value::as_int)? as u8;
            let kind = match kind.get(&row) {
                None => return Err(MmCifError::MissingField(kind.tag)),
                // A polymer without an `_entity_poly` row has no kind, so it fails here too.
                Some(v) => entity_kind(v, polymer_kinds.get(&id).copied()).ok_or(MmCifError::WrongType(kind.tag))?,
            };
            let entity = Entity { kind, name: field(&row, description, string)? };
            entities.get_or_insert_with(id, || entity);
        }

        // Each value below shadows the column it is read from.
        let mut instances = OrderedMap::new();
        let rows = block.category("_atom_site");
        let [label_asym_id, label_entity_id, auth_asym_id, auth_seq_id, pdbx_ins_code, label_comp_id, label_atom_id, label_seq_id, type_symbol, label_alt_id, cartn_x, cartn_y, cartn_z, occupancy, b_factor, model] =
            [
                "_atom_site.label_asym_id",
                "_atom_site.label_entity_id",
                "_atom_site.auth_asym_id",
                "_atom_site.auth_seq_id",
                "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.label_comp_id",
                "_atom_site.label_atom_id",
                "_atom_site.label_seq_id",
                "_atom_site.type_symbol",
                "_atom_site.label_alt_id",
                "_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z",
                "_atom_site.occupancy",
                "_atom_site.B_iso_or_equiv",
                "_atom_site.pdbx_PDB_model_num",
            ]
            .map(|tag| Column::new(&rows, tag));
        for row in rows.iter() {
            let label_asym_id = field(&row, label_asym_id, string)?;
            let label_entity_id = field(&row, label_entity_id, Value::as_int)? as u8;
            let auth_asym_id = field(&row, auth_asym_id, code)?;
            let auth_seq_id = field(&row, auth_seq_id, Value::as_int)?;
            let pdbx_ins_code = optional_field(&row, pdbx_ins_code, code)?;
            let label_comp_id = field(&row, label_comp_id, string)?;
            let label_atom_id = field(&row, label_atom_id, string)?;
            let label_seq_id = optional_field(&row, label_seq_id, Value::as_int)?;
            let atom_site = AtomSite {
                type_symbol: field(&row, type_symbol, element)?,
                label_alt_id: optional_field(&row, label_alt_id, string)?,
                position: Point3::new(
                    field(&row, cartn_x, Value::as_float)?,
                    field(&row, cartn_y, Value::as_float)?,
                    field(&row, cartn_z, Value::as_float)?,
                ),
                occupancy: field(&row, occupancy, Value::as_float)? as f32,
                b_factor: field(&row, b_factor, Value::as_float)? as f32,
                model: field(&row, model, Value::as_int)?,
            };
            let instance = instances.get_or_insert_with(label_asym_id, || InstanceRows {
                entity_id: label_entity_id,
                chain: auth_asym_id,
                residues: OrderedMap::new(),
            });
            let residue = instance.residues.get_or_insert_with((auth_seq_id, pdbx_ins_code), || ResidueRows {
                name: label_comp_id.clone(),
                seq_id: label_seq_id,
                atoms: OrderedMap::new(),
            });
            // Microheterogeneity: altlocs can model two residue types at one position (e.g.
            // PRO/SER). Their atoms differ, not just their positions, so frames can't represent
            // both; keep the first.
            if label_comp_id != residue.name {
                continue;
            }
            residue.atoms.get_or_insert_with(label_atom_id, Vec::new).push(atom_site);
        }

        let mut chem_comp_bonds: HashMap<String, Vec<ChemCompBond>> = HashMap::new();
        let rows = block.category("_chem_comp_bond");
        let [comp_id, atom_id_1, atom_id_2, value_order, aromatic_flag] = [
            "_chem_comp_bond.comp_id",
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
            "_chem_comp_bond.value_order",
            "_chem_comp_bond.pdbx_aromatic_flag",
        ]
        .map(|tag| Column::new(&rows, tag));
        for row in rows.iter() {
            chem_comp_bonds.entry(field(&row, comp_id, string)?).or_default().push(ChemCompBond {
                atom_1: field(&row, atom_id_1, string)?,
                atom_2: field(&row, atom_id_2, string)?,
                value_order: field(&row, value_order, bond_order)?,
                aromatic_flag: field(&row, aromatic_flag, flag)?,
            });
        }

        Ok(Self { entities, instances, chem_comp_bonds })
    }

    pub fn build_structure(&self) -> Result<(Structure, Vec<Frame>), MmCifError> {
        let mut atoms: Vec<Atom> = vec![];
        let mut residues: Vec<Residue> = vec![];
        let mut instances: Vec<EntityInstance> = vec![];
        let mut atom_rows: Vec<&[AtomSite]> = vec![];
        for (id, instance) in self.instances.iter() {
            let entity = self.entities.index_of(&instance.entity_id).ok_or(MmCifError::UnknownEntity(instance.entity_id))?;
            let residue_start = residues.len();
            for ((number, _), residue) in instance.residues.iter() {
                let atom_start = atoms.len();
                for (name, rows) in residue.atoms.iter() {
                    atoms.push(Atom { element: rows[0].type_symbol, name: name.clone() });
                    atom_rows.push(rows);
                }
                residues.push(Residue { name: residue.name.clone(), number: *number, seq_id: residue.seq_id, atoms: atom_start..atoms.len() });
            }
            instances.push(EntityInstance { id: id.clone(), chain: instance.chain.clone(), entity, residues: residue_start..residues.len() });
        }
        let frames = build_frames(&atom_rows)?;
        let templates = build_template_bonds(&self.chem_comp_bonds);
        let mut structure = Structure { atoms, residues, instances, entities: self.entities.values().cloned().collect(), bonds: vec![] };
        structure.bonds = Bond::bonds(&structure, &templates);
        Ok((structure, frames))
    }
}

fn build_template_bonds(chem_comp_bonds: &HashMap<String, Vec<ChemCompBond>>) -> HashMap<String, Vec<ResidueTemplateBond>> {
    let mut templates: HashMap<String, Vec<ResidueTemplateBond>> = HashMap::new();
    for (residue_name, bonds) in chem_comp_bonds {
        for bond in bonds {
            templates.entry(residue_name.to_string()).or_default().push(ResidueTemplateBond {
                atom_name_1: bond.atom_1.to_string(),
                atom_name_2: bond.atom_2.to_string(),
                aromatic: bond.aromatic_flag,
                order: bond.value_order,
            });
        }
    }
    templates
}

// One frame per (model, altloc). Altlocs are local to a few residues, so an atom without a row for
// the frame's altloc takes its row without an altloc, else its highest-occupancy row; that keeps
// every frame covering every atom.
fn build_frames(atom_rows: &[&[AtomSite]]) -> Result<Vec<Frame>, MmCifError> {
    let all_rows = || atom_rows.iter().flat_map(|rows| rows.iter());
    let models: BTreeSet<i64> = all_rows().map(|atom_site| atom_site.model).collect();
    let alt_ids: BTreeSet<&str> = all_rows().filter_map(|atom_site| atom_site.label_alt_id.as_deref()).collect();
    let alt_ids: Vec<Option<&str>> = if alt_ids.is_empty() { vec![None] } else { alt_ids.into_iter().map(Some).collect() };
    let mut frames = vec![];
    for &model in &models {
        for &alt_id in &alt_ids {
            let picked = atom_rows
                .iter()
                .enumerate()
                .map(|(atom, rows)| pick_row(rows, model, alt_id).ok_or(MmCifError::MissingAtom { atom, model }))
                .collect::<Result<Vec<_>, _>>()?;
            frames.push(Frame::new(
                frames.len(),
                picked.iter().map(|atom_site| atom_site.position).collect(),
                picked.iter().map(|atom_site| atom_site.occupancy).collect(),
                picked.iter().map(|atom_site| atom_site.b_factor).collect(),
            ));
        }
    }
    Ok(frames)
}

fn pick_row<'a>(rows: &'a [AtomSite], model: i64, alt_id: Option<&str>) -> Option<&'a AtomSite> {
    let in_model = || rows.iter().filter(move |atom_site| atom_site.model == model);
    in_model()
        .find(|atom_site| atom_site.label_alt_id.as_deref() == alt_id)
        .or_else(|| in_model().find(|atom_site| atom_site.label_alt_id.is_none()))
        .or_else(|| in_model().reduce(|best, atom_site| if atom_site.occupancy > best.occupancy { atom_site } else { best }))
}
