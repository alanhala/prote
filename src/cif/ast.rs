#[derive(Debug)]
pub struct Cif {
    pub(crate) blocks: Vec<DataBlock>,
}

#[derive(Debug)]
pub struct DataBlock {
    pub name: String,
    pub items: Vec<Member>,
}

#[derive(Debug)]
pub enum Member {
    Item { tag: String, value: Value },
    Loop(Loop),
    SaveFrame(SaveFrame),
}

#[derive(Debug)]
pub enum Value {
    Integer { value: i64, su: Option<u64> },
    Float { value: f64, su: Option<u64> },
    Str(String),
    Inapplicable,
    Unknown,
}

#[derive(Debug)]
pub struct Loop {
    pub(crate) tags: Vec<String>,
    pub(crate) rows: Vec<Vec<Value>>,
}

#[derive(Debug)]
pub struct SaveFrame {
    pub(crate) name: String,
    pub(crate) items: Vec<Member>,
}

#[derive(Debug)]
pub struct Row<'a> {
    tags: &'a [String],
    values: &'a [Value],
}

impl<'a> Row<'a> {
    pub fn get(&self, tag: &str) -> Option<&'a Value> {
        let i = self.tags.iter().position(|t| t.eq_ignore_ascii_case(tag))?;
        self.values.get(i)
    }
}

impl Cif {
    pub fn block(&self, name: &str) -> Option<&DataBlock> {
        self.blocks.iter().find(|b| b.name.eq_ignore_ascii_case(name))
    }
}

impl DataBlock {
    pub fn loop_by_category(&self, category: &str) -> Option<&Loop> {
        self.items.iter().find_map(|m| match m {
            Member::Loop(l) if l.category().is_some_and(|c| c.eq_ignore_ascii_case(category)) => Some(l),
            _ => None,
        })
    }

    pub fn get(&self, tag: &str) -> Option<&Value> {
        self.items.iter().find_map(|m| match m {
            Member::Item { tag: t, value } if t.eq_ignore_ascii_case(tag) => Some(value),
            _ => None,
        })
    }

    pub fn find_loop(&self, tag: &str) -> Option<&Loop> {
        self.items.iter().find_map(|m| match m {
            Member::Loop(l) if l.has_tag(tag) => Some(l),
            _ => None,
        })
    }

    pub fn column(&self, tag: &str) -> Option<impl Iterator<Item = &Value> + '_> {
        self.find_loop(tag)?.column(tag)
    }
}

impl Loop {
    fn tag_index(&self, tag: &str) -> Option<usize> {
        self.tags.iter().position(|t| t.eq_ignore_ascii_case(tag))
    }

    pub fn category(&self) -> Option<&str> {
        self.tags.first()?.split('.').next()
    }

    pub fn has_tag(&self, tag: &str) -> bool {
        self.tag_index(tag).is_some()
    }

    pub fn column(&self, tag: &str) -> Option<impl Iterator<Item = &Value> + '_> {
        let i = self.tag_index(tag)?;
        Some(self.rows.iter().map(move |row| &row[i]))
    }

    pub fn iter(&self) -> impl Iterator<Item = Row<'_>> + '_ {
        let tags: &[String] = &self.tags;
        self.rows.iter().map(move |values| Row { tags, values })
    }
}

impl Value {
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }
    pub fn as_int(&self) -> Option<i64> {
        match self {
            Value::Integer { value, .. } => Some(*value),
            _ => None,
        }
    }
    pub fn as_float(&self) -> Option<f64> {
        match self {
            Value::Float { value, .. } => Some(*value),
            Value::Integer { value, .. } => Some(*value as f64),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // 3.1 Lookups

    #[test]
    fn block_lookup_is_case_insensitive() {
        let cif = Cif::new("data_4d1m\n_a 1");
        assert!(cif.block("4D1M").is_some());
    }

    #[test]
    fn block_lookup_returns_none_for_missing_name() {
        let cif = Cif::new("data_4d1m\n_a 1");
        assert!(cif.block("nope").is_none());
    }

    #[test]
    fn item_get_is_case_insensitive() {
        let cif = Cif::new("data_x\n_entry.id 1ABC");
        let block = cif.block("x").unwrap();
        assert!(matches!(block.get("_Entry.Id"), Some(Value::Str(s)) if s == "1ABC"));
    }

    #[test]
    fn item_get_returns_first_match_when_tag_repeats() {
        let cif = Cif::new("data_x\n_a 1\n_a 2");
        let block = cif.block("x").unwrap();
        assert!(matches!(block.get("_a"), Some(Value::Integer { value: 1, .. })));
    }

    #[test]
    fn item_get_returns_none_for_a_loop_only_tag() {
        let cif = Cif::new("data_x\nloop_\n_a\n1\n2");
        let block = cif.block("x").unwrap();
        assert!(block.get("_a").is_none());
    }

    // 3.2 Loop access

    fn atom_site_block() -> Cif {
        Cif::new(
            "data_x\nloop_\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n1.0 2.0\n3.0 4.0",
        )
    }

    #[test]
    fn find_loop_is_case_insensitive() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        assert!(block.find_loop("_ATOM_SITE.CARTN_X").is_some());
    }

    #[test]
    fn column_yields_values_in_row_order() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        let l = block.find_loop("_atom_site.Cartn_x").unwrap();
        let ys: Vec<f64> = l.column("_atom_site.Cartn_y").unwrap().map(|v| v.as_float().unwrap()).collect();
        assert_eq!(ys, vec![2.0, 4.0]);
    }

    #[test]
    fn column_returns_none_for_tag_not_in_any_loop() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        let l = block.find_loop("_atom_site.Cartn_x").unwrap();
        assert!(l.column("_not.a_tag").is_none());
    }

    #[test]
    fn column_reads_the_correct_index_on_a_multi_column_loop() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        let l = block.find_loop("_atom_site.Cartn_x").unwrap();
        let xs: Vec<f64> = l.column("_atom_site.Cartn_x").unwrap().map(|v| v.as_float().unwrap()).collect();
        assert_eq!(xs, vec![1.0, 3.0]);
    }

    #[test]
    fn category_splits_the_first_tag_on_the_dot() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        let l = block.find_loop("_atom_site.Cartn_x").unwrap();
        assert_eq!(l.category(), Some("_atom_site"));
    }

    #[test]
    fn category_on_a_dotless_tag_returns_the_whole_tag() {
        let cif = Cif::new("data_x\nloop_\n_a\n1\n2");
        let block = cif.block("x").unwrap();
        let l = block.find_loop("_a").unwrap();
        assert_eq!(l.category(), Some("_a"));
    }

    #[test]
    fn loop_by_category_is_case_insensitive_and_none_when_missing() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        assert!(block.loop_by_category("_ATOM_SITE").is_some());
        assert!(block.loop_by_category("_not_a_category").is_none());
    }

    #[test]
    fn loop_iter_yields_rows_in_order_and_row_get_is_case_insensitive() {
        let cif = Cif::new("data_x\nloop_\n_a\n_b\n1 2\n3 4\n5 6");
        let block = cif.block("x").unwrap();
        let l = block.find_loop("_a").unwrap();
        let rows: Vec<Row> = l.iter().collect();

        assert_eq!(rows.len(), 3);
        assert!(matches!(rows[0].get("_a"), Some(Value::Integer { value: 1, .. })));
        assert!(matches!(rows[1].get("_b"), Some(Value::Integer { value: 4, .. })));
        assert!(matches!(rows[2].get("_A"), Some(Value::Integer { value: 5, .. })));
        assert!(rows[0].get("_missing").is_none());
    }

    #[test]
    fn loop_with_quoted_string_values_via_column_and_row_get() {
        let cif = Cif::new("data_x\nloop_\n_label\n_note\n'foo' 'bar'\n'baz' 'qux'");
        let block = cif.block("x").unwrap();

        let labels: Vec<&str> = block.column("_label").unwrap().map(|v| v.as_str().unwrap()).collect();
        assert_eq!(labels, vec!["foo", "baz"]);

        let l = block.find_loop("_label").unwrap();
        let second_row = l.iter().nth(1).unwrap();
        assert_eq!(second_row.get("_note").unwrap().as_str(), Some("qux"));
    }

    #[test]
    fn data_block_column_delegates_to_the_loops_column() {
        let cif = atom_site_block();
        let block = cif.block("x").unwrap();
        let via_block: Vec<f64> = block.column("_atom_site.Cartn_y").unwrap().map(|v| v.as_float().unwrap()).collect();
        let via_loop: Vec<f64> = block
            .find_loop("_atom_site.Cartn_y")
            .unwrap()
            .column("_atom_site.Cartn_y")
            .unwrap()
            .map(|v| v.as_float().unwrap())
            .collect();
        assert_eq!(via_block, via_loop);
    }

    // 3.3 Value accessors

    #[test]
    fn integer_as_int() {
        assert_eq!(Value::Integer { value: 5, su: None }.as_int(), Some(5));
    }

    #[test]
    fn integer_as_float_is_promoted() {
        assert_eq!(Value::Integer { value: 5, su: None }.as_float(), Some(5.0));
    }

    #[test]
    fn float_as_float() {
        assert_eq!(Value::Float { value: 5.5, su: None }.as_float(), Some(5.5));
    }

    #[test]
    fn float_as_int_is_none() {
        assert_eq!(Value::Float { value: 5.5, su: None }.as_int(), None);
    }

    #[test]
    fn as_str_on_non_string_variants_is_none() {
        assert_eq!(Value::Integer { value: 1, su: None }.as_str(), None);
        assert_eq!(Value::Float { value: 1.0, su: None }.as_str(), None);
        assert_eq!(Value::Inapplicable.as_str(), None);
        assert_eq!(Value::Unknown.as_str(), None);
    }

    #[test]
    fn accessors_on_inapplicable_and_unknown_are_none() {
        for v in [Value::Inapplicable, Value::Unknown] {
            assert_eq!(v.as_int(), None);
            assert_eq!(v.as_float(), None);
            assert_eq!(v.as_str(), None);
        }
    }

    #[test]
    fn standard_uncertainty_survives_on_the_value_but_accessors_drop_it() {
        let v = Value::Integer { value: 5, su: Some(2) };
        assert_eq!(v.as_int(), Some(5));
        assert!(matches!(v, Value::Integer { su: Some(2), .. }));
    }
}
