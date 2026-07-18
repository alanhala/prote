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

impl Cif {
    pub fn block(&self, name: &str) -> Option<&DataBlock> {
        self.blocks.iter().find(|b| b.name.eq_ignore_ascii_case(name))
    }
}

impl DataBlock {
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

    pub fn has_tag(&self, tag: &str) -> bool {
        self.tag_index(tag).is_some()
    }

    pub fn column(&self, tag: &str) -> Option<impl Iterator<Item = &Value> + '_> {
        let i = self.tag_index(tag)?;
        Some(self.rows.iter().map(move |row| &row[i]))
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
