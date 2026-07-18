use crate::cif::lexer::TokenKind;

use super::ast::{Cif, DataBlock, Loop, Member, Value};
use super::lexer::{Lexer, Token};

pub struct Parser<'a> {
    tokens: Vec<Token<'a>>,
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Parser {
            tokens: Lexer::new(input).map(Result::unwrap).collect(),
            pos: 0,
        }
    }

    fn peek_kind(&self) -> Option<TokenKind> {
        self.tokens.get(self.pos).map(|t| t.kind)
    }

    fn consume(&mut self) -> Token<'a> {
        let token = self.tokens[self.pos];
        self.pos += 1;
        token
    }

    fn parse_cif(&mut self) -> Cif {
        let mut blocks = Vec::new();
        while let Some(block) = self.try_data_block() {
            blocks.push(block);
        }
        Cif { blocks }
    }

    fn try_data_block(&mut self) -> Option<DataBlock> {
        let name = self.try_data_heading()?;
        let mut items = Vec::new();
        while let Some(member) = self.try_member() {
            items.push(member);
        }
        Some(DataBlock { name, items })
    }

    fn try_data_heading(&mut self) -> Option<String> {
        match self.peek_kind()? {
            TokenKind::DataHeading => Some(self.consume().lexeme[5..].to_string()), // drop "data_"
            _ => None,
        }
    }

    fn try_member(&mut self) -> Option<Member> {
        self.try_item().or_else(|| self.try_loop().map(Member::Loop))
    }

    fn try_item(&mut self) -> Option<Member> {
        let tag = self.try_tag()?;
        let value = self.try_value().expect("tag must be followed by a value");
        Some(Member::Item { tag, value })
    }

    fn try_loop(&mut self) -> Option<Loop> {
        match self.peek_kind()? {
            TokenKind::Loop => {
                self.consume();
                let mut tags: Vec<String> = vec![];
                while let Some(tag) = self.try_tag() {
                    tags.push(tag);
                }
                let mut values: Vec<Value> = vec![];
                while let Some(value) = self.try_value() {
                    values.push(value);
                }

                let cols = tags.len();
                assert!(cols > 0, "loop_ has no tags");
                assert!(!values.is_empty(), "loop_ has no values");
                assert!(values.len() % cols == 0, "loop values not a multiple of tags");

                let mut it = values.into_iter();
                let mut rows = Vec::new();
                loop {
                    let row: Vec<Value> = it.by_ref().take(cols).collect();
                    if row.is_empty() {
                        break;
                    }
                    rows.push(row);
                }

                Some(Loop { tags, rows })
            }
            _ => None,
        }
    }

    fn try_tag(&mut self) -> Option<String> {
        match self.peek_kind()? {
            TokenKind::Tag => Some(self.consume().lexeme.to_string()),
            _ => None,
        }
    }

    /// value ::= INTEGER | FLOAT | strings | '.' | '?'
    fn try_value(&mut self) -> Option<Value> {
        match self.peek_kind()? {
            TokenKind::Integer { value, su } => {
                self.consume();
                Some(Value::Integer { value, su })
            }
            TokenKind::Float { value, su } => {
                self.consume();
                Some(Value::Float { value, su })
            }
            TokenKind::UnquotedString => {
                let token = self.consume();
                Some(Value::Str(token.lexeme.to_string()))
            }
            TokenKind::QuotedString => {
                let token = self.consume();
                Some(Value::Str(token.lexeme[1..token.lexeme.len() - 1].to_string()))
            }
            TokenKind::TextField => {
                let token = self.consume();
                Some(Value::Str(token.lexeme[1..token.lexeme.len() - 2].to_string()))
            }
            TokenKind::Inapplicable => {
                self.consume();
                Some(Value::Inapplicable)
            }
            TokenKind::Unknown => {
                self.consume();
                Some(Value::Unknown)
            }
            _ => None,
        }
    }
}

impl Cif {
    pub fn new(input: &str) -> Cif {
        Parser::new(input).parse_cif()
    }
}
