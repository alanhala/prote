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
                assert!(values.len().is_multiple_of(cols), "loop values not a multiple of tags");

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

#[cfg(test)]
mod tests {
    use super::*;

    // 2.1 Data blocks and items

    #[test]
    fn single_block_single_item() {
        let cif = Cif::new("data_x\n_a 1");
        assert_eq!(cif.blocks.len(), 1);
        assert_eq!(cif.blocks[0].name, "x");
        assert_eq!(cif.blocks[0].items.len(), 1);
        assert!(matches!(
            &cif.blocks[0].items[0],
            Member::Item { tag, value: Value::Integer { value: 1, su: None } } if tag == "_a"
        ));
    }

    #[test]
    fn data_prefix_is_stripped_from_block_name() {
        let cif = Cif::new("data_4d1m\n_a 1");
        assert_eq!(cif.blocks[0].name, "4d1m");
    }

    #[test]
    fn multiple_blocks_parsed_in_order() {
        let cif = Cif::new("data_a\n_x 1\ndata_b\n_y 2");
        assert_eq!(cif.blocks.len(), 2);
        assert_eq!(cif.blocks[0].name, "a");
        assert_eq!(cif.blocks[1].name, "b");
    }

    #[test]
    fn item_value_types_round_trip() {
        let cif = Cif::new("data_x\n_i 5\n_f 1.5\n_s foo\n_dot .\n_q ?");
        let block = &cif.blocks[0];
        assert!(matches!(block.get("_i"), Some(Value::Integer { value: 5, su: None })));
        assert!(matches!(block.get("_f"), Some(Value::Float { value, su: None }) if *value == 1.5));
        assert!(matches!(block.get("_s"), Some(Value::Str(s)) if s == "foo"));
        assert!(matches!(block.get("_dot"), Some(Value::Inapplicable)));
        assert!(matches!(block.get("_q"), Some(Value::Unknown)));
    }

    #[test]
    fn quoted_string_value_strips_quotes() {
        let cif = Cif::new("data_x\n_a 'hi'");
        assert!(matches!(cif.blocks[0].get("_a"), Some(Value::Str(s)) if s == "hi"));
    }

    #[test]
    fn text_field_value_strips_delimiters() {
        let cif = Cif::new("data_x\n_a\n;text\n;");
        assert!(matches!(cif.blocks[0].get("_a"), Some(Value::Str(s)) if s == "text"));
    }

    #[test]
    fn content_before_first_data_heading_is_silently_dropped() {
        let cif = Cif::new("_a 1\ndata_x\n_b 2");
        assert!(cif.blocks.is_empty());
    }

    // 2.2 Loops

    #[test]
    fn well_formed_loop_two_tags_two_rows() {
        let cif = Cif::new("data_x\nloop_\n_a\n_b\n1 2\n3 4");
        let Member::Loop(l) = &cif.blocks[0].items[0] else {
            panic!("expected a loop");
        };
        assert_eq!(l.tags, vec!["_a".to_string(), "_b".to_string()]);
        assert_eq!(l.rows.len(), 2);
        assert!(matches!(l.rows[0][0], Value::Integer { value: 1, .. }));
        assert!(matches!(l.rows[0][1], Value::Integer { value: 2, .. }));
        assert!(matches!(l.rows[1][0], Value::Integer { value: 3, .. }));
        assert!(matches!(l.rows[1][1], Value::Integer { value: 4, .. }));
    }

    #[test]
    fn loop_with_one_tag_and_n_values() {
        let cif = Cif::new("data_x\nloop_\n_a\n1\n2\n3");
        let Member::Loop(l) = &cif.blocks[0].items[0] else {
            panic!("expected a loop");
        };
        assert_eq!(l.tags.len(), 1);
        assert_eq!(l.rows.len(), 3);
        assert!(l.rows.iter().all(|row| row.len() == 1));
    }

    #[test]
    #[should_panic(expected = "loop values not a multiple of tags")]
    fn loop_value_count_not_a_multiple_of_tag_count_panics() {
        Cif::new("data_x\nloop_\n_a\n_b\n1 2 3");
    }

    #[test]
    #[should_panic(expected = "loop_ has no values")]
    fn loop_with_tags_but_no_values_panics() {
        Cif::new("data_x\nloop_\n_a\n_b\n");
    }

    #[test]
    #[should_panic(expected = "loop_ has no tags")]
    fn loop_with_no_tags_panics() {
        Cif::new("data_x\nloop_\n1 2 3");
    }

    #[test]
    fn two_consecutive_loops_do_not_bleed_into_each_other() {
        let cif = Cif::new("data_x\nloop_\n_a\n1\n2\nloop_\n_b\n3\n4\n5");
        let items = &cif.blocks[0].items;
        assert_eq!(items.len(), 2);

        let Member::Loop(first) = &items[0] else {
            panic!("expected the first member to be a loop");
        };
        assert_eq!(first.tags, vec!["_a".to_string()]);
        assert_eq!(first.rows.len(), 2);
        assert!(matches!(first.rows[0][..], [Value::Integer { value: 1, .. }]));
        assert!(matches!(first.rows[1][..], [Value::Integer { value: 2, .. }]));

        let Member::Loop(second) = &items[1] else {
            panic!("expected the second member to be a loop");
        };
        assert_eq!(second.tags, vec!["_b".to_string()]);
        assert_eq!(second.rows.len(), 3);
        assert!(matches!(second.rows[0][..], [Value::Integer { value: 3, .. }]));
        assert!(matches!(second.rows[1][..], [Value::Integer { value: 4, .. }]));
        assert!(matches!(second.rows[2][..], [Value::Integer { value: 5, .. }]));
    }

    #[test]
    fn loop_ends_at_the_next_tag() {
        let cif = Cif::new("data_x\nloop_\n_a\n1\n2\n_b 5");
        let items = &cif.blocks[0].items;
        assert_eq!(items.len(), 2);
        let Member::Loop(l) = &items[0] else {
            panic!("expected a loop first");
        };
        assert_eq!(l.rows.len(), 2);
        assert!(matches!(l.rows[0][..], [Value::Integer { value: 1, su: None }]));
        assert!(matches!(l.rows[1][..], [Value::Integer { value: 2, su: None }]));
        assert!(matches!(
            &items[1],
            Member::Item { tag, value: Value::Integer { value: 5, su: None } } if tag == "_b"
        ));
    }

    // 2.3 Fail-loud behaviour on invalid input

    #[test]
    #[should_panic(expected = "tag must be followed by a value")]
    fn tag_with_no_following_value_panics() {
        Cif::new("data_x\n_a");
    }

    #[test]
    #[should_panic(expected = "UnterminatedString")]
    fn lexer_error_in_input_panics_parser_construction() {
        Cif::new("data_x\n_a 'oops");
    }

    #[test]
    fn save_frame_is_not_parsed_and_silently_ends_the_block() {
        let cif = Cif::new("data_x\n_a 1\nsave_s\n_b 2\nsave_\n_c 3");
        let block = &cif.blocks[0];
        assert_eq!(block.items.len(), 1);
        assert!(matches!(block.get("_a"), Some(Value::Integer { value: 1, .. })));
        assert!(block.get("_b").is_none());
        assert!(block.get("_c").is_none());
    }

    #[test]
    fn stop_token_silently_ends_the_block() {
        let cif = Cif::new("data_x\n_a 1\nstop_\n_b 2");
        let block = &cif.blocks[0];
        assert_eq!(block.items.len(), 1);
        assert!(block.get("_b").is_none());
    }
}
