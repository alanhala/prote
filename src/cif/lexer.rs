const SP: char = ' ';
const HT: char = '\t';
const LF: char = '\n';
const CR: char = '\r';

#[derive(Debug, Clone, Copy)]
pub struct Token<'a> {
    pub lexeme: &'a str,
    pub kind: TokenKind,
}

#[derive(Debug, Clone, Copy)]
pub enum TokenKind {
    Tag,
    Loop,
    Stop,
    Global,
    Inapplicable,
    Unknown,
    DataHeading,
    SaveEnd,
    SaveHeading,
    QuotedString,
    UnquotedString,
    TextField,
    Integer { value: i64, su: Option<u64> },
    Float { value: f64, su: Option<u64> },
}

#[derive(Debug)]
pub enum LexError {
    UnterminatedString { quote: char, at: usize },
    UnterminatedTextField { at: usize },
    EmptyTag { at: usize },
    BareDataHeading { at: usize },
}

impl std::fmt::Display for LexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            LexError::UnterminatedString { quote, at } => {
                write!(f, "unterminated {quote} string starting at byte {at}")
            }
            LexError::UnterminatedTextField { at } => {
                write!(f, "unterminated text field starting at byte {at}")
            }
            LexError::EmptyTag { at } => write!(f, "tag has no name at byte {at}"),
            LexError::BareDataHeading { at } => {
                write!(f, "bare `data_` with no block code at byte {at}")
            }
        }
    }
}

impl std::error::Error for LexError {}

pub struct Lexer<'a> {
    whole: &'a str,
    rest: &'a str,
    byte: usize,
}

impl<'a> Lexer<'a> {
    pub fn new(input: &'a str) -> Self {
        Self {
            whole: input,
            rest: input,
            byte: 0,
        }
    }

    fn classify_number(lexeme: &str) -> Option<TokenKind> {
        let (num, su) = match lexeme.split_once('(') {
            Some((num, rest)) => match rest.strip_suffix(')') {
                Some(su) => (num, Some(su)),
                None => return None,
            },
            None => (lexeme, None),
        };

        let su = match su {
            Some(s) => Some(s.parse::<u64>().ok()?),
            None => None,
        };

        if let Ok(value) = num.parse::<i64>() {
            return Some(TokenKind::Integer { value, su });
        }
        if let Ok(value) = num.parse::<f64>() {
            if value.is_finite() {
                return Some(TokenKind::Float { value, su });
            }
        }
        None
    }

    fn is_ordinary_char(c: char) -> bool {
        matches!(c,
            '!' | '%' | '&' | '(' | ')' | '*' | '+' | ',' | '-' | '.' | '/'
          | '0'..='9'
          | ':' | '<' | '=' | '>' | '?' | '@'
          | 'A'..='Z'
          | '\\' | '^' | '`'
          | 'a'..='z'
          | '{' | '|' | '}' | '~'
        )
    }

    fn is_non_blank_char(c: char) -> bool {
        Self::is_ordinary_char(c) || matches!(c, '"' | '#' | '$' | '\'' | '_' | ';' | '[' | ']')
    }

    fn at_line_start(&self, c_at: usize) -> bool {
        c_at == 0 || matches!(self.whole.as_bytes()[c_at - 1], b'\n' | b'\r')
    }

    fn take(&mut self, c_onwards: &'a str, c_at: usize, len: usize) -> &'a str {
        self.byte = c_at + len;
        self.rest = &c_onwards[len..];
        &c_onwards[..len]
    }

    fn skip_comment(&mut self, c_onwards: &'a str, c_at: usize) {
        let len = c_onwards.find('\n').unwrap_or(c_onwards.len());
        self.take(c_onwards, c_at, len);
    }

    fn lex_quoted(&mut self, quote: char, c_onwards: &'a str, c_at: usize) -> Result<Token<'a>, LexError> {
        let q = quote as u8;
        let bytes = self.rest.as_bytes();
        let close = (0..bytes.len())
            .find(|&i| bytes[i] == q && bytes.get(i + 1).is_none_or(|&n| matches!(n as char, SP | HT | LF | CR)));
        match close {
            None => Err(LexError::UnterminatedString { quote, at: c_at }),
            Some(rel) => {
                // token = opening quote + `rel` content bytes + closing quote
                let lexeme = self.take(c_onwards, c_at, rel + 2);
                Ok(Token {
                    lexeme,
                    kind: TokenKind::QuotedString,
                })
            }
        }
    }

    fn lex_tag(&mut self, c_onwards: &'a str, c_at: usize) -> Result<Token<'a>, LexError> {
        let end = c_onwards
            .find(|ch| !Self::is_non_blank_char(ch))
            .unwrap_or(c_onwards.len());
        if end == 1 {
            return Err(LexError::EmptyTag { at: c_at });
        }
        let lexeme = self.take(c_onwards, c_at, end);
        Ok(Token {
            lexeme,
            kind: TokenKind::Tag,
        })
    }

    fn lex_text_field(&mut self, c_onwards: &'a str, c_at: usize) -> Result<Token<'a>, LexError> {
        match c_onwards.find("\n;") {
            None => Err(LexError::UnterminatedTextField { at: c_at }),
            Some(nl) => {
                let lexeme = self.take(c_onwards, c_at, nl + "\n;".len());
                Ok(Token {
                    lexeme,
                    kind: TokenKind::TextField,
                })
            }
        }
    }

    fn lex_value(&mut self, c_onwards: &'a str, c_at: usize) -> Result<Token<'a>, LexError> {
        let end = c_onwards
            .find(|ch: char| !Self::is_non_blank_char(ch))
            .unwrap_or(c_onwards.len());
        let lexeme = self.take(c_onwards, c_at, end);
        let kind = if lexeme.eq_ignore_ascii_case("loop_") {
            TokenKind::Loop
        } else if lexeme.eq_ignore_ascii_case("stop_") {
            TokenKind::Stop
        } else if lexeme.eq_ignore_ascii_case("global_") {
            TokenKind::Global
        } else if lexeme == "." {
            TokenKind::Inapplicable
        } else if lexeme == "?" {
            TokenKind::Unknown
        } else if lexeme.eq_ignore_ascii_case("data_") {
            return Err(LexError::BareDataHeading { at: c_at });
        } else if lexeme.len() >= 5 && lexeme[..5].eq_ignore_ascii_case("data_") {
            TokenKind::DataHeading
        } else if lexeme.eq_ignore_ascii_case("save_") {
            TokenKind::SaveEnd
        } else if lexeme.len() >= 5 && lexeme[..5].eq_ignore_ascii_case("save_") {
            TokenKind::SaveHeading
        } else if let Some(kind) = Self::classify_number(lexeme) {
            kind
        } else {
            TokenKind::UnquotedString
        };
        Ok(Token { lexeme, kind })
    }
}

impl<'a> Iterator for Lexer<'a> {
    type Item = Result<Token<'a>, LexError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mut chars = self.rest.chars();
            let c = chars.next()?;
            let c_at = self.byte;
            let c_onwards = self.rest;
            self.rest = chars.as_str();
            self.byte += c.len_utf8();

            match c {
                HT | SP | LF | CR => continue,
                '#' => self.skip_comment(c_onwards, c_at),
                '\'' | '"' => return Some(self.lex_quoted(c, c_onwards, c_at)),
                '_' => return Some(self.lex_tag(c_onwards, c_at)),
                ';' if self.at_line_start(c_at) => return Some(self.lex_text_field(c_onwards, c_at)),
                _ => return Some(self.lex_value(c_onwards, c_at)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokens(input: &str) -> Vec<Token<'_>> {
        Lexer::new(input).map(Result::unwrap).collect()
    }

    fn results(input: &str) -> Vec<Result<Token<'_>, LexError>> {
        Lexer::new(input).collect()
    }

    // 1.1 Whitespace, comments, and empty input

    #[test]
    fn empty_input_yields_no_tokens() {
        assert_eq!(tokens("").len(), 0);
    }

    #[test]
    fn whitespace_only_yields_no_tokens() {
        assert_eq!(tokens("  \t\n\r\n ").len(), 0);
    }

    #[test]
    fn comment_line_is_skipped() {
        assert_eq!(tokens("# just a comment\n").len(), 0);
    }

    #[test]
    fn comment_terminated_by_eof_is_skipped() {
        assert_eq!(tokens("# comment, no trailing newline").len(), 0);
    }

    #[test]
    fn trailing_comment_after_token() {
        let toks = tokens("foo # bar");
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, "foo");
        assert!(matches!(toks[0].kind, TokenKind::UnquotedString));
    }

    #[test]
    fn hash_mid_token_is_not_a_comment() {
        let toks = tokens("ab#cd");
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, "ab#cd");
        assert!(matches!(toks[0].kind, TokenKind::UnquotedString));
    }

    // 1.2 Number classification

    #[test]
    fn plain_integer() {
        let toks = tokens("42");
        assert!(matches!(toks[0].kind, TokenKind::Integer { value: 42, su: None }));
    }

    #[test]
    fn signed_integers() {
        assert!(matches!(tokens("-7")[0].kind, TokenKind::Integer { value: -7, su: None }));
        assert!(matches!(tokens("+7")[0].kind, TokenKind::Integer { value: 7, su: None }));
    }

    #[test]
    fn plain_float() {
        assert!(matches!(tokens("3.14")[0].kind, TokenKind::Float { value, su: None } if value == 3.14));
    }

    #[test]
    fn leading_and_trailing_dot_floats() {
        for input in [".5", "2.", "+.5"] {
            let toks = tokens(input);
            assert_eq!(toks.len(), 1, "input {input:?} should be a single token");
            assert!(
                matches!(toks[0].kind, TokenKind::Float { .. }),
                "input {input:?} should classify as Float, got {:?}",
                toks[0].kind
            );
        }
    }

    #[test]
    fn exponent_floats() {
        assert!(matches!(tokens("1e10")[0].kind, TokenKind::Float { .. }));
        assert!(matches!(tokens("-2.5E-3")[0].kind, TokenKind::Float { .. }));
    }

    #[test]
    fn integer_with_standard_uncertainty() {
        assert!(matches!(tokens("12(3)")[0].kind, TokenKind::Integer { value: 12, su: Some(3) }));
    }

    #[test]
    fn float_with_standard_uncertainty() {
        assert!(matches!(
            tokens("1.5(3)")[0].kind,
            TokenKind::Float { value, su: Some(3) } if value == 1.5
        ));
    }

    #[test]
    fn non_numeric_su_falls_back_to_unquoted_string() {
        let toks = tokens("1(x)");
        assert_eq!(toks[0].lexeme, "1(x)");
        assert!(matches!(toks[0].kind, TokenKind::UnquotedString));
    }

    #[test]
    fn unmatched_paren_falls_back_to_unquoted_string() {
        let toks = tokens("1(2");
        assert_eq!(toks[0].lexeme, "1(2");
        assert!(matches!(toks[0].kind, TokenKind::UnquotedString));
    }

    #[test]
    fn trailing_letters_are_not_a_number() {
        let toks = tokens("12abc");
        assert_eq!(toks[0].lexeme, "12abc");
        assert!(matches!(toks[0].kind, TokenKind::UnquotedString));
    }

    #[test]
    fn number_fallback_cases_are_unquoted_strings() {
        for input in ["1.2.3", "+", "-", "1-2"] {
            let toks = tokens(input);
            assert_eq!(toks.len(), 1, "input {input:?} should be a single token");
            assert!(
                matches!(toks[0].kind, TokenKind::UnquotedString),
                "input {input:?} should classify as UnquotedString, got {:?}",
                toks[0].kind
            );
        }
    }

    #[test]
    fn integer_overflow_promotes_to_float() {
        assert!(matches!(tokens("99999999999999999999")[0].kind, TokenKind::Float { .. }));
    }

    #[test]
    fn inf_and_nan_are_not_finite_numbers() {
        for input in ["inf", "NaN"] {
            let toks = tokens(input);
            assert!(
                matches!(toks[0].kind, TokenKind::UnquotedString),
                "input {input:?} should be rejected by the is_finite guard, got {:?}",
                toks[0].kind
            );
        }
    }

    #[test]
    fn lone_dot_and_question_mark_are_not_numbers() {
        assert!(matches!(tokens(".")[0].kind, TokenKind::Inapplicable));
        assert!(matches!(tokens("?")[0].kind, TokenKind::Unknown));
    }

    // 1.3 Quoted strings

    #[test]
    fn simple_single_quoted_string() {
        let toks = tokens("'hello'");
        assert_eq!(toks[0].lexeme, "'hello'");
        assert!(matches!(toks[0].kind, TokenKind::QuotedString));
    }

    #[test]
    fn primed_atom_name_in_double_quotes() {
        let toks = tokens("\"O5'\"");
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, "\"O5'\"");
        assert!(matches!(toks[0].kind, TokenKind::QuotedString));
    }

    #[test]
    fn non_terminating_apostrophe_inside_single_quotes() {
        let toks = tokens("'it's fine'");
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, "'it's fine'");
        assert!(matches!(toks[0].kind, TokenKind::QuotedString));
    }

    #[test]
    fn closing_quote_at_eof() {
        let toks = tokens("'end'");
        assert_eq!(toks[0].lexeme, "'end'");
        assert!(matches!(toks[0].kind, TokenKind::QuotedString));
    }

    #[test]
    fn double_quoted_string() {
        let toks = tokens("\"a b c\"");
        assert_eq!(toks[0].lexeme, "\"a b c\"");
        assert!(matches!(toks[0].kind, TokenKind::QuotedString));
    }

    #[test]
    fn unterminated_quote_is_an_error() {
        let results = results("'unterminated");
        assert!(matches!(results[0], Err(LexError::UnterminatedString { quote: '\'', .. })));
    }

    #[test]
    fn empty_quoted_string() {
        let toks = tokens("''");
        assert_eq!(toks[0].lexeme, "''");
        assert!(matches!(toks[0].kind, TokenKind::QuotedString));
    }

    // 1.4 Tags

    #[test]
    fn plain_tag() {
        let toks = tokens("_atom_site.Cartn_x");
        assert_eq!(toks[0].lexeme, "_atom_site.Cartn_x");
        assert!(matches!(toks[0].kind, TokenKind::Tag));
    }

    #[test]
    fn bare_underscore_before_whitespace_is_empty_tag() {
        assert!(matches!(results("_ ")[0], Err(LexError::EmptyTag { .. })));
    }

    #[test]
    fn bare_underscore_at_eof_is_empty_tag() {
        assert!(matches!(results("_")[0], Err(LexError::EmptyTag { .. })));
    }

    #[test]
    fn tag_terminated_by_whitespace_newline_or_eof() {
        assert_eq!(tokens("_a b")[0].lexeme, "_a");
        assert_eq!(tokens("_a\nb")[0].lexeme, "_a");
        assert_eq!(tokens("_a")[0].lexeme, "_a");
    }

    #[test]
    fn tag_with_bracket_characters() {
        let toks = tokens("_foo[1]");
        assert_eq!(toks[0].lexeme, "_foo[1]");
        assert!(matches!(toks[0].kind, TokenKind::Tag));
    }

    // 1.5 Text fields

    #[test]
    fn multiline_text_field() {
        let toks = tokens(";line one\nline two\n;");
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, ";line one\nline two\n;");
        assert!(matches!(toks[0].kind, TokenKind::TextField));
    }

    #[test]
    fn semicolon_mid_line_is_not_a_text_field_opener() {
        let toks = tokens("x ;not a field");
        assert!(!toks.iter().any(|t| matches!(t.kind, TokenKind::TextField)));
        assert_eq!(toks[0].lexeme, "x");
    }

    #[test]
    fn unterminated_text_field_is_an_error() {
        let results = results(";opened but never closed");
        assert!(matches!(results[0], Err(LexError::UnterminatedTextField { .. })));
    }

    #[test]
    fn content_line_with_leading_space_before_semicolon_is_not_a_terminator() {
        // Only a bol `;` (no preceding space) closes the field, per the tokens
        // EBNF; `find("\n;")` naturally enforces that, so " ;still content"
        // stays inside the field instead of ending it early.
        let input = ";first\n ;still content\nreal end\n;";
        let toks = tokens(input);
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, input);
        assert!(matches!(toks[0].kind, TokenKind::TextField));
    }

    #[test]
    fn text_field_captures_blank_lines_and_hashes_verbatim() {
        let input = ";alpha\n\n# not a comment\nbeta\n;";
        let toks = tokens(input);
        assert_eq!(toks.len(), 1);
        assert_eq!(toks[0].lexeme, input);
        assert!(matches!(toks[0].kind, TokenKind::TextField));
    }

    #[test]
    fn text_field_at_byte_zero_is_recognized() {
        let toks = tokens(";x\n;");
        assert!(matches!(toks[0].kind, TokenKind::TextField));
    }

    // 1.6 Keyword / reserved-word classification

    #[test]
    fn loop_keyword_is_case_insensitive() {
        for input in ["loop_", "LOOP_", "Loop_"] {
            assert!(matches!(tokens(input)[0].kind, TokenKind::Loop), "input {input:?}");
        }
    }

    #[test]
    fn stop_and_global_keywords() {
        assert!(matches!(tokens("stop_")[0].kind, TokenKind::Stop));
        assert!(matches!(tokens("global_")[0].kind, TokenKind::Global));
    }

    #[test]
    fn bare_data_heading_is_an_error() {
        assert!(matches!(results("data_")[0], Err(LexError::BareDataHeading { .. })));
    }

    #[test]
    fn data_heading_with_block_code() {
        let toks = tokens("data_4d1m");
        assert_eq!(toks[0].lexeme, "data_4d1m");
        assert!(matches!(toks[0].kind, TokenKind::DataHeading));
    }

    #[test]
    fn bare_save_end() {
        assert!(matches!(tokens("save_")[0].kind, TokenKind::SaveEnd));
    }

    #[test]
    fn save_heading_with_frame_code() {
        let toks = tokens("save_foo");
        assert_eq!(toks[0].lexeme, "save_foo");
        assert!(matches!(toks[0].kind, TokenKind::SaveHeading));
    }

    #[test]
    fn keyword_prefix_matching_is_case_insensitive() {
        assert!(matches!(tokens("Data_Foo")[0].kind, TokenKind::DataHeading));
        assert!(matches!(tokens("SAVE_Bar")[0].kind, TokenKind::SaveHeading));
    }

    #[test]
    fn tag_and_loop_keyword_do_not_collide() {
        assert!(matches!(tokens("_loop")[0].kind, TokenKind::Tag));
        assert!(matches!(tokens("loop_")[0].kind, TokenKind::Loop));
    }

    #[test]
    fn reserved_word_lookalikes_are_unquoted_strings() {
        for input in ["loopy_", "datax", "data", "savely_"] {
            let toks = tokens(input);
            assert_eq!(toks.len(), 1, "input {input:?}");
            assert!(
                matches!(toks[0].kind, TokenKind::UnquotedString),
                "input {input:?} should be UnquotedString, got {:?}",
                toks[0].kind
            );
        }
    }

    #[test]
    fn embedded_question_mark_and_dot_stay_unquoted_strings() {
        for input in ["a?b", "1.2.3"] {
            let toks = tokens(input);
            assert!(matches!(toks[0].kind, TokenKind::UnquotedString), "input {input:?}");
        }
    }

    #[test]
    fn error_byte_offsets_match_real_position() {
        let results = results("abc _ ");
        match &results[1] {
            Err(LexError::EmptyTag { at }) => assert_eq!(*at, 4),
            other => panic!("expected EmptyTag at byte 4, got {other:?}"),
        }
    }
}
