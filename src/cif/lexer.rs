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
        let close = (0..bytes.len()).find(|&i| {
            bytes[i] == q
                && bytes
                    .get(i + 1)
                    .map_or(true, |&n| matches!(n as char, SP | HT | LF | CR))
        });
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
