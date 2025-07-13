use std::fs;

fn main() {
    let cif_file = fs::read_to_string("4d1m.cif").unwrap();
    let mut lexer = Lexer::new(cif_file);
    for result in lexer {
        match result {
            Ok(token) => println!("Token {:?}", token),
            Err(lexer_error) => println!("Lexer Error"),
        }
    }
}

#[derive(Debug)]
struct Token {
    lexeme: String,
    kind: TokenKind,
}

#[derive(Debug)]
enum TokenKind {
    Comment,
    Integer,
    Float,
}

#[derive(Debug)]
struct Cursor {
    input: String,
    line: usize,
    position: usize,
    offset: usize,
    buffer: Vec<char>,
}

#[derive(Debug)]
struct StateInput {
    char: char,
    buffer: Vec<char>,
}

impl Cursor {
    fn new(input: String) -> Cursor {
        Cursor {
            input,
            line: 0,
            position: 0,
            offset: 0,
            buffer: vec![],
        }
    }

    fn peek(&self) -> Option<char> {
        self.input.chars().nth(self.position + self.offset)
    }

    fn advance(&mut self) {
        self.buffer.push(self.peek().unwrap());
        self.offset += 1;
    }

    fn new_line(&mut self) {
        self.line += 1;
    }

    fn start(&mut self) {
        self.position += self.offset;
        self.offset = 0;
        self.buffer = vec![];
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
enum State {
    StartLine,
    Start,
    Comment,
    Integer,
    DecimalPoint,
    Float,
}

impl State {
    fn transition(self, input: &StateInput) -> Result<(State, LexerAction), LexerError> {
        // println!("Processing {:?}", input);
        match self {
            State::StartLine => match input.char {
                ';' => todo!(),
                _ => Ok((State::Start, LexerAction::Advance { consume: false })),
            },

            State::Start => match input.char {
                '#' => Ok((State::Comment, LexerAction::Advance { consume: true })),
                '\n' => Ok((State::StartLine, LexerAction::Advance { consume: true })),
                ' ' => Ok((State::Start, LexerAction::Advance { consume: true })),
                char if char.is_ascii_digit() => {
                    Ok((State::Integer, LexerAction::Advance { consume: true }))
                }
                _ => Err(LexerError),
            },

            State::Comment => match input.char {
                c if c.is_alphanumeric()
                    || c.is_ascii_punctuation()
                    || (c.is_whitespace() && c != '\n') =>
                {
                    Ok((State::Comment, LexerAction::Advance { consume: true }))
                }
                '\n' => Ok((
                    State::StartLine,
                    LexerAction::Emit {
                        token: Token {
                            lexeme: String::from_iter(input.buffer.clone()),
                            kind: TokenKind::Comment,
                        },
                        consume: false,
                    },
                )),
                _ => Err(LexerError),
            },

            State::Integer => match input.char {
                char if char.is_ascii_digit() => {
                    Ok((State::Integer, LexerAction::Advance { consume: true }))
                }
                char if char.is_ascii_whitespace() => Ok((
                    State::Start,
                    LexerAction::Emit {
                        token: Token {
                            lexeme: String::from_iter(input.buffer.clone()),
                            kind: TokenKind::Integer,
                        },
                        consume: false,
                    },
                )),
                '.' => Ok((State::DecimalPoint, LexerAction::Advance { consume: true })),
                _ => Err(LexerError),
            },

            State::DecimalPoint => match input.char {
                char if char.is_ascii_digit() => {
                    Ok((State::Float, LexerAction::Advance { consume: true }))
                }
                _ => Err(LexerError),
            },

            State::Float => match input.char {
                char if char.is_ascii_digit() => {
                    Ok((State::Float, LexerAction::Advance { consume: true }))
                }
                char if char.is_ascii_whitespace() => Ok((
                    State::Start,
                    LexerAction::Emit {
                        token: Token {
                            lexeme: String::from_iter(input.buffer.clone()),
                            kind: TokenKind::Float,
                        },
                        consume: false,
                    },
                )),
                _ => Err(LexerError),
            },
        }
    }

    fn finish(self, acc: String) -> Result<Token, LexerError> {
        match self {
            State::Comment => Ok(Token {
                lexeme: acc,
                kind: TokenKind::Comment,
            }),
            _ => Err(LexerError),
        }
    }
}

#[derive(Debug)]
struct LexerError;

enum LexerAction {
    Emit { token: Token, consume: bool },
    Advance { consume: bool },
}

struct Lexer {
    cursor: Cursor,
    state: State,
    end: bool,
}

impl Lexer {
    fn new(input: String) -> Lexer {
        Lexer {
            cursor: Cursor::new(input),
            state: State::StartLine,
            end: false,
        }
    }
}

impl Iterator for Lexer {
    type Item = Result<Token, LexerError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.end {
            return None;
        }

        loop {
            // println!("Cursor: {:?}", self.cursor);
            // println!("State: {:?}", self.state);
            if self.state == State::StartLine {
                self.cursor.new_line();
            }
            if self.state == State::Start {
                self.cursor.start();
            }
            match self.cursor.peek() {
                None => {
                    self.end = true;
                    return Some(
                        self.state
                            .finish(String::from_iter(self.cursor.buffer.clone())),
                    );
                }
                Some(ch) => {
                    let input = StateInput {
                        char: ch,
                        buffer: self.cursor.buffer.clone(), // TODO: Try to not use clone,
                    };
                    match self.state.transition(&input) {
                        Ok((new_state, action)) => {
                            self.state = new_state;
                            match action {
                                LexerAction::Emit { token, consume } => {
                                    if consume {
                                        self.cursor.advance();
                                    }
                                    self.cursor.start();

                                    return Some(Ok(token));
                                }
                                LexerAction::Advance { consume } => {
                                    if consume {
                                        self.cursor.advance();
                                    }
                                }
                            }
                        }
                        Err(e) => {
                            self.end = true;
                            return Some(Err(e));
                        }
                    }
                }
            }
        }
    }
}
