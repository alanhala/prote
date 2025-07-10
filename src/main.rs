use std::{error::Error, fmt::Display, fs};

fn main() {
    let cif_file = fs::read_to_string("4d1m.cif").expect("Could not open the file");
    let mut lexer = Lexer::new(cif_file);
    for next in lexer.next() {
        match next {
            Ok(token) => println!("Token {:?}", token),
            Err(error) => println!("Lexer error at {}", error.position),
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
}

struct Cursor {
    input: String,
    position: usize,
    offset: usize,
}

impl Cursor {
    fn new(input: String) -> Cursor {
        Cursor {
            input: input,
            position: 0,
            offset: 0,
        }
    }

    fn peek(&self) -> Option<char> {
        self.input.chars().nth(self.offset)
    }

    fn lexeme(&self) -> String {
        self.input[self.position..self.offset].to_string()
    }

    fn align(&mut self) {
        self.position = self.offset;
        self.offset = 0;
    }

    fn advace(&mut self) {
        self.offset += 1
    }
}

struct Lexer {
    cursor: Cursor,
    state: Box<dyn State>,
    end: bool,
}

impl Iterator for Lexer {
    type Item = Result<Token, LexerError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.end {
            return None;
        }

        loop {
            let transition = self.state.consume(&self.cursor);
            match transition.transition_type {
                TransitionType::Advance => {
                    self.state = transition.state;
                    self.cursor.advace();
                }
                TransitionType::EmitToken(token) => {
                    self.state = transition.state;
                    self.cursor.align();
                    return Some(Ok(token));
                }
                TransitionType::End(token) => {
                    self.state = transition.state;
                    self.end = true;

                    return match token {
                        None => None,
                        Some(token) => Some(Ok(token)),
                    };
                }
                TransitionType::Error(error) => {
                    self.state = transition.state;
                    self.end = true;

                    return Some(Err(error));
                }
            }
        }
    }
}

impl Lexer {
    fn new(input: String) -> Lexer {
        Lexer {
            cursor: Cursor::new(input),
            state: Box::new(StartState),
            end: false,
        }
    }
}

trait State {
    fn consume(&self, cursor: &Cursor) -> Transition;
}

struct StartState;
struct CommentState;
struct ErrorState;

struct Transition {
    state: Box<dyn State>,
    transition_type: TransitionType,
}
enum TransitionType {
    EmitToken(Token),
    Advance,
    End(Option<Token>),
    Error(LexerError),
}

#[derive(Debug)]
struct LexerError {
    position: usize,
}

impl State for ErrorState {
    fn consume(&self, cursor: &Cursor) -> Transition {
        Transition {
            state: Box::new(ErrorState),
            transition_type: TransitionType::End(None),
        }
    }
}

impl State for CommentState {
    fn consume(&self, cursor: &Cursor) -> Transition {
        match cursor.peek() {
            None => Transition {
                state: Box::new(CommentState),
                transition_type: TransitionType::End(Some(Token {
                    lexeme: cursor.lexeme(),
                    kind: TokenKind::Comment,
                })),
            },
            Some(c) => match c {
                char if char.is_alphanumeric()
                    || char.is_ascii_punctuation()
                    || char.is_whitespace() =>
                {
                    Transition {
                        state: Box::new(CommentState),
                        transition_type: TransitionType::Advance,
                    }
                }
                '\n' => Transition {
                    state: Box::new(StartState),
                    transition_type: TransitionType::EmitToken(Token {
                        lexeme: cursor.lexeme(),
                        kind: TokenKind::Comment,
                    }),
                },
                _ => Transition {
                    state: Box::new(ErrorState),
                    transition_type: TransitionType::Error(LexerError {
                        position: cursor.position,
                    }),
                },
            },
        }
    }
}

impl State for StartState {
    fn consume(&self, cursor: &Cursor) -> Transition {
        match cursor.peek() {
            None => Transition {
                state: Box::new(StartState),
                transition_type: TransitionType::End(None),
            },
            Some(c) => match c {
                '#' => Transition {
                    state: Box::new(CommentState),
                    transition_type: TransitionType::Advance,
                },
                _ => todo!(),
            },
        }
    }
}
