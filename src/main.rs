use std::fs;

fn main() {
    let cif_file = fs::read_to_string("4d1m.cif").expect("Could not open the file");
    let mut lexer = Lexer::new(cif_file);
    let asd = lexer.next();
    println!("{:?}", asd);
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
    type Item = Token;

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
                TransitionType::EmitToken { token, end } => {
                    self.state = transition.state;
                    self.cursor.align();
                    self.end = end;
                    return Some(token);
                }
                TransitionType::End => {
                    self.state = transition.state;
                    return None;
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
struct ErrorState {
    error: String,
}
struct EndState;

struct Transition {
    state: Box<dyn State>,
    transition_type: TransitionType,
}
enum TransitionType {
    EmitToken { token: Token, end: bool },
    Advance,
    End,
}

// TODO: Define how to end the token
impl State for CommentState {
    fn consume(&self, cursor: &Cursor) -> Transition {
        match cursor.peek() {
            None => Transition {
                state: Box::new(CommentState),
                transition_type: TransitionType::EmitToken {
                    token: Token {
                        lexeme: cursor.lexeme(),
                        kind: TokenKind::Comment,
                    },
                    end: true,
                },
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
                    transition_type: TransitionType::EmitToken {
                        token: Token {
                            lexeme: cursor.lexeme(),
                            kind: TokenKind::Comment,
                        },
                        end: false,
                    },
                },
                _ => todo!(), // # Error
            },
        }
    }
}

impl State for StartState {
    fn consume(&self, cursor: &Cursor) -> Transition {
        match cursor.peek() {
            None => Transition {
                state: Box::new(StartState),
                transition_type: TransitionType::End,
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

// impl State for EndState {
//     fn consume(&self, char: char) -> Transition {
//         Transition {
//             state: Box::new(EndState),
//             transition_type: TransitionType::End,
//         }
//     }
// }
