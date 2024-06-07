use std::fmt::{ self, Display };

#[derive(Clone, Debug)]
pub struct NotSatisfiedError;

impl Display for NotSatisfiedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Not Satisfied")
    }
}

impl std::error::Error for NotSatisfiedError {}
