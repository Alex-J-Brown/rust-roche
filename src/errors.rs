use pyo3::prelude::*;
use pyo3::exceptions::PyOSError;
use std::fmt;

#[derive(Debug)]
pub enum RocheError {
    // Dbrent failed to bracket minimum
    DbrentError,
    // Parameter error
    ParameterError(String),
}

impl std::error::Error for RocheError {}

impl fmt::Display for RocheError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RocheError::DbrentError => write!(f, "failed to bracket minimum with dbrent"),
            RocheError::ParameterError(msg) => write!(f, "parameter error: {}", msg),
        }
    }
}

impl std::convert::From<RocheError> for PyErr {
    fn from(err: RocheError) -> PyErr {
        PyOSError::new_err(err.to_string())
    }
}