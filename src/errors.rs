use pyo3::prelude::*;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use std::fmt;

#[derive(Debug)]
pub enum RocheError {
    // Dbrent failed to bracket minimum
    DbrentError,
    // Parameter error
    ParameterError(String),
    // error in Face function
    FaceError(String),
    // error in rtsafe function
    RtsafeError(String)
}

impl std::error::Error for RocheError {}

impl fmt::Display for RocheError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RocheError::DbrentError => write!(f, "failed to bracket minimum with dbrent"),
            RocheError::ParameterError(msg) => write!(f, "{}", msg),
            RocheError::FaceError(msg) => write!(f, "{}", msg),
            RocheError::RtsafeError(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::convert::From<RocheError> for PyErr {
    fn from(err: RocheError) -> PyErr {
        match err {
            RocheError::DbrentError => PyRuntimeError::new_err(err.to_string()),
            RocheError::ParameterError(_) => PyValueError::new_err(err.to_string()),
            RocheError::FaceError(_) => PyRuntimeError::new_err(err.to_string()),
            RocheError::RtsafeError(_) => PyRuntimeError::new_err(err.to_string()), 
        }
    }
}