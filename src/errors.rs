use pyo3::prelude::*;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use std::fmt;

#[derive(Debug)]
pub enum RocheError {
    // error in Dbrent function
    DbrentError(String),
    // error in lin_min function
    LinminError(String),
    // error in pot_min function
    PotminError(String),
    // Parameter error
    ParameterError(String),
    // error in Face function
    FaceError(String),
    // error in rtsafe function
    RtsafeError(String),
    // error in wd_phases function
    WdphasesError(String)
}

impl std::error::Error for RocheError {}

impl fmt::Display for RocheError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RocheError::DbrentError(msg) => write!(f, "{}", msg),
            RocheError::LinminError(msg) => write!(f, "{}", msg),
            RocheError::PotminError(msg) => write!(f, "{}", msg),
            RocheError::ParameterError(msg) => write!(f, "{}", msg),
            RocheError::FaceError(msg) => write!(f, "{}", msg),
            RocheError::RtsafeError(msg) => write!(f, "{}", msg),
            RocheError::WdphasesError(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::convert::From<RocheError> for PyErr {
    fn from(err: RocheError) -> PyErr {
        match err {
            RocheError::DbrentError(_) => PyRuntimeError::new_err(err.to_string()),
            RocheError::LinminError(_) => PyRuntimeError::new_err(err.to_string()),
            RocheError::PotminError(_) => PyRuntimeError::new_err(err.to_string()),
            RocheError::ParameterError(_) => PyValueError::new_err(err.to_string()),
            RocheError::FaceError(_) => PyRuntimeError::new_err(err.to_string()),
            RocheError::RtsafeError(_) => PyRuntimeError::new_err(err.to_string()), 
            RocheError::WdphasesError(_) => PyRuntimeError::new_err(err.to_string()), 
        }
    }
}