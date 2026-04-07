use crate::Vec3;
use pyo3::prelude::*;
use std::f64::consts::PI;

///
/// set_earth_iangle computes the earth vector given an inclination and orbital phase.
///
/// Arguments:
///
/// * `iangle`: orbital inclination
/// * `phase`: orbital phase (1 = one orbit)
///
/// Returns:
///
/// * the earth vector (unit vector)
///
#[pyfunction]
pub fn set_earth_iangle(iangle: f64, phase: f64) -> Vec3 {
    let iangle_rad: f64 = iangle.to_radians();
    let phase_rad: f64 = 2.0 * PI * phase;
    let (sini, cosi) = iangle_rad.sin_cos();
    let (sinp, cosp) = phase_rad.sin_cos();
    Vec3 {
        x: sini * cosp,
        y: -sini * sinp,
        z: cosi,
    }
}

///
/// set_earth computes the earth vector given an inclination and orbital phase.
///
/// Arguments:
///
/// * `cosi`: cosine of orbital inclination
/// * `sini`: sine of orbital inclination
/// * `phase`: orbital phase (1 = one orbit)
///
/// Returns:
///
///  * the earth vector (unit vector)
///
#[pyfunction]
pub fn set_earth(cosi: f64, sini: f64, phase: f64) -> Vec3 {
    let phase_rad: f64 = 2.0 * PI * phase;
    let (sinp, cosp) = phase_rad.sin_cos();
    Vec3 {
        x: sini * cosp,
        y: -sini * sinp,
        z: cosi,
    }
}
