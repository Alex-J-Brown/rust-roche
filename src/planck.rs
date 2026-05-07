use pyo3::prelude::*;

// Speed of light, MKS, (exact value)
pub const C: f64 = 2.99792458e8;

// Planck's constant, MKS
pub const H: f64 = 6.6262e-34;

// Boltzmann's constant, MKS
pub const K: f64 = 1.3806e-23;

///
/// Computes the Planck function Bnu = (2 h \nu^3/c^2)/(exp(h \nu/kT) - 1)
///  as a function of wavelength and temperature. Output units are W/m**2/Hz/sr.
///
/// Arguments:
///
/// * `wave`: wavelength in nanometres
/// * `temp`: temperature in K
///
#[pyfunction]
pub fn planck(wave: f64, temp: f64) -> f64 {
    //
    // Uses Planck function in B-lambda form.
    // Returns specific intensity per unit wavelength of a blackbody
    // Units: W / m**3 / sr  (Wavelength converted to meters for mks unit system)
    //
    let lambda: f64 = wave * 1e-9;
    let x: f64 = (H * C) / (lambda * K * temp);
    let prefactor: f64 = (2.0 * H * C * C) / lambda.powi(5);

    if x > 40.0 {
        prefactor * (-x).exp()
    } else {
        prefactor / (x.exp() - 1.0)
    }
}

///
/// Computes the logarithmic derivative of the Planck function Bnu wrt
/// wavelength (i.e. d ln(Bnu) / d ln(lambda)) as a function of wavelength and temperature
///
/// Arguments:
///
/// * `wave`: wavelength in nanometres
/// * `temp`: temperature in K
///
#[pyfunction]
pub fn dplanck(wave: f64, temp: f64) -> f64 {
    //
    // Derivative of B-lambda wrt lambda.
    //
    let lambda: f64 = wave * 1e-9;
    let x: f64 = H * C / (lambda * K * temp);

    x / (1.0 - (-x).exp()) - 5.0
}

///
/// Computes the logarithmic derivative of the Planck function Bnu wrt
/// T (i.e. d ln(Bnu) / d ln(T)) as a function of wavelength and temperature
///
/// Arguments:
///
/// * `wave`: wavelength in nanometres
/// * `temp`: temperature in K
#[pyfunction]
pub fn dlpdlt(wave: f64, temp: f64) -> f64 {
    let fac2: f64 = 1.0e9 * H * C / K;

    let exponent: f64 = fac2 / (wave * temp);
    exponent / (1.0 - (-exponent).exp())
}
