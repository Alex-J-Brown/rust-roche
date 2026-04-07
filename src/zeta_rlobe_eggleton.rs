use crate::errors::RocheError;
use pyo3::prelude::*;
///
/// zeta_rlobe_eggleton returns d log(rl) / d log (m2)
/// where rl is Peter Eggleton's formula for the volume-averaged
/// Roche lobe radius divided by the orbital separation. This assumes
/// that m1+m2 = constant.
///
/// Arguments:
///
/// * `q`: mass ratio = M2/M1
///
/// Returns:
///
/// * d log(rl) / d log (m2) where rl = Roche lobe radius of
/// secondary star divided by separation according to Eggleton's formula.
///
#[pyfunction]
pub fn zeta_rlobe_eggleton(q: f64) -> Result<f64, RocheError> {
    if q <= 0. {
        let message = format!("q = {} <= 0", q);
        return Err(RocheError::ParameterError(message));
    }
    let q1 = q.powf(1.0 / 3.0);
    let loneq = (1.0 + q1).ln();
    Ok((1. + q) / 3.0 * (2.0 * loneq - q1 / (1. + q1)) / (0.6 * q1 * q1 * loneq))
}

///
/// dzetadq_rlobe_eggleton returns d zeta / d q where zeta is the result
/// of zeta_rlobe_eggleton(double q). This has been tested successfully
/// against finite difference value.
///
/// Arguments:
///
/// * `q`: mass ratio = M2/M1
///
/// Returns:
///
/// * d zeta d q
///
#[pyfunction]
pub fn dzetadq_rlobe_eggleton(q: f64) -> Result<f64, RocheError> {
    if q <= 0. {
        let message = format!("q = {} <= 0", q);
        return Err(RocheError::ParameterError(message));
    }
    let q1 = q.powf(1.0 / 3.0);
    let q2 = q1 * q1;
    let opq1 = 1. + q1;
    let loneq = opq1.ln();
    let denom = 0.6 * q2 + loneq;
    let numer = 2. * loneq - q1 / opq1;
    Ok(numer / denom / 3.0
        + (1.0 + q) / 3.0
            * ((1. + 2. * q1) / 3.0 / (q1 * opq1).powi(2)
                - numer * (0.4 / q1 + 1. / (3. * q2 * (1. + q1))) / denom)
            / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zeta_rlobe_eggleton_test() -> Result<(), RocheError> {
        // Values from trm.roche.zeta_rlobe_eggleton
        assert_eq!(zeta_rlobe_eggleton(0.2)?, 2.3365106916200284);
        assert!(zeta_rlobe_eggleton(-0.2).is_err());
        Ok(())
    }

    #[test]
    fn dzetadq_rlobe_eggleton_test() -> Result<(), RocheError> {
        // Values from trm.roche.zeta_rlobe_eggleton
        assert_eq!(dzetadq_rlobe_eggleton(0.2)?, 0.135112859605613);
        assert!(dzetadq_rlobe_eggleton(-0.2).is_err());
        Ok(())
    }
}
