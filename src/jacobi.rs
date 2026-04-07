use crate::Vec3;
use pyo3::prelude::*;

///
/// returns the Jacobi constant corresponding to a particular
/// position and velocity.
///
/// Arguments:
///
/// * `q`: mass ratio M2/M1
/// * `r`: vector coordinates of position
/// * `v`: velocity vector
///
/// Returns:
///
/// Jacobi constant
///
#[pyfunction]
pub fn jacobi(q: f64, r: &Vec3, v: &Vec3) -> f64 {
    let f1: f64 = 1.0 / (1.0 + q);
    let f2: f64 = f1 * q;

    let yz_sqr: f64 = r.y.powi(2) + r.z.powi(2);

    (v.x.powi(2) + v.y.powi(2) + v.z.powi(2) - r.y.powi(2) - (r.x - f2).powi(2)) / 2.0
        - f1 / (r.x.powi(2) + yz_sqr).sqrt()
        - f2 / ((r.x - 1.0).powi(2) + yz_sqr).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_test() -> () {
        // Values from trm.roche.jacobi
        let r = Vec3::new(0.2, 0.3, 0.05);
        let v = Vec3::new(0.1, 0.1, 0.0);
        assert_eq!(jacobi(0.2, &r, &v), -2.5196337006460348);
    }
}
