use crate::errors::RocheError;
use crate::{Star, Vec3};
use crate::{rpot1, rpot2, x_l1_1, x_l1_2};
use pyo3::prelude::*;

///
/// ref_sphere computes the radius of a reference sphere just outside a Roche distorted star
/// along the line of centres and centred upon its centre of mass. This sphere, which is guaranteed
/// to enclose the equipotential in question can then be used to define regions for searching for
/// equipotential crossing when computing eclipses. The Roche-distorted star is defined by the mass
/// ratio and the (linear) filling factor defined as the distance from the centre of mass of the
/// star to its surface in the direction of the L1 point divided by the distance to the L1 point.
/// A filling factor = 1 is Roche filling. Note that the size of the Roche lobe is calculated as
/// appropriate given any asynchronism.
///
/// Arguments:
///
/// * `q`: the mass ratio = M2/M1
/// * `star`: specifies which star, primary or secondary is under consideration
/// * `spin`: ratio spin/orbital frequencies to allow for asynchronism
/// * `ffac`: linear filling factor.
///
/// Returns:
///
/// * the radius of the reference sphere. This will be 0.1% expanded above the minimum
/// size to avoid round off bugs, if it remains within Roche lobe.
/// * the reference potential. Roche potential on surface of distorted star.
///
#[pyfunction]
pub fn ref_sphere(q: f64, star: Star, spin: f64, ffac: f64) -> Result<(f64, f64), RocheError> {
    let tref: f64;
    let rref: f64;
    let pref: f64;

    if star == Star::Primary {
        tref = x_l1_1(q, spin)?;
        rref = tref * 1.0_f64.min(1.001 * ffac);
        pref = rpot1(
            q,
            spin,
            &Vec3 {
                x: ffac * tref,
                y: 0.0,
                z: 0.0,
            },
        )?;
        Ok((rref, pref))
    } else if star == Star::Secondary {
        tref = 1.0 - x_l1_2(q, spin)?;
        rref = tref * 1.0_f64.min(1.001 * ffac);
        pref = rpot2(
            q,
            spin,
            &Vec3 {
                x: 1.0 - ffac * tref,
                y: 0.0,
                z: 0.0,
            },
        )?;
        Ok((rref, pref))
    } else {
        let message = format!("{:?} is not and instance of Star.", star);
        return Err(RocheError::ParameterError(message));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ref_sphere_test() -> Result<(), RocheError> {
        // Values from trm.roche.ref_sphere
        assert_eq!(
            ref_sphere(0.2, Star::Secondary, 1.0, 0.8)?,
            (0.27342861229381593, -2.3996722470168605)
        );
        assert!(ref_sphere(-0.2, Star::Secondary, 1.0, 0.8).is_err());
        Ok(())
    }
}
