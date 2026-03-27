use std::error::Error;
use crate::{Vec3, Star};
use crate::{rpot_val, rpot_grad, ref_sphere, sphere_eclipse_vector, dbrent};
use pyo3::prelude::*;
use pyo3::exceptions::PyOSError;
use std::fmt;

#[derive(Debug)]
pub struct FblinkError;

impl std::error::Error for FblinkError {}

impl fmt::Display for FblinkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "failed to bracket minimum with dbrent")
    }
}

impl std::convert::From<FblinkError> for PyErr {
    fn from(err: FblinkError) -> PyErr {
        PyOSError::new_err(err.to_string())
    }
}

/// 
/// fblink works out whether or not a given point is eclipsed by a Roche-distorted star by searching along the line
/// of sight to the point to see if the Roche potential ever drops below the value at the stellar surface.
///
/// \param q      mass ratio = M2/M1
/// \param star   star concerned
/// \param spin   ratio of spin to orbital frequency
/// \param star   which star is doing the eclipsing, primary or secondary
/// \param ffac   the filling factor of the star
/// \param acc    accuracy of location of minimum potential, units of separation. The accuracy in height relative to the
/// Roche potential is acc*acc/(2*R) where R is the radius of curvature of the Roche potential surface, so don't be too
/// picky. 1.e-4 would be more than good enough in most cases.
/// \param earth  vector pointing towards earth
/// \param p      point of interest
/// \return true if minimum potential is below the potential at stellar surface
///
#[pyfunction]
pub fn fblink(q: f64, star: Star, spin: f64, ffac: f64, acc: f64, earth: &Vec3, p: &Vec3) -> Result<bool, FblinkError> {

    let (rref, pref) = ref_sphere(q, star, spin, ffac);

    let cofm: Vec3 = match star {
        Star::Primary => Vec3::cofm1(),
        Star::Secondary => Vec3::cofm2(),
    };
    
    // First compute the multipliers cutting the reference sphere (if any)
    let mut lam1 = 0.0;
    let mut lam2 = 0.0;
    if !sphere_eclipse_vector(earth, p, &cofm, rref, &mut lam1, &mut lam2) {
        return Ok(false);
    }
    if lam1 == 0.0 {
        return Ok(true);
    }

    // Create function objects for 1D minimisation in lambda direction
    let func = |lam: f64| {
        rpot_val(q, star, spin, earth, p, lam)
    };

    // Now try to bracket a minimum. We just crudely compute function at regularly spaced intervals filling in the
    // gaps until the step size between the points drops below the threshold. Take every opportunity to jump out early
    // either if the potential is below the threshold or if we have bracketed a minimum.
    let mut nstep: i32 = 1;
    let mut step: f64 = lam2 - lam1;

    let mut f1: f64 = 0.0;
    let mut f2: f64 = 0.0;
    let mut flam: f64 = 1.0;
    let mut lam: f64 = lam1;

    while step > acc {

        lam = lam1 + step/2.0;

        for _ in 0..nstep {

            flam = func(lam);
            if flam <= pref {
                return Ok(true);
            }

            // Calculate these as late as possible because they may often not be needed
            if nstep == 1 {
                f1 = func(lam1);
                f2 = func(lam2);
            }

            if flam < f1 && flam < f2 {
                break;
            }

            lam += step;
        }
        if flam < f1 && flam < f2 {
            break;
        }
        step /= 2.0;
        nstep *= 2;
    }

    if flam < f1 && flam < f2 {

        // OK, minimum bracketted, so finally pin it down accurately
        // Possible that multiple minima could cause problems but I have
        // never seen this in practice.
        let dfunc = |lam: f64| {
            let (_dp, dl) = rpot_grad(q, star, spin, earth, p, lam);
            dl
        };

        // this may fail
        let (_xmin, flam) =  match dbrent(lam1, lam, lam2, |x| func(x), |x| dfunc(x), acc, true, pref){
            Ok(res) => res,
            Err(_) =>  return Err(FblinkError),
        };

        Ok(flam < pref)
    } else {
        // Not bracketted even after a detailed search, and we have not jumped 
        // out either, so assume no eclipse
        Ok(false)
    }
}