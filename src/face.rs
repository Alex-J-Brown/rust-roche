use crate::{Vec3, Star};
use crate::{rpot1, rpot2, drpot1, drpot2};
use crate::errors::RocheError;
use pyo3::prelude::*;

/// 
/// 'face' computes the position and orientation of a face on either star in a binary assuming Roche geometry given
/// a direction, a reference radius and a potential.
///
/// \param q    the mass ratio = M2/M1.
/// \param star specifies which star, primary or secondary is under consideration.
/// \param spin ratio of star in questions spin to the orbital frequency
/// \param dirn the direction (unit) vector from the centre of mass of the secondary to the face in question.
/// \param rref reference radius. This is a radius large enough to guarantee crossing of the reference potential. See ref_sphere
/// \param pref reference potential. This defines the precise location of the face.
/// \param acc  location accuracy (units of separation)
/// \return pvec position vector of centre of face (position vector in standard binary coordinates), returned
/// \return dvec orientation vector perpendicular to face, returned
/// \return r    distance from centre of mass of star, returned
/// \return g    magnitude of gravity at face, returned
/// \exception The routine throws exceptions if it cannot bracket the reference potential. This can occur if the reference radius fails to enclose
/// the face in question, or if the face is so deep in the potential that the initial search fails to reach it. Finally if acc is set too low an
/// exception may be thrown if too many binary chops occur. The behaviour at the L1 point is undefined so do not try to call it there.
/// 
#[pyfunction]
pub fn face(q: f64, star: Star, spin: f64, direction: Vec3, rref: f64, pref: f64, acc: f64) -> Result<(Vec3, Vec3, f64, f64), RocheError> {

    let mut pvec: Vec3;
    let mut r: f64;

    let cofm: Vec3 = match star {
        Star::Primary => Vec3::cofm1(),
        Star::Secondary => Vec3::cofm2(),
    };

    let rp: fn(f64, f64, &Vec3) -> Result<f64, RocheError> = match star {
        Star::Primary => rpot1,
        Star::Secondary => rpot2,
    };

    let drp: fn(f64, f64, &Vec3) -> Result<Vec3, RocheError> = match star {
        Star::Primary => drpot1,
        Star::Secondary => drpot2,
    };

    let mut tref: f64 = rp(q, spin, &(cofm + rref*direction))?;
    if tref < pref {
        let message = format!(
            "point at reference radius {} appears to be at lower potential {} than the reference potential {}", rref, tref, pref
        );
        return Err(RocheError::FaceError(message));
    }

    let mut r1: f64 = rref/2.;
    let mut r2: f64 = rref;
    tref = pref + 1.;

    const MAXSEARCH: i32 = 30;
    let mut i: i32 = 0;
    while i < MAXSEARCH && tref > pref {
        r1 = r2/2.;
        tref = rp(q, spin, &(cofm + r1*direction))?;
        if tref > pref {
            r2 = r1;
        }
        i+=1;
    }
    if tref > pref {
        let message = "could not find a radius with a potential below the reference potential; probably bad inputs.";
        return Err(RocheError::FaceError(message.to_string()));
    }

    const MAXCHOP: i32 = 100;
    let mut nchop: i32 = 0;
    while r2 - r1 > acc && nchop < MAXCHOP {
        r = (r1 + r2)/2.;
        pvec = cofm + r*direction;
        if rp(q, spin, &pvec)? < pref {
            r1 = r;
        }else {
            r2 = r;
        }
        nchop += 1;
    }
    if nchop == MAXCHOP {
        return Err(RocheError::FaceError("reached maximum number of binary chops".to_string()));
    }
    r = (r1 + r2)/2.;
    pvec = cofm + r*direction;
    let mut dvec: Vec3 = drp(q, spin, &pvec)?;
    let g = dvec.length();
    dvec /= g;
    Ok((pvec, dvec, r, g))
}