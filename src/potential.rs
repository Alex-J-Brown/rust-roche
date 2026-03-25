use std::f64::consts::PI;
use crate::{Vec3, Star};


pub fn rpot_val(q: f64, star: Star, spin: f64, earth: &Vec3, p: &Vec3, lam: f64) -> f64 {

    let r: Vec3 = *p + lam* *earth;
    match star {
        Star::Primary => rpot1(q, spin, &r),
        Star::Secondary => rpot2(q, spin, &r),
        }
}


pub fn rpot_val_grad(q: f64, star: Star, spin: f64, earth: &Vec3, p: &Vec3, lam: f64) -> (f64, f64, f64) {

        let r: Vec3 = *p + lam* *earth;
        let d: Vec3 = match star {
            Star::Primary => drpot1(q, spin, &r),
            Star::Secondary => drpot2(q, spin, &r),
        };
        let rpot: f64 = match star {
            Star::Primary => rpot1(q, spin, &r),
            Star::Secondary => rpot2(q, spin, &r),
        };
        let ed: Vec3 = Vec3::new(earth.y, -earth.x, 0.);
        let dphi: f64 = 2.*PI*lam*d.dot(&ed);
        let dlam: f64 = d.dot(earth);

        (rpot, dphi, dlam)
    }


pub fn rpot_grad(
    q: f64,
    star: Star,
    spin: f64,
    earth: &Vec3,
    p: &Vec3,
    lam: f64
) -> (f64, f64) {

        let r: Vec3 = *p + lam * *earth;

        let d: Vec3 = match star {
            Star::Primary => drpot1(q, spin, &r),
            Star::Secondary => drpot2(q, spin, &r),
        };

        let ed: Vec3 = Vec3::new(earth.y, -earth.x, 0.);

        let dphi: f64 = 2.*PI * lam*d.dot(&ed);
        let dlam: f64 = d.dot(earth);

        (dphi, dlam)
}


///
/// rpot computes the Roche potential at a given point. This is for the standard synchronised Roche geometry
/// \param q mass ratio = M2/M1
/// \param p the point in question (units scaled by separation)
/// \return the Roche potential.
///
pub fn rpot(q: f64, p: &Vec3) -> f64 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let mu: f64 = q/(1. + q);
    let comp: f64 = 1. - mu;
    let x2y2: f64 = p.x*p.x + p.y*p.y;
    let z2: f64 = p.z*p.z;
    let r1sq: f64 = x2y2+z2;
    let r1: f64 = r1sq.sqrt();
    let r2: f64 = (r1sq + 1. - 2.*p.x).sqrt();
    -comp/r1 - mu/r2 - (x2y2 + mu*(mu - 2.*p.x))/2.
}


///
/// rpot1 computes the Roche potential at a given point allowing for non-synchronous rotation of the primary.
/// \param q mass ratio = M2/M1
/// \paran spin ratio spin/orbital frequencies
/// \param p the point in question (units scaled by separation)
/// \return the Roche potential.
///
pub fn rpot1(q: f64, spin: f64, p: &Vec3) -> f64 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let mu: f64 = q/(1. + q);
    let comp: f64 = 1. - mu;
    let x2y2: f64 = p.x*p.x + p.y*p.y;
    let z2: f64 = p.z*p.z;
    let r1sq: f64 = x2y2+z2;
    let r1: f64 = r1sq.sqrt();
    let r2: f64 = (r1sq + 1. - 2.*p.x).sqrt();
    -comp/r1 - mu/r2 - spin*spin*x2y2/2. + mu*p.x
}


///
/// rpot2 computes the Roche potential at a given point allowing for non-synchronous rotation of the secondary.
/// \param q mass ratio = M2/M1
/// \paran spin ratio spin/orbital frequencies
/// \param p the point in question (units scaled by separation)
/// \return the Roche potential.
///
pub fn rpot2(q: f64, spin: f64, p: &Vec3) -> f64 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let mu: f64 = q/(1. + q);
    let comp: f64 = 1. - mu;
    let x2y2: f64 = p.x*p.x + p.y*p.y;
    let z2: f64 = p.z*p.z;
    let r1sq: f64 = x2y2+z2;
    let r1: f64 = r1sq.sqrt();
    let r2: f64 = (r1sq + 1. - 2.*p.x).sqrt();
    -comp/r1 - mu/r2 - spin*spin*(0.5 + 0.5*x2y2 - p.x) - comp*p.x
}


///
/// drpot computes partial derivatives of Roche potential with respect
/// to position at p for mass ratio q.
/// \param q mass ratio  = M2/M1
/// \param p position (units of separation)
///
pub fn drpot(q: f64, p: &Vec3) -> Vec3 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let r1sq: f64 = p.sqr();
    let r1: f64 = r1sq.sqrt();
    let r2sq: f64 = r1sq + 1. - 2.*p.x;
    let r2: f64 = r2sq.sqrt();
    let mu: f64 = q/(1. + q);
    let mu1: f64 = mu/r2/r2sq;
    let comp: f64 = (1. - mu)/r1/r1sq;
    Vec3::new(comp*p.x + mu1*(p.x - 1.) - p.x + mu,
              comp*p.y + mu1*p.y - p.y,
              comp*p.z + mu1*p.z)
}


///
/// drpot computes partial derivatives of Roche potential with respect
/// to position at p for mass ratio q.
/// \param q mass ratio  = M2/M1
/// \param p position (units of separation)
///
pub fn drpot1(q: f64, spin: f64, p: &Vec3) -> Vec3 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let r1sq: f64 = p.sqr();
    let r1: f64 = r1sq.sqrt();
    let r2sq: f64 = r1sq + 1. - 2.*p.x;
    let r2: f64 = r2sq.sqrt();
    let mu: f64 = q/(1. + q);
    let mu1: f64 = mu/r2/r2sq;
    let comp: f64 = (1. - mu)/r1/r1sq;
    let ssq: f64 = spin*spin;
    Vec3::new(comp*p.x + mu1*(p.x - 1.) - ssq*p.x + mu,
              comp*p.y + mu1*p.y - ssq*p.y,
              comp*p.z + mu1*p.z)
}


///
/// drpot2 computes partial derivatives of asynchronous Roche potential with respect
/// to position at p for mass ratio q, star 2.
/// \param q mass ratio  = M2/M1
/// \param spin ratio of spin to orbital frequency.
/// \param p position (units of separation)
///
pub fn drpot2(q: f64, spin: f64, p: &Vec3) -> Vec3 {
    if q <= 0. {
        panic!("q = {} <= 0", q);
    }
    let r1sq: f64 = p.sqr();
    let r1: f64 = r1sq.sqrt();
    let r2sq: f64 = r1sq + 1. - 2.*p.x;
    let r2: f64 = r2sq.sqrt();
    let mu: f64 = q/(1. + q);
    let mu1: f64 = mu/r2/r2sq;
    let comp: f64 = (1. - mu)/r1/r1sq;
    let ssq: f64 = spin*spin;
    Vec3::new(comp*p.x + mu1*(p.x - 1.) - ssq*(p.x - 1.) + mu - 1.,
              comp*p.y + mu1*p.y - ssq*p.y,
              comp*p.z + mu1*p.z)
}
