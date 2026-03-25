use std::f64::consts::TAU;
use crate::{Vec3, Star};
use crate::{ref_sphere, sphere_eclipse, set_earth, pot_min, fblink};

/// 
/// ingress_egress tests for whether a given point is eclipsed by a Roche-distorted star. If 
/// it is, it computes the ingress and egress phases using a binary chop. The accuracy on the 
/// phase should be set to be below the expected uncertainties of the phases of your data.
///
/// \param q       the mass ratio = M2/M1.
/// \param star    which star, primary or secondary, is doing the eclipsing
/// \param spin    ratio of spin to orbit of eclipsing star
/// \param ffac    linear filling factor
/// \param iangle  inclination angle
/// \param delta   the accuracy in phase wanted.
/// \param r       position vector of point of interest.
/// \param ingress ingress phase (if eclipsed)
/// \param egress  egress phase 
/// \return false = not eclipsed; true = eclipsed.
/// 
pub fn ingress_egress(q: f64, star: Star, spin: f64, ffac: f64, iangle: f64, delta: f64, r: &Vec3, ingress: &mut f64, egress: &mut f64) -> bool {
    let rref: f64;
    let pref: f64;
    (rref, pref) = ref_sphere(q, star, spin, ffac);
    let ri: f64 = iangle.to_radians();
    let (sini, cosi) = ri.sin_cos();

    let cofm: Vec3 = match star {
        Star::Primary => Vec3::cofm1(),
        Star::Secondary => Vec3::cofm2(),
    };

    let mut phi1: f64 = 0.0;
    let mut phi2: f64 = 0.0;
    let mut lam1: f64 = 0.0;
    let mut lam2: f64 = 0.0;
    let mut phi: f64 = 0.0;
    let mut lam: f64 = 0.0;

    if sphere_eclipse(cosi, sini, r, &cofm, rref, &mut phi1, &mut phi2, &mut lam1, &mut lam2) {
        
        let acc: f64 = 2.*(2.0*TAU*(lam2 - lam1)*delta).sqrt();

        if pot_min(q, star, spin, cosi, sini, r, phi1, phi2, lam1, lam2, rref, pref, acc, &mut phi, &mut lam) {

            let mut pin: f64 = phi;
            let mut pout: f64 = phi1;
            let mut pmid: f64;

            while (pin - pout).abs() > delta {
                pmid = (pin + pout)/2.0;
                if fblink(q, star, spin, ffac, acc, &set_earth(cosi, sini, pmid), r).unwrap() {
                    pin = pmid;
                } else {
                    pout = pmid;
                }
            }
            *ingress = (pin+pout)/2.0;
            *ingress = *ingress - ingress.floor();

            pin = phi;
            pout = phi2;
            while (pin-pout).abs() > delta {
                pmid = (pin+pout)/2.;
                if fblink(q, star, spin, ffac, acc, &set_earth(cosi, sini, pmid), r).unwrap() {
                    pin = pmid;
                } else {
                    pout = pmid;
                }
            }
            *egress = (pin+pout)/2.0;
            *egress = *egress - egress.floor();
            if *egress < *ingress {
                *egress += 1.0;
            }
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }

}