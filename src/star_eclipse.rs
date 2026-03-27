use crate::{Vec3, Star, Etype};
use crate::{ingress_egress, sphere_eclipse};

pub fn star_eclipse(q: f64, spin:f64, r: f64, ffac: f64, iangle: f64, posn: &Vec3, delta: f64, roche: bool, star: Star, eclipses: &mut Etype) -> () {
    let ri = iangle.to_radians();
    let (sini, cosi) = ri.sin_cos();
    let cofm = match star {
        Star::Primary => Vec3::cofm1(),
        Star::Secondary => Vec3::cofm2(),
    };
    let mut lam1: f64 = 0.0;
    let mut lam2: f64 = 0.0;
    let mut ingress: f64 = 0.0;
    let mut egress: f64 = 0.0;
    // let mut eclipses = Etype::new();
    if (roche && ingress_egress(q, star, spin, ffac, iangle, delta, &posn, &mut ingress, &mut egress)) ||
        (!roche && sphere_eclipse(cosi, sini, &posn, &cofm, r, &mut ingress, &mut egress, &mut lam1, &mut lam2)) {
        eclipses.push((ingress, egress));
        }
}