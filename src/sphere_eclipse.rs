use crate::Vec3;
use pyo3::prelude::*;
use std::f64::consts::TAU;

///
/// sphere_eclipse tells you whether or not a given sphere will eclipse a
/// given point or not. If the answer is yes, it will also return with four
/// parameters to define the phase range and the multiplier range delimiting the
/// region within which the spheres surface is crossed.
/// These can then be used as the starting point for later computation.
/// (The line of sight is described as the point in question plus a scalar multiplier
/// times a unit vector pointing towards Earth -- this is the "multiplier" referred to above
/// and below). The multiplier must be positive: in other words the routine does not
/// project backwards. If the point in inside the sphere, phi1 will be set = 0, phi2 = 1,
/// lam1 = 0, and lam2 = the largest value of the multiplier lambda
///
/// Arguments:
///
/// * `cosi`: cosine of orbital inclination
/// * `sini`: sine of orbital inclination
/// * `r`: the position vector of the point in question (units of binary separation)
/// * `c`: the centre of the sphere enclosing the star (units of binary separation)
/// * `rsphere`: the radius defining the sphere enclosing the star (units of binary separation)
/// * `phi1`: if eclipse by the sphere, this is the start of the phase range (0 -1)
/// * `phi2`: if eclipse by the sphere, this is the end of the phase range (least amount > phi1)
/// * `lam1`: if eclipse by the sphere, this is the start of the multiplier range (>=0)
/// * `lam2`: if eclipse by the sphere, this is the end of the multiplier range (>lam1)
///
/// Returns:
///
/// * false = not eclipsed; true = eclipsed.
///
pub fn sphere_eclipse(
    cosi: f64,
    sini: f64,
    r: &Vec3,
    c: &Vec3,
    rsphere: f64,
    phi1: &mut f64,
    phi2: &mut f64,
    lam1: &mut f64,
    lam2: &mut f64,
) -> bool {
    let d: Vec3 = *r - *c;

    let pdist: f64 = (d.x * d.x + d.y * d.y).sqrt();
    let bquad: f64 = d.z * cosi - pdist * sini;
    if bquad >= 0. {
        return false;
    }
    let cquad: f64 = d.sqr() - rsphere * rsphere;

    let mut fac: f64 = bquad * bquad - cquad;
    if fac <= 0.0 {
        return false;
    }
    fac = fac.sqrt();

    *lam2 = -bquad + fac;
    *lam1 = 0_f64.max(cquad / (*lam2));

    if cquad < 0. {
        *phi1 = 0.;
        *phi2 = 1.
    } else {
        let delta: f64 = ((cosi * d.z + cquad.sqrt()) / (sini * pdist)).acos();
        let phi: f64 = d.y.atan2(-d.x);
        *phi1 = (phi - delta) / TAU;
        *phi1 -= phi1.floor();
        *phi2 = *phi1 + 2. * delta / TAU;
    }
    return true;
}

// wrapper of the above function to use with Python (e.g avoiding mutable arguments)

///
/// sphere_eclipse tells you whether or not a given sphere will eclipse a
/// given point or not. If the answer is yes, it will also return with four
/// parameters to define the phase range and the multiplier range delimiting the
/// region within which the spheres surface is crossed.
/// These can then be used as the starting point for later computation.
/// (The line of sight is described as the point in question plus a scalar multiplier
/// times a unit vector pointing towards Earth -- this is the "multiplier" referred to above
/// and below). The multiplier must be positive: in other words the routine does not
/// project backwards. If the point in inside the sphere, phi1 will be set = 0, phi2 = 1,
/// lam1 = 0, and lam2 = the largest value of the multiplier lambda
///
/// Arguments:
///
/// * `cosi`:     cosine of orbital inclination
/// * `sini`:     sine of orbital inclination
/// * `r (Vec3)`: the position vector of the point in question (units of binary separation)
/// * `c (Vec3)`: the centre of the sphere enclosing the star (units of binary separation)
/// * `rsphere`:  the radius defining the sphere enclosing the star (units of binary separation)
///
/// Returns:
///
///  * (eclipsed, phi1, phi2, lam1, lam2)
///
#[pyfunction]
#[pyo3(name = "sphere_eclipse")]
pub fn sphere_eclipse_wrapper(
    cosi: f64,
    sini: f64,
    r: &Vec3,
    c: &Vec3,
    rsphere: f64,
) -> (bool, f64, f64, f64, f64) {
    let mut phi1 = 0.0;
    let mut phi2 = 0.0;
    let mut lam1 = 0.0;
    let mut lam2 = 0.0;

    let eclipsed = sphere_eclipse(
        cosi, sini, r, c, rsphere, &mut phi1, &mut phi2, &mut lam1, &mut lam2,
    );

    (eclipsed, phi1, phi2, lam1, lam2)
}

///
/// This version of sphere_eclipse tells you whether or not a given sphere will eclipse
/// a given point at a particular phase or not. If the answer is yes,
/// it will also return with the multiplier values giving the cut points. These can then
/// be used as starting points for Roche lobe computations. These can then be used as the
/// starting point for later computation. Points inside the sphere are regarded as being
/// eclipsed with the lower mulitplier set = 0
///
/// Arguments:
///
/// * `earth`:   vector towards Earth
/// * `r`:       the position vector of the point in question (units of binary separation)
/// * `c`:       the centre of the sphere (units of binary separation)
/// * `rsphere`: the radius defining the sphere enclosing the star (units of binary separation)
/// * `lam1`:    if eclipse by the sphere, this is the start of the multiplier range
/// * `lam2`:    if eclipse by the sphere, this is the end of the multiplier range
///
/// Returns:
///
/// * false = not eclipsed; true = eclipsed.
///
pub fn sphere_eclipse_vector(
    earth: &Vec3,
    r: &Vec3,
    c: &Vec3,
    rsphere: f64,
    lam1: &mut f64,
    lam2: &mut f64,
) -> bool {
    let d: Vec3 = *r - *c;

    let bquad: f64 = earth.dot(&d);
    if bquad >= 0. {
        return false;
    }
    let cquad: f64 = d.sqr() - rsphere * rsphere;

    let mut fac: f64 = bquad * bquad - cquad;
    if fac <= 0. {
        return false;
    }
    fac = fac.sqrt();

    *lam2 = -bquad + fac;
    *lam1 = 0_f64.max(cquad / (*lam2));
    return true;
}

// wrapper of the above function to use with Python (e.g avoiding mutable arguments)

///
/// This version of sphere_eclipse tells you whether or not a given sphere will eclipse
/// a given point at a particular phase or not. If the answer is yes,
/// it will also return with the multiplier values giving the cut points. These can then
/// be used as starting points for Roche lobe computations. These can then be used as the
/// starting point for later computation. Points inside the sphere are regarded as being
/// eclipsed with the lower multiplier set = 0.
/// The multiplier along the line of sight is lambda, with the smallest and largest values
/// returned as lam1 and lam2, respectively.
///
/// Arguments:
///
/// * `earth`:   vector towards Earth
/// * `r`:       the position vector of the point in question (units of binary separation)
/// * `c`:       the centre of the sphere (units of binary separation)
/// * `rsphere`: the radius defining the sphere enclosing the star (units of binary separation)
///
/// Returns:
///
/// * (eclipsed, lam1, lam2)
///
#[pyfunction]
#[pyo3(name = "sphere_eclipse_vector")]
pub fn sphere_eclipse_vector_wrapper(
    earth: &Vec3,
    r: &Vec3,
    c: &Vec3,
    rsphere: f64,
) -> (bool, f64, f64) {
    let mut lam1 = 0.0;
    let mut lam2 = 0.0;

    let eclipsed = sphere_eclipse_vector(earth, r, c, rsphere, &mut lam1, &mut lam2);

    (eclipsed, lam1, lam2)
}
