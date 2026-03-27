use crate::Vec3;
use crate::set_earth;
use crate::fblink;
use crate::Star;
use crate::ingress_egress;
use crate::errors::RocheError;
use pyo3::prelude::*;

// q, i, delta_phi make a triad of values, any two of which can be used to find the third.


/// Given a mass ratio and a phase width of the eclipse, findi finds the inclination at which the eclipse has that width. 
/// It returns -1 if the eclipse is always wider than the specified width, -2 if it is always narrower. 
/// Otherwise it returns the inclination in degrees. 
/// \param q the mass ratio = M2/M1
/// \param dphi the phase width of the eclipse
/// \param acc the accuracy to use in fblink.
/// \param delta_i the accuracy in inclination to return. The returned inclination is accurate to within delta_i degrees.
///
#[pyfunction]
#[pyo3(signature = (q, dphi, acc=1.0e-4, delta_i=1.0e-5))]
pub fn findi(q: f64, dphi: f64, acc: f64, delta_i: f64) -> Result<f64, RocheError> {

    if q <= 0.0 {
        return Err(RocheError::ParameterError("mass ratio must be positive".to_string()));
    }
    if dphi <= 0.0 || dphi >= 0.25 {
        return Err(RocheError::ParameterError("phase width must be between 0 and 0.25".to_string()));
    }
    if acc <= 0.0 || acc >= 0.1 {
        return Err(RocheError::ParameterError("accuracy must be between 0 and 0.1".to_string()));
    }
    if delta_i <= 0.0 || delta_i >= 10 {
        return Err(RocheError::ParameterError("delta_i must be between 0 and 10 degrees".to_string()));
    }

    let mut ilo: f64 = 65.0;
    let mut ihi: f64 = 90.0;
    let phi: f64 = dphi/2.0;
    let earth1: Vec3 = set_earth::set_earth_iangle(ilo, phi);
    let earth2: Vec3 = set_earth::set_earth_iangle(ihi, phi);
    let r: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    
    // eclipsed at ilo?
    let elo: bool = match fblink::fblink(q, Star::Secondary, 1.0, 1.0, acc, &earth1, &r) {
        Ok(result) => result,
        Err(_) => return Err(RocheError::DbrentError),
    };
    let ehi: bool = match fblink::fblink(q, Star::Secondary, 1.0, 1.0, acc, &earth2, &r) {
        Ok(result) => result,
        Err(_) => return Err(RocheError::DbrentError),
    };
    if elo && ehi {
        return Ok(-2.0);
    } else if !elo && !ehi {
        return Ok(-1.0);
    }
    while (ihi - ilo) > delta_i {
        let imid: f64 = (ilo + ihi)/2.0;
        let earth_mid: Vec3 = set_earth::set_earth_iangle(imid, phi);
        let emid: bool = match fblink::fblink(q, Star::Secondary, 1.0, 1.0, acc, &earth_mid, &r) {
            Ok(result) => result,
            Err(_) => return Err(RocheError::DbrentError),
        };
        if emid {
            ihi = imid;
        } else {
            ilo = imid;
        }
    }
    Ok((ilo + ihi)/2.0)
}


/// Given an inclination and a phase width of the eclipse, findq finds the mass ratio at which the eclipse has that width. 
/// It returns -1 if the eclipse is always wider than the specified width, -2 if it is always narrower. 
/// Otherwise it returns the mass ratio. 
/// \param i the inclination in degrees
/// \param dphi the phase width of the eclipse
/// \param acc the accuracy to use in fblink.
/// \param delta_q the accuracy in mass ratio to return. The returned mass ratio is accurate to within delta_q.
///
#[pyfunction]
#[pyo3(signature = (i, dphi, acc=1.0e-4, delta_q=1.0e-5))]
pub fn findq(i: f64, dphi: f64, acc: f64, delta_q: f64) -> Result<f64, RocheError> {

    if i <= 0.0 || i >= 90.0 {
        return Err(RocheError::ParameterError("inclination must be between 0 and 90 degrees".to_string()));
    }
    if dphi <= 0.0 || dphi >= 0.25 {
        return Err(RocheError::ParameterError("phase width must be between 0 and 0.25".to_string()));
    }
    if acc <= 0.0 || acc >= 0.1 {
        return Err(RocheError::ParameterError("accuracy must be between 0 and 0.1".to_string()));
    }
    if delta_q <= 0.0 || delta_q >= 0.1 {
        return Err(RocheError::ParameterError("delta_q must be between 0 and 0.1".to_string()));
    }

    let mut qlo: f64 = 0.001;
    let mut qhi: f64 = 2.0;
    let phi: f64 = dphi/2.0;
    let earth: Vec3 = set_earth::set_earth_iangle(i, phi);
    let r: Vec3 = Vec3::new(0.0, 0.0, 0.0);

    let elo: bool = match fblink::fblink(qlo, Star::Secondary, 1.0, 1.0, acc, &earth, &r) {
        Ok(result) => result,
        Err(_) => return Err(RocheError::DbrentError),
    };
    let ehi: bool = match fblink::fblink(qhi, Star::Secondary, 1.0, 1.0, acc, &earth, &r) {
        Ok(result) => result,
        Err(_) => return Err(RocheError::DbrentError),
    };
    if elo && ehi {
        return Ok(-2.0);
    } else if !elo && !ehi {
        return Ok(-1.0);
    }
    while (qhi - qlo) > delta_q {
        let qmid: f64 = (qlo + qhi)/2.0;
        let emid: bool = match fblink::fblink(qmid, Star::Secondary, 1.0, 1.0, acc, &earth, &r) {
            Ok(result) => result,
            Err(_) => return Err(RocheError::DbrentError),
        };
        if emid {
            qhi = qmid;
        } else {
            qlo = qmid;
        }
    }
    Ok((qlo + qhi)/2.0)
}

#[pyfunction]
#[pyo3(signature = (q, iangle, delta=1.0e-6))]
/// Given an mass ration and inclination, findphi finds the phase width of the white dwarf eclipse. 
/// It returns -1 if white dwarf is not eclipsed.
///
/// \param q the mass ratio = M2/M1
/// \param i the inclination in degrees
/// \param acc the accuracy to use in fblink.
/// \param delta the accuracy to use in ingress_egress.
/// \return the phase width of the white dwarf eclipse, or -1 if the white dwarf is not eclipsed.
///
pub fn findphi(q: f64, iangle: f64, delta: f64) -> f64 {

    let r: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    let mut ingress: f64 = 0.0;
    let mut egress: f64 = 0.0;
    //q: f64, star: Star, spin: f64, ffac: f64, iangle: f64, delta: f64, r: &Vec3, ingress: &mut f64, egress: &mut f64
    let status: bool = ingress_egress::ingress_egress(
        q, Star::Secondary, 1.0, 1.0, iangle, delta, &r, &mut ingress, &mut egress
    );
    if !status {
        return -1.0;
    }
    egress-ingress
}