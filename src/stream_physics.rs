use bulirsch::{self, Integrator};
use crate::Vec3;
use crate::errors::RocheError;
use crate::x_l1;
use pyo3::prelude::*;

/// 
/// strinit sets a particle just inside the L1 point with the 
/// correct velocity as given in Lubow and Shu.
///
/// Arguments:
/// 
/// * `q`: mass ratio = M2/M1
/// 
/// Returns:
/// 
/// * start position
/// * start velocity
///
#[pyfunction]
pub fn strinit(q: f64) -> Result<(Vec3, Vec3), RocheError> {
    
    const SMALL: f64 = 1.0e-5;
    let rl1: f64 = x_l1(q)?;
    let mu: f64 = q/(1.0+q);
    let a: f64 = (1.0-mu)/rl1.powi(3)+mu/(1.0-rl1).powi(3);
    let lambda1: f64 = (((a-2.0) + (a*(9.0*a-8.0)).sqrt())/2.0).sqrt();
    let m1: f64 = (lambda1*lambda1-2.0*a-1.0)/2.0/lambda1;

    let r: Vec3 = Vec3::new(rl1-SMALL, -m1*SMALL, 0.0);
    let v: Vec3 = Vec3::new(-lambda1*SMALL, -lambda1*m1*SMALL, 0.0);

    Ok((r, v))

}


/// 
/// stradv advances a particle of given position and velocity until
/// it reaches a specified radius. It then returns with updated position and
/// velocity. It is up to the user not to request a value that cannot be reached.
///
/// Arguments:
/// 
/// * `q`:    mass ratio = M2/M1
/// * `r`:    Initial and final position
/// * `v`:    Initial and final velocity
/// * `rad`:  Radius to aim for
/// * `acc`:  Accuracy with which to place output point at rad.
/// * `smax`: Largest time step allowed. It is possible that the
/// routine could take such a large step that it misses
/// the point when the stream is inside the requested
/// radius. This allows one to control this. Typical
/// value = 1.e-3.
///
/// Returns:
/// 
/// * time step taken
///
pub fn stradv(q: f64, r: &mut Vec3, v: &mut Vec3, rad: f64, acc: f64, smax: f64) -> f64 {

    const TMAX: f64 = 10.0;
    let t_next: f64 = 1.0e-2;

    let mut time: f64 = 0.0;

    // let to: f64;
    let mut ro = *r;
    let mut vo = *v;
    
    // Store initial radius
    let rinit: f64 = r.length();
    let mut rnow: f64 = rinit;

    // set up Bulirsch-Stoer integrator
    let system = OrbitalSystem{ q: q };
    let mut integrator = Integrator::default().with_abs_tol(1.0e-8).with_rel_tol(1.0e-8).into_adaptive();
    // Initialise arrays
    let mut y = ndarray::array![r.x, r.y, r.z, v.x, v.y, v.z];
    let mut y_next = ndarray::Array::zeros(y.raw_dim());
    
    let mut yo = y.clone();
    let mut delta_t = t_next.min(smax);
    // Step until radius crossed
    while (rinit > rad && rnow > rad) || (rinit < rad && rnow < rad) {
        ro = *r;
        vo = *v;
        yo = y.clone();
        integrator
            .step(&system, delta_t, y.view(), y_next.view_mut())
            .unwrap();
        y.assign(&y_next);
        r.set(y[0], y[1], y[2]);
        v.set(y[3], y[4], y[5]);
        rnow = r.length();
        time += delta_t;
        
        if time > TMAX {
            panic!("roche::stradv taken too long without crossing given radius.")
        }
    }

    // Now refine by reinitialising and binary chopping until
    // close enough to requested radius.

    let mut lo: f64 = 0.0;
    let mut hi: f64 = delta_t;
    let mut rlo: f64 = ro.length();
    let mut rhi: f64 = rnow;
    let to: f64 = time;

    while (rhi-rlo).abs() > acc {
        delta_t = (lo+hi)/2.0;
        y = yo.clone();
        *r = ro;
        *v = vo;
        time = to;

        integrator
            .step(&system, delta_t, y.view(), y_next.view_mut())
            .unwrap();
        y.assign(&y_next);

        r.set(y[0], y[1], y[2]);
        v.set(y[3], y[4], y[5]);
        rnow = r.length();

        if (rhi > rad && rnow > rad) || (rhi < rad && rnow < rad) {
            rhi = rnow;
            hi = delta_t;
        } else {
            rlo = rnow;
            lo = delta_t;
        }
    }

    time

}

// wrapper for python library, avoiding mutable references

/// 
/// stradv advances a particle of given position and velocity until
/// it reaches a specified radius. It then returns with updated position and
/// velocity. It is up to the user not to request a value that cannot be reached.
///
/// \param q    mass ratio = M2/M1
/// \param r    Initial position
/// \param v    Initial velocity
/// \param rad  Radius to aim for
/// \param acc  Accuracy with which to place output point at rad.
/// \param smax Largest time step allowed. It is possible that the
/// routine could take such a large step that it misses
/// the point when the stream is inside the requested
/// radius. This allows one to control this. Typical
/// value = 1.e-3.
/// \returns (timestep, new position, new velocity)
///
#[pyfunction]
#[pyo3(name = "stradv")]
pub fn stradv_wrapper(q: f64, r: &Vec3, v: &Vec3, rad: f64, acc: f64, smax: f64) -> (f64, Vec3, Vec3) {
    let mut r_mut = *r;
    let mut v_mut = *v;
    let timestep = stradv(q, &mut r_mut, &mut v_mut, rad, acc, smax);
    (timestep, r_mut, v_mut)
}

///
/// rocacc calculates and returns the acceleration (in the rotating frame)
/// in a Roche potential of a particle of given position and velocity.
///
/// \param q mass ratio = M2/M1
/// \param r position, scaled in units of separation.
/// \param v velocity, scaled in units of separation
///
#[pyfunction]
pub fn rocacc(q: f64, r: &Vec3, v: &Vec3) -> (f64, f64, f64) {


    let f1: f64 = 1.0 / (1.0+q);
    let f2: f64 = f1*q;

    let yzsq: f64 = r.y*r.y + r.z*r.z;
    let r1sq: f64 = r.x*r.x + yzsq;
    let r2sq: f64 = (r.x-1.0)*(r.x-1.0) + yzsq;
    let fm1: f64 = f1/(r1sq*(r1sq.sqrt()));
    let fm2: f64 = f2/(r2sq*(r2sq.sqrt()));
    let fm3 = fm1+fm2;

    let x: f64 = -fm3*r.x + fm2 + 2.0*v.y + r.x - f2;
    let y: f64 = -fm3*r.y       - 2.0*v.x + r.y;
    let z: f64 = -fm3*r.z;
    (x, y, z)
}


struct OrbitalSystem {
    q: f64,
}

impl bulirsch::System for OrbitalSystem {
    type Float = f64;
    
    fn system(&self, y: bulirsch::ArrayView1<Self::Float>, mut dydt: bulirsch::ArrayViewMut1<Self::Float>) {
        dydt[[0]] = y[[3]];
        dydt[[1]] = y[[4]];
        dydt[[2]] = y[[5]];
        let r = Vec3::new(y[[0]], y[[1]], y[[2]]);
        let v = Vec3::new(y[[3]], y[[4]], y[[5]]);
        (dydt[[3]], dydt[[4]], dydt[[5]]) = rocacc(self.q, &r, &v);
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strinit_stradv_test() -> Result<(), RocheError> {
        // Values from trm.roche.bspot
        let (mut r, mut v) = strinit(0.2)?;
        let _time = stradv(0.2, &mut r, &mut v, 0.3, 1.0e-7, 1.0e-3);
        assert!((r - Vec3::new(0.2660591412807423, 0.13860932478255575, 0.0)).length() < 1.0e-7);
        assert!((v - Vec3::new(-1.4769457229627583, 0.31712381217252994, 0.0)).length() < 1.0e-7);

        Ok(())
    }
}