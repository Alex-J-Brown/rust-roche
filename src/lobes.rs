use crate::{Vec3, Star};
use crate::x_lagrange::x_l1;
use crate::errors::RocheError;
use crate::potential::{rpot, drpot};
use pyo3::prelude::*;

// structure to find specific roche potential along a line
struct LineRoche{
    // mass ratio
    q: f64,
    // which star are we concerned with? (primary or secondary)
    star: Star,
    // direction of line in x
    dx: f64,
    // direction of line in y
    dy: f64,
    // critical potential to solve for
    cpot: f64
    
}

impl LineRoche {
    pub fn new(q: f64, star: Star, dx: f64, dy: f64, cpot: f64) -> Self {
        Self { q, star, dx, dy, cpot }
    }

    fn cost(&self, lam: f64) -> Result<(f64, f64), RocheError>{
        let p = match self.star {
            Star::Primary => Vec3::new(lam*self.dx, lam*self.dy, 0.),
            Star::Secondary => Vec3::new(1.0+lam*self.dx, lam*self.dy, 0.),
        };
        // how far are we from root?
        let f =  rpot(self.q, &p)? - self.cpot;
        // gradient of potential at point
        let dp = drpot(self.q, &p)?;
        // dot product of gradient with line direction - gives scalar gradient in direction of line
        let d = self.dx*dp.x + self.dy*dp.y;
        Ok((f, d))
    }
}

#[pyfunction]
pub fn lobe1(q: f64, n: i32) -> Result<(Vec<f64>, Vec<f64>), RocheError> {
    // Accuracy of location of surface in terms of binary separation
    const FRAC: f64 = 1.0e-6;

    // Compute the potential at the inner Lagrange point
    let rl1 = x_l1(q)?;
    let p: Vec3 = Vec3::new(rl1, 0.0, 0.0);
    let cpot: f64 = rpot(q, &p)?;

    let mut xarr: Vec<f64> = Vec::with_capacity(n as usize);
    let mut yarr: Vec<f64> = Vec::with_capacity(n as usize);
    for i in 0..n {
        if i==0 || i==n-1 {
            // special case as derivative is zero at L1
            xarr.push(rl1);
            yarr.push(0.0);
        } else {
            let theta: f64 = (i as f64) * std::f64::consts::PI * 2.0 / ((n as f64) - 1.0);
            let dx = theta.cos();
            let dy = theta.sin();
            let line = LineRoche::new(q, Star::Primary, dx, dy, cpot);
            let lam = rtsafe(rl1/4.0, rl1, |lam| line.cost(lam), FRAC)?;
            xarr.push(lam*dx);
            yarr.push(lam*dy);
        }
    }
    Ok((xarr, yarr))
}


#[pyfunction]
pub fn lobe2(q: f64, n: i32) -> Result<(Vec<f64>, Vec<f64>), RocheError> {
    // Accuracy of location of surface in terms of binary separation
    const FRAC: f64 = 1.0e-6;

    // Compute the potential at the inner Lagrange point
    let rl1 = x_l1(q)?;
    let p: Vec3 = Vec3::new(rl1, 0.0, 0.0);
    let cpot: f64 = rpot(q, &p)?;
    let upper = 1.0 - rl1;
    let lower = upper/4.0;
    let mut xarr: Vec<f64> = Vec::with_capacity(n as usize);
    let mut yarr: Vec<f64> = Vec::with_capacity(n as usize);
    for i in 0..n {
        if i==0 || i==n-1 {
            // special case as derivative is zero at L1
            xarr.push(rl1);
            yarr.push(0.0);
        } else {
            let theta: f64 = (i as f64) * std::f64::consts::PI * 2.0 / ((n as f64) - 1.0);
            let dx = -theta.cos();
            let dy = theta.sin();
            let line = LineRoche::new(q, Star::Secondary, dx, dy, cpot);
            let lam = match rtsafe(lower, upper, |lam| line.cost(lam), FRAC) {
                Ok(lam) => lam,
                Err(e) => return Err(e),
            };
            xarr.push(1.0+lam*dx);
            yarr.push(lam*dy);
        }
    }
    Ok((xarr, yarr))
}

/// rtsafe is a Numerical Recipes-based routine to find roots
/// of a function using bisection or Newton-Raphson as appropriate.
/// \param func function object. Returns a tuple of (function value, derivative) at given x.
/// \param x1 value to the left of the root
/// \param x2 value to the right of the root
/// \param xacc minimum accuracy in returned root
/// \return Returns the x value of the root.
pub fn rtsafe<F>(
    x1: f64,
    x2: f64,
    func: F,
    xacc: f64,
) -> Result<f64, RocheError>
where
    F: Fn(f64) -> Result<(f64, f64), RocheError>,
{
    let mut xlo = x1;
    let mut xhi = x2;
    let mut fl;
    let mut fh;
    let mut df;
    const MAXITER: i32 = 100;
    (fl, _) = func(xlo)?;
    (fh, _) = func(xhi)?;

    if (fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0) {
        return Err(RocheError::RtsafeError("Root must be bracketed in rtsafe".to_string()));
    }

    // return if any of the endpoints is a root
    if fl == 0.0 {
        return Ok(xlo);
    } else if fh == 0.0 {
        return Ok(xhi);
    }

    // If fhi < 0.0, set things up so that xlo is below the root and xhi is above
    if fh < 0.0 {
        std::mem::swap(&mut xlo, &mut xhi);
        std::mem::swap(&mut fl, &mut fh);
    }

    let mut rts = 0.5 * (xlo + xhi);
    let mut dxold = (xhi - xlo).abs();
    let mut dx = dxold;
    let mut f;
    (f, df) = func(rts)?;
    let mut iter = 0;
    while iter < MAXITER {
        if ((rts - xhi) * df - f) * ((rts - xlo) * df - f) >= 0.0 
            || ((2.0*f).abs() > (dxold * df).abs()){
            // Bisect if Newton-Raphson is out of range or not decreasing fast enough
            dxold = dx;
            dx = 0.5 * (xhi - xlo);
            rts = xlo + dx;
            if xlo == rts {
                return Ok(rts);
            }
        } else {
            // Newton-Raphson step
            dxold = dx;
            dx = f / df;
            let temp = rts;
            rts -= dx;
            if temp == rts {
                return Ok(rts);
            }
        }

        if dx.abs() < xacc {
            return Ok(rts);
        }

        (f, df) = func(rts)?;
        if f < 0.0 {
            xlo = rts;
        } else {
            xhi = rts;
        }
        iter += 1;
    }

    return Err(RocheError::RtsafeError("Maximum number of iterations exceeded in rtsafe".to_string()));
}