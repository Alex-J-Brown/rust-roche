use std::f64::consts::TAU;
use crate::errors::RocheError;
use crate::{Vec3, Star};
use crate::{rpot_grad, rpot_val, set_earth};


///
///  linmin minimises along a line in phase, lambda space. It returns at the
/// minimum or as soon as the potential drops below a reference value. It is
/// assumed that the potential is dropping with x at the starting point.
///
/// Arguments:
/// 
/// * `q`: mass ratio = M2/M1
/// * `star`: star in question
/// * `spin`: spin to orbital frequency ratio
/// * `cosi`: cosine of orbital inclination (both passed to speed computations)
/// * `sini`: sine of orbital inclination
/// * `p`: point of origin
/// * `phi`: value of phase at start (0 - 1) and at end (returned)
/// * `lam`: value of lambda at start and at end (returned)
/// * `dphi`: rate at which phi changes
/// * `dlam`: rate at which lambda changes
/// * `phi1`: minimum value of phi
/// * `phi2`: maximum value of phi
/// * `lam1`: minimum lambda
/// * `lam2`: maximum lambda
/// * `pref`: reference potential.
/// * `acc`: accuracy in position
/// 
/// Returns:
/// 
/// * `pmin`: value of function at minimum
/// * `jammed`: true if minimum is on a boundary
/// 
pub fn linmin(q: f64, star: Star, spin: f64, cosi: f64, sini: f64, p: &Vec3, phi: &mut f64, lam: &mut f64, mut dphi: f64, mut dlam: f64, phi1: f64, phi2: f64, lam1: f64, lam2: f64, pref: f64, acc: f64) -> Result<(f64, bool), RocheError> {

    let mut jammed = false;

    let make_func = |phi0: f64, lam0: f64, dphi0: f64, dlam0: f64| {
        
        let func = move |x: f64| {
            let earth: Vec3 = set_earth(cosi, sini, phi0 + dphi0 * x);
            Ok(rpot_val(q, star, spin, &earth, p, lam0 + dlam0 * x)?)
        };

        let dfunc = move |x: f64| {
            let earth: Vec3 = set_earth(cosi, sini, phi0 + dphi0 * x);
            let (dp, dl) = rpot_grad(q, star, spin, &earth, p, lam0 + dlam0 * x)?;
            Ok(dp * dphi0 + dl * dlam0)
        };

        (func, dfunc)
    };

    let (mut func, mut dfunc) = make_func(*phi, *lam, dphi, dlam);

    // --- determine xmax from boundaries ---
    let mut xmax: f64 = 1.0e30;
    let mut nbound: i32 = 0;

    let mut check = |bound: f64, val: f64, d: f64, id: i32| {
        if d != 0.0 {
            let x: f64 = (bound - val) / d;
            if x > 0.0 && x < xmax {
                xmax = x;
                nbound = id;
            }
        }
    };

    check(phi1, *phi, dphi, 1);
    check(phi2, *phi, dphi, 2);
    check(lam1, *lam, dlam, 3);
    check(lam2, *lam, dlam, 4);

    // --- initial bracketing ---
    // Now the aim is to bracket the minimum, while accounting for the maximum
    // possible step so that we can then apply dbrent.
    let xa: f64 = 0.0;
    let fa: f64 = func(xa)?;

    let mut xb: f64 = 1e-8 * xmax;
    let mut fb: f64 = func(xb)?;

    // for _ in 0..7 {
    //     if fb < fa && xa != xb {
    //         break;
    //     }
    //     xb *= 10.0;
    //     fb = func(xb);
    // }
    let mut nten: i32 = 0;
    const NTEN: i32 = 7;
    while (fb >= fa || xa == xb) && nten < NTEN {
        nten += 1;
        xb *= 10.0;
        fb = func(xb)?;
    }

    if fb <= pref {
        *phi += dphi * xb;
        *lam += dlam * xb;
        return Ok((fb, jammed));
    }

    if fb >= fa {
        // Let's hope that we have not stepped past the minimum without
        // knowing it
        return Ok((fa, jammed));
    }

    // --- bracket other side ---
    // OK, so fb < fa so we are heading downhill at least. Now try
    // to find other side starting from xb, looking for a point when
    // we go up or dip below the critical potential
    let mut bracketted = false;
    
    let mut xc: f64 = 0.0;
    let mut fc: f64 = 0.0;

    let xbold: f64 = xb;
    let fbold: f64 = fb;
    xmax -= xb;

    const NTRY: i32 = 5;

    let xmax_rem: f64 = xmax - xb;

    for n in 1..=NTRY {
        xc = xbold + xmax_rem * (n as f64) / (NTRY as f64);
        fc = func(xc)?;

        if fc <= pref {
            *phi += dphi * xc;
            *lam += dlam * xc;
            return Ok((fc, jammed));
        }

        if fc < fb {
            xb = xc;
            fb = fc;
        } else {
            bracketted = true;
            break;
        }
    }

    jammed = false;

    // --- boundary logic ---
    if !bracketted {

        let dc: f64 = dfunc(xc)?;

        if dc > 0.0 {

            // We have crashed into the end stop without crossing the minimum
            // but the derivative says that we are going up. Damn!  Go back to
            // old xb and check that derivative was going down there it really
            // ought to have been ...
            xb = xbold;
            let db: f64 = dfunc(xb)?;
            if db < 0.0 {

                // OK, let's try to zero in on the point at which the derivative
                // switches sign.
                let mut xm: f64 = (xb+xc)/2.0;
                while (xc-xb) > 1.0e-8*xc {
                    xm = (xb+xc)/2.0;
                    if dfunc(xm)? > 0.0 {
                        xc = xm;
                    } else {
                        xb = xm;
                    }
                }
                let fm: f64 = func(xm)?;
                if fm <= pref {
                    *phi += dphi*xm;
                    *lam += dlam*xm;
                    return Ok((fm, jammed));
                }
                if fm < fc && fm < fbold {
                    xb = xm;
                } else {
                    return Err(RocheError::LinminError("Failed to bracket minimum, error 3.".to_string()))
                }
            } else {
                return Err(RocheError::LinminError("Failed to bracket minimum, error 1.".to_string()))
            }
        } else {
            // We are trapped on a boundary; re-define line minimisation functions.
            jammed = true;
            *phi += dphi*xmax;
            *lam += dlam*xmax;
            xmax = 1.0;
            (dphi, dlam) = rpot_grad(q, star, spin, &set_earth(cosi, sini, *phi), p, *lam)?;

            match nbound {
                1 => {
                    *phi = phi1;
                    dphi = 0.0;
                    if dlam > 0.0 {
                        dlam = lam1 - *lam;
                    } else {
                        dlam = lam2 - *lam;
                    }
                }
                2 => {
                    *phi = phi2;
                    dphi = 0.0;
                    if dlam > 0.0 {
                        dlam = lam1 - *lam;
                    } else {
                        dlam = lam2 - *lam;
                    }
                }
                3 => {
                    *lam = lam1;
                    dlam = 0.0;
                    if dphi > 0.0 {
                        dphi = phi1 - *phi;
                    } else {
                        dphi = phi2 - *phi;
                    }
                }
                4 => {
                    *lam = lam1;
                    dlam = 0.0;
                    if dphi > 0.0 {
                        dphi = phi1 - *phi;
                    } else {
                        dphi = phi2 - *phi;
                    }
                }
                _ => {}
            }
            (func, dfunc) = make_func(*phi, *lam, dphi, dlam);

            // Again try to bracket minimum
            nten = 0;
            xb = 1.0e-6;
            fb = func(xb)?;

            while fb >= fa && nten < NTEN {
                xb *= 10.0;
                nten += 1;
                fb = func(xb)?;
            }
            if fb <= pref {
                *phi += dphi*xb;
                *lam += dlam*xb;
                return Ok((fb, jammed));
            }
            if fb >= fa {
                return Ok((fa, jammed));
            }

            // Now the bracketting steps.
            bracketted = false;
            for n in 1..=NTRY {
                xc = xmax*n as f64/NTRY as f64;
                fc = func(xc)?;
                if fc <= pref {
                    *phi += dphi*xc;
                    *lam += dlam*xc;
                    return Ok((fc, jammed))
                } else {
                    bracketted = true;
                    break;
                }
            }

            if !bracketted {

                if dfunc(xc)? > 0.0 {
                    return Err(RocheError::LinminError("Failed to bracket minimum, error 2.".to_string()))
                } else {
                    *phi += dphi*xmax;
                    *lam += dlam*xmax;
                    return Ok((fc, jammed));
                }
            }

        }

    }

    // --- Brent refinement ---

    let xacc: f64 = acc / ((TAU*dphi).powi(2) + dlam*dlam).sqrt();

    let (xmin, pmin) = dbrent(xa, xb, xc, &func, &dfunc, xacc, true, pref)?;

    *phi += dphi * xmin;
    *lam += dlam * xmin;

    Ok((pmin, jammed))
}


/// 
/// The line of sight to any fixed point in a binary sweeps out a cone at the
/// binary rotates. Positions on the cone can be parameterised by the orbital
/// phase phi and the multiplier ('lambda') needed to get from the fixed point.
/// The question pot_min tries to solve is "does the cone intersect a surface
/// of fixed Roche potential lying within a Roche lobe?". It does so by
/// minimisation over a region of phi and lambda. It stops as soon as any
/// potential below a critical value is found. The initial range of phi and
/// lambda can be determined using sphere_eclipse which calculates them for a
/// sphere.
///
/// Arguments:
/// 
/// * `q`:  mass ratio = M2/M1.
/// * `cosi`: cosine orbital inclination.
/// * `sini`: sine orbital inclination.
/// * `star`: which star (needed for asynchronous case).
/// * `spin`: ratio of spin/orbital.
/// * `p`: point of origin.
/// * `phi1`: minimum phase within which eclipse may occur (0 - 1).
/// * `phi2`: maximum phase within which an eclipse may occur (> phi1).
/// * `lam1`: minimum multiplier line of sight crosses eclipsing star.
/// * `lam2`: maximum multiplier line of sight crosses eclipsing star.
/// * `rref`: reference radius.
/// * `pref`: reference potential.
/// * `acc`:  absolute accuracy in position to go for.
/// * `phi`:  phi at minimum potential. Ingress occurs between phi1 and phi.
///           if there is an eclipse. Egress occurs between phi and phi2.
/// * `lam`:  lambda at minimum potential.
/// 
/// Returns:
/// 
/// * true if minimum potential is below the reference.
///
pub fn pot_min(q: f64, star: Star, spin: f64, cosi: f64, sini: f64, p: &Vec3, phi1: f64, phi2: f64, lam1: f64, lam2: f64, rref: f64, pref: f64, acc: f64, phi: &mut f64, lam: &mut f64) -> Result<bool, RocheError> {

    *phi = (phi1 + phi2)/2.;
    *lam = (lam1 + lam2)/2.;

    let mut rp: f64 = TAU * *phi;
    let (sinp, cosp) = rp.sin_cos();

    let mut earth: Vec3 = Vec3::new(sini*cosp, -sini*sinp, cosi);
    let mut pot: f64;
    let mut dphi: f64;
    let mut dlam: f64;
    pot = rpot_val(q, star, spin, &earth, p, *lam)?;
    (dphi, dlam) = rpot_grad(q, star, spin, &earth, p, *lam)?;
    if pot <= pref {
        return Ok(true)
    };

    let mut gdphi: f64 = -dphi;
    let mut gdlam: f64 = -dlam;
    dphi = gdphi;
    let mut hdphi: f64 = gdphi;
    dlam = gdlam;
    let mut hdlam: f64 = gdlam;

    const ITMAX: i32 = 200;

    let delphi: f64 = q/(1.0 + q)*((acc*acc)/(rref*rref))/2.0;
    let mut pmin: f64;
    let mut gam: f64;
    let mut dgg: f64;
    let mut gg: f64;
    let mut jammed: bool;
    for _ in 0..ITMAX {

        (pmin, jammed) = linmin(q, star, spin, cosi, sini, p, phi, lam, dphi, dlam, phi1, phi2, lam1, lam2, pref, acc)?;

        if pmin <= pref {
            return Ok(true)
        }
        if jammed || (pmin - pot).abs() < delphi {
            return Ok(false)
        }

        pot = pmin;
        rp = TAU * *phi;
        let (sinp, cosp) = rp.sin_cos();
        earth.set(sini*cosp, -sini*sinp, cosi);
        (dphi, dlam) = rpot_grad(q, star, spin, &earth, p, *lam)?;

        gg = gdphi*gdphi + gdlam*gdlam;
        if gg == 0. {
            return Ok(false)
        }

        dgg = (dphi+gdphi)*dphi + (dlam+gdlam)*dlam;
        gam = dgg/gg;

        gdphi = -dphi;
        gdlam = -dlam;
        dphi = gdphi + gam*hdphi;
        hdphi = gdphi + gam*hdphi;
        dlam = gdlam + gam*hdlam;
        hdlam = gdlam + gam*hdlam;

    }
    Err(RocheError::PotminError("Too many iterations.".to_string()))
}


/// 
/// Given a bracketted minimum, this routine refines the minimum, making use of derivatives. 
/// It comes from NR, with the added option of a quick bail out. ax to cx must bracket a 
/// minimum otherwise you will be in trouble.
/// 
/// Arguments:
/// 
/// * `ax`: one extreme of x.
/// * `bx`: mid x value.
/// * `cx`: other extreme x value.
/// * `func`: function object returning function value at x.
/// * `dfunc`: function object retunring derivative value at x.
/// * `acc`: (absolute) accuracy in x.
/// * `stopfast`: retrun as soon as function value goes below a reference value or not.
/// * `fref`: reference function value if stopfast = true.
/// 
/// Returns:
/// 
/// * x at minimum.
/// * function value at minimum.
///
pub fn dbrent<F, G>(
    ax: f64,
    bx: f64,
    cx: f64,
    func: F,
    dfunc: G,
    acc: f64,
    stopfast: bool,
    fref: f64,
) -> Result<(f64, f64), RocheError>
where
    F: Fn(f64) -> Result<f64, RocheError>,
    G: Fn(f64) -> Result<f64, RocheError>,
{
    const ITMAX: usize = 100;

    let mut a: f64 = ax.min(cx);
    let mut b: f64 = ax.max(cx);

    let mut x: f64 = bx;
    let mut w: f64 = bx;
    let mut v: f64 = bx;

    let mut fx: f64 = func(x)?;
    let mut fw: f64 = fx;
    let mut fv: f64 = fx;

    if stopfast && fx < fref {
        return Ok((x, fx));
    }

    let mut dx: f64 = dfunc(x)?;
    let mut dw: f64 = dx;
    let mut dv: f64 = dx;

    let mut e: f64 = 0.0;
    let mut d: f64 = 0.0;

    for _ in 0..ITMAX {

        let xm: f64 = 0.5 * (a + b);
        let tol1: f64 = acc;
        let tol2: f64 = 2.0 * tol1;

        if (x - xm).abs() <= (tol2 - 0.5 * (b - a)) {
            return Ok((x, fx));
        }

        let mut d1: f64;
        let mut d2: f64;

        if e.abs() > tol1 {

            d1 = 2.0 * (b - a);
            d2 = d1;

            if dw != dx {
                d1 = (w - x) * dx / (dx - dw);
            }

            if dv != dx {
                d2 = (v - x) * dx / (dx - dv);
            }

            let u1: f64 = x + d1;
            let u2: f64 = x + d2;

            let ok1: bool = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
            let ok2: bool = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;

            let olde: f64 = e;
            e = d;

            if ok1 || ok2 {

                d = if ok1 && ok2 {
                    if d1.abs() < d2.abs() { d1 } else { d2 }
                } else if ok1 {
                    d1
                } else {
                    d2
                };

                if d.abs() <= 0.5 * olde.abs() {

                    let u: f64 = x + d;

                    if (u - a) < tol2 || (b - u) < tol2 {
                        d = tol1.copysign(xm - x);
                    }

                } else {
                    e = if dx >= 0.0 { a - x } else { b - x };
                    d = 0.5 * e;
                }

            } else {

                e = if dx >= 0.0 { a - x } else { b - x };
                d = 0.5 * e;

            }

        } else {

            e = if dx >= 0.0 { a - x } else { b - x };
            d = 0.5 * e;

        }

        let u: f64;
        let fu: f64;

        if d.abs() >= tol1 {

            u = x + d;
            fu = func(u)?;

            if stopfast && fu < fref {
                return Ok((u, fu));
            }

        } else {

            u = x + tol1.copysign(d);
            fu = func(u)?;

            if stopfast && fu < fref {
                return Ok((u, fu));
            }

            if fu > fx {
                return Ok((x, fx));
            }
        }

        let du: f64 = dfunc(u)?;

        if fu <= fx {

            if u >= x { a = x } else { b = x };

            v = w; fv = fw; dv = dw;
            w = x; fw = fx; dw = dx;
            x = u; fx = fu; dx = du;

        } else {

            if u < x { a = u } else { b = u };

            if fu <= fw || w == x {

                v = w; fv = fw; dv = dw;
                w = u; fw = fu; dw = du;

            } else if fu < fv || v == x || v == w {

                v = u; fv = fu; dv = du;

            }
        }
    }

    return Err(RocheError::DbrentError("too many iterations.".to_string()))
}

