use crate::Vec3;
use crate::errors::RocheError;
use crate::x_l1;

///
/// blink tests for the occultation of a point in a semi-detached
/// binary. blink = true when point is eclipsed, otherwise = false. It works 
/// by first trying to eliminate as many paths that go nowhere near the 
/// Roche lobe before getting down to the hard work of stepping
/// along the path.
///
/// NB Note that this version of !!emph{blink} works in units of binary separation
/// NOT the inner lagrangian point distance.
///
/// Arguments:
/// 
/// * `q`: the mass ratio = M2/M1.
/// * `r`: the position vector (in units of binary separation)
/// * `e`: unit vector pointing towards Earth. Standardly this is (sini(i)*cos(phi),-sin(i)*sin(phi),cos(i))
/// * `acc`: step size parameter. acc specifies the size of steps taken when trying to see if the 
/// photon path goes inside the Roche lobe. The step size is roughly acc times the radius of the lobe filling star. 
/// This means that the photon path could mistakenly not be occulted if it passed less than about (acc**2)/8 of 
/// the radius below the surface of the Roche lobe. acc of order 0.1 should therefore do the job.
/// Returns:
/// 
/// * false = not eclipsed; true = eclipsed.
///
pub fn blink(q: f64, r: &Vec3, e: &Vec3, acc: f64) -> Result<bool, RocheError> {

    let mut x_cofm: f64;
    let c1: f64;
    let c2: f64;
    let step: f64;
    let pp: f64;
    let crit: f64;
    let mut p: f64;
    let p1: f64;
    let p2: f64;


    let mut r1: f64;
    let mut r2: f64;
    // Compute q dependent quantities
    if q <= 0. {
        let message = format!("q = {} <= 0", q);
        return Err(RocheError::ParameterError(message));
    }

    if acc <= 0.0 {
        let message = format!("Invalid accuracy parameter. {} <= 0.0.", acc);
        return Err(RocheError::ParameterError(message));
    }
 
    x_cofm = 1.0 / (1.0 + q);
    c1 = 2.0 * x_cofm;
    x_cofm *= q;
    c2 = 2.0 * x_cofm;

    // Locate the inner Lagrangian point (L1)
    let rl1: f64 = x_l1(q)?;

    // Evaluate Roche potential at L1 point.
    r1 = rl1;
    r2 = 1.0 - rl1;
    let xc: f64 = rl1 - x_cofm;
    crit = c1/r1 + c2/r2 + xc*xc;

    // The red star lies entirely within the sphere centred on its
    // centre of mass and reaching the inner Lagrangian point.

    let rsphere: f64 = 1.0 - rl1;
    pp = rsphere*rsphere;
    step = rsphere*acc;

    // From now on computations are done every call. Main point is
    // to try to bail out as soon as possible to save time 

    // evaluate closest approach distance to sphere.

    let xt: f64 = r.x - 1.0;
    let mut b: f64 = e.x*xt + e.y*r.y + e.z*r.z;
    let mut c: f64 = xt*xt + r.y*r.y + r.z*r.z - pp;

    // Photon path crosses sphere at two points given by quadratic
    // equation l*l + 2*b*l + c = 0. First check that this has
    // real roots. If not, there is no eclipse. This test should
    // eliminate most cases.

    let mut fac: f64 = b*b - c;

    if fac <= 0.0 {
        return Ok(false);
    }

    fac = fac.sqrt();

    // If the larger root is negative, the photon starts after the sphere.
    // This should get rid of another whole stack
    

    let par2: f64 = -b + fac;

    if par2 <= 0.0 {
        return Ok(false);
    }

    //  Now for the hard work. The photon's path does go
    //  inside the sphere ... groan. First evaluate smaller root 
    //  but limit to >= 0 to ensure that only photon path
    //  and not its extrapolation backwards in included.

    let mut par1: f64 = -b - fac;
    par1 = if par1 > 0.0 {
        par1
    } else {
        0.0
    };

    // Now follow the photon's path in finite steps to see if 
    // it intercepts the red star. We start off at closest
    // approach point to increase chance of only needing one 
    // computation of the roche potential

    let par: f64 =  if b < 0.0 {
        -b
    } else {
        0.0
    };

    let mut x1: f64 = r.x + par*e.x;
    let mut x2: f64 = r.y + par*e.y;
    let mut x3: f64 = r.z + par*e.z;

    // Test roche potential for an occultation
    let mut xm: f64 = x1 - 1.0;
    let mut yy: f64 = x2*x2;
    let mut rr: f64 = yy + x3*x3;
    let rs2: f64 = xm*xm + rr;

    // Point at c of m of red star. Definitely eclipsed, avoids
    // division by 0 later

    if rs2 <= 0.0 {
        return Ok(true);
    }

    let rs1: f64 = x1*x1 + rr;
    r1 = rs1.sqrt();
    r2 = rs2.sqrt();

    // Deeper in well than inner Lagrangian, therefore eclipsed.

    let mut xc: f64 = x1 - x_cofm;
    c = c1/r1 + c2/r2 + xc*xc + yy;

    if c > crit {
        return Ok(true);
    }

    // Now we need to step. Determine step direction by 
    // evaluating the first derivative

    let a: f64 = x1*e.x + x2*e.y;
    b = a + x3*e.z;

    if -c1*b/(rs1*r1) - c2*(b - e.x)/(rs2*r2) + 2.0*(a - x_cofm*e.x) > 0.0 {
        p1 = par;
        p2 = par2;
    } else {
        p1 = par;
        p2 = par1;
    }

    // Loop while the Roche potential increases in depth

    let mut nstep: i32 = ((p2-p1).abs() / step + 0.5).floor() as i32;

    nstep = nstep.max(2);
    
    let dp: f64 = (p2-p1)/(nstep as f64);
    let mut cmax: f64 = c - 1.0;
    let mut i: i32 = 0;
    while c > cmax && i < nstep {
        i += 1;
        cmax = c;
        p = p1 + dp*(i as f64);
        x1 = r.x + p*e.x;
        x2 = r.y + p*e.y;
        x3 = r.z + p*e.z;

        // test roche potential for an occultation, again guarding
        // against division by 0

        xm = x1 - 1.0;
        yy = x2*x2;
        rr = yy + x3*x3;
        r2 = (xm*xm + rr).sqrt();
        if r2 <= 0.0 {
            return Ok(true);
        }
        r1 = (x1*x1 + rr).sqrt();
        xc = x1 - x_cofm;
        c = c1/r1 + c2/r2 + xc*xc + yy;
        if c > crit {
            return Ok(true);
        }
    }

    return Ok(false)

}