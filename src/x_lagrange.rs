///
/// x_l1 calculates the x coordinate of the L1 Lagrangian point in terms
/// of the orbital separation of the two stars which have a mass ratio, q.
/// It works by solving for the root of a quintic polynomial by Newton-Raphson
/// iteration. L1 is the point in between the two stars and so will be between
/// 0 and 1.
/// 
///  \param q mass ratio = M2/M1
/// 
pub fn x_l1(q: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let mu: f64 = q/(1. + q);
    let a1: f64 = -1.0 + mu;
    let a2: f64 =  2.0 - 2.*mu;
    let a3: f64 = -1.0 + mu;
    let a4: f64 =  1.0 + 2.*mu;
    let a5: f64 = -2.0 - mu;
    let a6: f64 = 1.0;
    let d1: f64 = 1.0*a2;
    let d2: f64 = 2.0*a3;
    let d3: f64 = 3.0*a4;
    let d4: f64 = 4.0*a5;
    let d5: f64 = 5.0*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1./(1. + q);
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


///
/// x_l1_1 calculates the x coordinate of the L1 Lagrangian point in terms
/// of the orbital separation of the two stars which have a mass ratio, q,
/// accounting for asynchronous rotation of the primary.
/// It works by solving for the root of a quintic polynomial by Newton-Raphson
/// iteration. L1 is the point in between the two stars and so will be between
/// 0 and 1.
/// 
///  \param q mass ratio = M2/M1
///  \param spin
/// 
pub fn x_l1_1(q: f64, spin: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let spin_squared: f64 = spin*spin;
    let mu: f64 = q/(1. + q);
    let a1: f64 = -1. + mu;
    let a2: f64 =  2. - 2.*mu;
    let a3: f64 = -1. + mu;
    let a4: f64 =  spin_squared + 2.*mu;
    let a5: f64 = -2.*spin_squared - mu;
    let a6: f64 = spin_squared;
    let d1: f64 = 1.*a2;
    let d2: f64 = 2.*a3;
    let d3: f64 = 3.*a4;
    let d4: f64 = 4.*a5;
    let d5: f64 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1./(1. + q);
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


///
/// x_l1_2 calculates the x coordinate of the L1 Lagrangian point in terms
/// of the orbital separation of the two stars which have a mass ratio, q,
/// accounting for asynchronous rotation of the secondary.
/// It works by solving for the root of a quintic polynomial by Newton-Raphson
/// iteration. L1 is the point in between the two stars and so will be between
/// 0 and 1.
/// 
///  \param q mass ratio = M2/M1
///  \param spin
/// 
pub fn x_l1_2(q: f64, spin: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let spin_squared: f64 = spin*spin;
    let mu: f64 = q/(1. + q);
    let a1: f64 = -1. + mu;
    let a2: f64 =  2. - 2.*mu;
    let a3: f64 = -spin_squared + mu;
    let a4: f64 =  3.*spin_squared + 2.*mu - 2.;
    let a5: f64 = 1. - mu - 3.*spin_squared;
    let a6: f64 = spin_squared;
    let d1: f64 = 1.*a2;
    let d2: f64 = 2.*a3;
    let d3: f64 = 3.*a4;
    let d4: f64 = 4.*a5;
    let d5: f64 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1./(1. + q);
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


///
/// x_l2 calculates the x coordinate of the L2 Lagrangian point in terms
/// of the orbital separation of the two stars which have a mass ratio, q,
/// accounting for asynchronous rotation of the primary.
/// It works by solving for the root of a quintic polynomial by Newton-Raphson
/// iteration. L2 is the point on the side of the secondary opposite the primary,
/// ands so x_l2 > 1.
/// 
///  \param q mass ratio = M2/M1
/// 
pub fn x_l2(q: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let mu = q/(1. + q);
    let a1 = -1. + mu;
    let a2 =  2. - 2.*mu;
    let a3 = -1. - mu;
    let a4 =  1. + 2.*mu;
    let a5 = -2. - mu;
    let a6 = 1.;
    let d1 = 1.*a2;
    let d2 = 2.*a3;
    let d3 = 3.*a4;
    let d4 = 4.*a5;
    let d5 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = 1.5;
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}


///
/// x_l3 calculates the x coordinate of the L3 Lagrangian point in terms
/// of the orbital separation of the two stars which have a mass ratio, q,
/// accounting for asynchronous rotation of the primary.
/// It works by solving for the root of a quintic polynomial by Newton-Raphson
/// iteration. L3 is the point on the side of the Primary opposite the secondary,
/// ands so x_l3 < 0.
/// 
///  \param q mass ratio = M2/M1
/// 
pub fn x_l3(q: f64) -> f64 {

    const NMAX: i32 = 1000;
    const EPS: f64 = 1.0e-12;

    if q <= 0. {
        panic!("q = {} <= 0", q);
    }

    let mu = q/(1. + q);
    let a1 = -1. + mu;
    let a2 =  -2. + 2.*mu;
    let a3 = 1. - mu;
    let a4 =  1. + 2.*mu;
    let a5 = -2. - mu;
    let a6 = 1.;
    let d1 = 1.*a2;
    let d2 = 2.*a3;
    let d3 = 3.*a4;
    let d4 = 4.*a5;
    let d5 = 5.*a6;

    let mut n: i32 = 0;
    let mut xold: f64 = 0.;
    let mut x: f64 = -1.;
    let mut f: f64;
    let mut df: f64;
    while n < NMAX && (x - xold).abs() > EPS*x.abs() {
        xold = x;
        f    = x*(x*(x*(x*(x*a6+a5)+a4)+a3)+a2)+a1;
        df   = x*(x*(x*(x*d5+d4)+d3)+d2)+d1;
        x   -= f/df;
        n+=1;
    }
    return x
}