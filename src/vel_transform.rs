use crate::errors::RocheError;


///
/// vtrans computes two velocity transforms, (1) a straight transform 
/// from rotating to inertial frame and (2) an inertial frame velocity 
/// in the disc.
/// 
/// Arguments:
/// 
/// * `q`: mass ratio M2/M1
/// * `transform_type`: integer representing velocity transform.
/// 1: rotating frame -> inertial frame
/// 2: rotating frame to inertial frame velocity of disc
/// 3: rotating frame (i.e. gived out what you entered)
/// * `x`: x-coordinate where transform is performed for
/// * `y`: y-coordinate where transform is performed for
/// * `vx`: velocity in x-axis in the rotating frame to transform
/// * `vy`: velocity in y-axis in the rotating frame to transform
/// 
/// Returns:
/// 
/// * `vx`: velocity in x-axis transformed to the chosen frame
/// * `vy`: velocity in y-axis tranformed to the chosen frame
/// 
pub fn vel_transform(q: f64, transform_type: i32, x: f64, y: f64, vx: f64, vy: f64) -> Result<(f64, f64), RocheError> {

    let mu: f64 = q/(1.0+q);
    let rad: f64 = (x*x + y*y).sqrt();
    let vkep: f64 = 1.0/((1.0 + q)*rad).sqrt();

    match transform_type {
        1 => Ok((vx - y, vy + x - mu)),
        2 => Ok((-vkep*y/rad, vkep*x/rad - mu)),
        3 => Ok((vx, vy)),
        _ => Err(RocheError::ParameterError(format!("{} is not a valid transform_type.", transform_type))),
    }
}