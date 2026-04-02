use crate::errors::RocheError;



pub fn vel_transform(q: f64, transform_type: i32, x: f64, y: f64, vx: f64, vy: f64) -> Result<(f64, f64), RocheError> {

    let mu: f64 = q/(1.0+q);
    let rad: f64 = (x*x + y*y).sqrt();
    let vkep: f64 = 1.0/((1.0 + q)*rad).sqrt();

    match transform_type {
        1 => Ok((vx - y, vy + x - mu)),
        2 => Ok((-vkep*y/rad, vkep*x/(rad-mu))),
        3 => Ok((vx, vy)),
        _ => Err(RocheError::ParameterError(format!("{} is not a valid transform_type.", transform_type))),
    }
}