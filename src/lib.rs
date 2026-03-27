pub mod roche_context;
pub use roche_context::RocheContext;

pub mod vec3;
pub use vec3::Vec3;

pub mod point;
pub use point::Point;

pub mod planck;
pub use planck::*;

pub mod blink;
pub use blink::blink;

pub mod x_lagrange;
pub use x_lagrange::*;

pub mod pot_min;
pub use pot_min::*;

pub mod set_earth;
pub use set_earth::*;

pub mod potential;
pub use potential::*;

pub mod zeta_rlobe_eggleton;
pub use zeta_rlobe_eggleton::*;

pub mod sphere_eclipse;
pub use sphere_eclipse::*;

pub mod ref_sphere;
pub use ref_sphere::ref_sphere;

pub mod stream_physics;
pub use stream_physics::*;

pub mod fblink;
pub use fblink::fblink;

pub mod ingress_egress;
pub use ingress_egress::ingress_egress;

pub mod disc_eclipse;
pub use disc_eclipse::disc_eclipse;

pub mod star_eclipse;
pub use star_eclipse::star_eclipse;

pub mod face;
pub use face::face;

pub mod vel_transform;
pub use vel_transform::vel_transform;


#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Star {
    Primary,
    Secondary,
}

pub type Etype = Vec<(f64, f64)>;


pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
