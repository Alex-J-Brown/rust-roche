pub mod roche_context;
pub mod vec3;
pub mod blink;
pub mod x_lagrange;
pub mod pot_min;
pub mod set_earth;
pub mod potential;
pub mod zeta_rlobe_eggleton;
pub mod stream_physics;
pub mod sphere_eclipse;
pub mod ref_sphere;
pub mod fblink;
pub mod ingress_egress;
pub mod disc_eclipse;
pub mod face;
pub mod vel_transform;

pub use blink::blink;
pub use vec3::Vec3;
pub use x_lagrange::*;
pub use set_earth::*;
pub use potential::*;
pub use zeta_rlobe_eggleton::*;
pub use sphere_eclipse::*;
pub use ref_sphere::ref_sphere;
pub use pot_min::*;
pub use fblink::fblink;
pub use ingress_egress::ingress_egress;
pub use roche_context::RocheContext;
pub use disc_eclipse::disc_eclipse;
pub use face::face;

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
