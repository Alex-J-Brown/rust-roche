use pyo3::prelude::*;

pub mod errors;

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
pub use ref_sphere::*;

pub mod stream_physics;
pub use stream_physics::*;

pub mod vstream_physics;
pub use vstream_physics::*;

pub mod fblink;
pub use fblink::*;

pub mod ingress_egress;
pub use ingress_egress::ingress_egress;

pub mod disc_eclipse;
pub use disc_eclipse::disc_eclipse;

pub mod star_eclipse;
pub use star_eclipse::star_eclipse;

pub mod face;
pub use face::*;

pub mod vel_transform;
pub use vel_transform::vel_transform;

pub mod solve_triads;
pub use solve_triads::findi;

pub mod lobes;
pub use lobes::lobe1;
pub use lobes::lobe2;

pub mod phases;
pub use phases::*;

pub mod jacobi;
pub use jacobi::*;

#[pyclass(from_py_object, eq, eq_int)]
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Star {
    Primary = 1,
    Secondary = 2,
}

pub type Etype = Vec<(f64, f64)>;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

// Python module

#[pymodule]
fn roche(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<vec3::Vec3>()?;
    m.add_class::<Star>()?;
    m.add_function(wrap_pyfunction!(set_earth::set_earth_iangle, m)?)?;
    m.add_function(wrap_pyfunction!(set_earth::set_earth, m)?)?;
    m.add_function(wrap_pyfunction!(x_lagrange::x_l1, m)?)?;
    m.add_function(wrap_pyfunction!(x_lagrange::x_l1_1, m)?)?;
    m.add_function(wrap_pyfunction!(x_lagrange::x_l1_2, m)?)?;
    m.add_function(wrap_pyfunction!(x_lagrange::x_l2, m)?)?;
    m.add_function(wrap_pyfunction!(x_lagrange::x_l3, m)?)?;
    m.add_function(wrap_pyfunction!(zeta_rlobe_eggleton::zeta_rlobe_eggleton, m)?)?;
    m.add_function(wrap_pyfunction!(zeta_rlobe_eggleton::dzetadq_rlobe_eggleton, m)?)?;
    m.add_function(wrap_pyfunction!(sphere_eclipse::sphere_eclipse_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(sphere_eclipse::sphere_eclipse_vector_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(ref_sphere::ref_sphere, m)?)?;
    m.add_function(wrap_pyfunction!(potential::rpot, m)?)?;
    m.add_function(wrap_pyfunction!(potential::rpot1, m)?)?;
    m.add_function(wrap_pyfunction!(potential::rpot2, m)?)?;
    m.add_function(wrap_pyfunction!(potential::drpot, m)?)?;
    m.add_function(wrap_pyfunction!(potential::drpot1, m)?)?;
    m.add_function(wrap_pyfunction!(potential::drpot2, m)?)?;
    m.add_function(wrap_pyfunction!(potential::rpot_val, m)?)?;
    m.add_function(wrap_pyfunction!(potential::rpot_val_grad, m)?)?;
    m.add_function(wrap_pyfunction!(ingress_egress::ingress_egress_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(fblink::fblink, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::stradv_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::rocacc, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::strinit, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::stream, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::streamr, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::strmnx_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(stream_physics::brightspot_position, m)?)?;
    m.add_function(wrap_pyfunction!(vstream_physics::vstream_reg, m)?)?;
    m.add_function(wrap_pyfunction!(face::face, m)?)?;
    m.add_function(wrap_pyfunction!(solve_triads::findi, m)?)?;
    m.add_function(wrap_pyfunction!(solve_triads::findq, m)?)?;
    m.add_function(wrap_pyfunction!(solve_triads::findphi, m)?)?;
    m.add_function(wrap_pyfunction!(lobes::lobe1, m)?)?;
    m.add_function(wrap_pyfunction!(lobes::lobe2, m)?)?;
    m.add_function(wrap_pyfunction!(lobes::vlobe1, m)?)?;
    m.add_function(wrap_pyfunction!(lobes::vlobe2, m)?)?;
    m.add_function(wrap_pyfunction!(phases::wdradius, m)?)?;
    m.add_function(wrap_pyfunction!(phases::wdphases, m)?)?;
    m.add_function(wrap_pyfunction!(phases::bsphases, m)?)?;
    m.add_function(wrap_pyfunction!(jacobi::jacobi, m)?)?;
    m.add_function(wrap_pyfunction!(planck::planck, m)?)?;
    m.add_function(wrap_pyfunction!(planck::dplanck, m)?)?;
    m.add_function(wrap_pyfunction!(planck::dlpdlt, m)?)?;
    Ok(())
}
