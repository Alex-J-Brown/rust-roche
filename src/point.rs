use crate::Vec3;
use crate::Etype;



///
/// Struct defining a point on the surface of a model grid (e.g. of a star or disc etc.)
/// A `Point` has a position, a direction which is the surface normal, an area,
/// a relative gravity, a vector of phase pairs defining when the point is eclipsed by
/// another model component, and a flux.
/// 
#[derive(Clone, Debug)]
pub struct Point {
    pub position: Vec3,
    pub direction: Vec3,
    pub area: f32,
    pub gravity: f32,
    pub eclipse: Etype,
    pub flux: f32,
}

impl Point {

    ///
    /// Creates a new Point.
    /// 
    pub fn new(position: Vec3, direction: Vec3, area: f64, gravity: f64, eclipse: Etype) -> Self {
        Self {
            position,
            direction,
            area: area as f32,
            gravity: gravity as f32,
            eclipse,
            flux: 0.0,
        }
    }
    ///
    /// sets the point's flux.
    /// 
    pub fn set_flux(&mut self, flux: f32) -> () {
        self.flux = flux;
    }
    
    ///
    ///checks that the given phase is not during one of the
    /// phase ranges when the point is eclipsed.
    ///  
    pub fn is_visible(&self, phase: f64) -> bool {
        let phi: f64 = phase - phase.floor();
        for &(p1, p2) in &self.eclipse {
            if (phi >= p1 && phi <= p2) || phi <= p2 - 1.0 {
                return false
            }
        }
        true
    }

    
    ///
    /// This version of is_visible will not correct for phases outside
    /// of expected range to speed up large loops.
    /// run phase = phase - phase.floor();
    /// outside of loop beforehand
    /// 
    pub fn is_visible_phase_normed(&self, phase: f64) -> bool {
        for &(p1, p2) in &self.eclipse {
            if (phase >= p1 && phase <= p2) || phase <= p2 - 1.0 {
                return false
            }
        }
        true
    }
}


impl Default for Point {
    fn default() -> Self {
        Self::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
            vec![(0.0, 0.0)],
            )
    }
}