use crate::math;

/// Struct for representing a participant in the ORCA algorithm.
#[derive(Debug, Clone)]
pub struct Participant {
    pub position: ndarray::Array1<f64>,
    pub velocity: ndarray::Array1<f64>,
    pub radius: f64,
    pub confidence: f64,
    pub vmax: f64,
    pub target: ndarray::Array1<f64>,
    pub in_obstacle: bool,
}

impl Participant {
    pub fn update_position(&mut self, position: &ndarray::Array1<f64>) {
        self.position = position.clone();
        self.velocity = &self.target - &self.position;
        if math::norm(&self.velocity) > self.vmax {
            self.velocity = math::normalize(&self.velocity) * self.vmax
        }
    }

    pub fn in_obstacle(&mut self) {
        self.in_obstacle = true;
    }
}
