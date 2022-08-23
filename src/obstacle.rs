use ndarray::Array1;

/// struct representing a static obstacle
#[derive(Debug)]
pub struct Obstacle {
    pub start: Array1<f64>,
    pub end: Array1<f64>,
    pub radius: f64,
}
