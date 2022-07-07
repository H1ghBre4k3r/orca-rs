use ndarray::Array1;

/// struct representing a static obstacle
pub struct Obstacle {
    pub start: Array1<f64>,
    pub end: Array1<f64>,
    pub radius: f64,
}
