use ndarray::Array1;

/// struct representing a Halfplane in Hessian normal form
#[derive(Debug)]
pub struct Halfplane {
    pub u: Array1<f64>,
    pub n: Array1<f64>,
}

impl Clone for Halfplane {
    fn clone(&self) -> Self {
        Halfplane {
            u: self.u.clone(),
            n: self.n.clone(),
        }
    }
}
