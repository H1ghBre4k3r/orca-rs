use ndarray::{arr1, Array1};

pub fn norm(vector: &ndarray::Array1<f64>) -> f64 {
    return vector.dot(vector).sqrt();
}

pub fn r2d(rad: f64) -> f64 {
    return rad.to_degrees();
}

pub fn d2r(deg: f64) -> f64 {
    return deg.to_radians();
}

pub fn vec2angle(vector: &ndarray::Array1<f64>) -> f64 {
    let re = &arr1(&[1.0, 0.0]);
    let re_unit = re / norm(re);

    let vec_len = norm(&vector);
    let mut x_unit = arr1(&[0.0, 0.0]);
    if vec_len != 0.0 {
        x_unit = vector / vec_len;
    }
    let mut ang = r2d(re_unit.dot(&x_unit).acos());
    if vector[1] < 0.0 {
        ang *= -1.0
    }
    return ang;
}

pub fn angle2vec(angle: f64) -> ndarray::Array1<f64> {
    return arr1(&[d2r(angle).cos(), d2r(angle).sin()]);
}

pub fn arcsin(gegen_kat: f64, hypo: f64) -> f64 {
    return r2d((gegen_kat / hypo).asin());
}

pub fn angle_diff(a: f64, b: f64) -> f64 {
    let x = (a + 360.0) % 360.0;
    let y = (b + 360.0) % 360.0;
    return (x - y)
        .abs()
        .min((x - y - 360.0).abs())
        .min((x - y + 360.0).abs());
}

pub fn dist(x: &ndarray::Array1<f64>, y: &ndarray::Array1<f64>) -> f64 {
    return norm(&(x - y));
}

pub fn mix(
    a: &ndarray::Array1<f64>,
    b: &ndarray::Array1<f64>,
    amount: f64,
) -> ndarray::Array1<f64> {
    return a.clone() + (b - a) * amount;
}

pub fn closest_point_on_line(
    l0: &ndarray::Array1<f64>,
    l1: &ndarray::Array1<f64>,
    tar: &ndarray::Array1<f64>,
    left: f64,
    right: f64,
) -> ndarray::Array1<f64> {
    let c = tar - l0;
    let mut v = l1 - l0;
    let distance = norm(&v);
    if distance == 0.0 {
        return l0.clone();
    }
    v /= distance;
    let d = norm(&(l0 - l1));
    let t = c.dot(&v) / d;
    return mix(l0, l1, t.clamp(left, right));
}

pub fn normalize(vec: &ndarray::Array1<f64>) -> ndarray::Array1<f64> {
    let len = norm(vec);
    if len == 0.0 {
        return ndarray::arr1(&[0.0, 0.0]);
    }
    return vec / len;
}

pub fn cross(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    return a[0] * b[1] - a[1] * b[0];
}

pub fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}
