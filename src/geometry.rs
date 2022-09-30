use log::debug;
use ndarray::{arr1, Array1};

use crate::{
    halfplane::Halfplane,
    math::{angle2vec, closest_point_on_line, cross, norm, vec2angle},
};

/// Calculate the fastest way out of a disk.
pub fn out_of_disk(
    disk_center: &Array1<f64>,
    disk_r: f64,
    velocity: &Array1<f64>,
) -> (Array1<f64>, Array1<f64>) {
    debug!("out_of_disk({:?}, {}, {:?})", disk_center, disk_r, velocity);
    let rel_vec = velocity - disk_center;
    let w_length = norm(&rel_vec);
    // rotate vector to outside of disk by 10 degrees
    let n = angle2vec(vec2angle(&rel_vec) + 10.0);
    // calculate length of "u" (i.e., the way out of the disk)
    let u = &n * (disk_r - w_length);
    (u, n)
}

/// Intersect all given halfplanes and find the best velocity
pub fn halfplane_intersection(
    halfplanes_u: &Vec<Halfplane>,
    current_velocity: &Array1<f64>,
    optimal_point: &Array1<f64>,
) -> Option<Array1<f64>> {
    debug!(
        "halfplane_intersection({:?}, {:?}, {:?})",
        halfplanes_u, current_velocity, optimal_point
    );
    let mut halfplanes = Vec::new();
    for plane in halfplanes_u {
        halfplanes.push(Halfplane {
            u: current_velocity.clone() + plane.u.clone(),
            n: plane.n.clone(),
        });
    }
    let mut new_point = optimal_point.clone();

    for (i, plane) in halfplanes.iter().enumerate() {
        if (&new_point - &plane.u).dot(&plane.n) < 0.0 {
            // intersect this halfplane with all other halfplanes
            let (left, right) = intersect_halfplane_with_other_halfplanes(plane, &halfplanes[..i]);
            // check for empty intersection
            if left.is_none() || right.is_none() {
                return None;
            }
            // get the closes point on the new line (this is our new "preferred velocity")
            new_point = closest_point_on_line(
                &plane.u,
                &(&plane.u + &arr1(&[plane.n[1], -plane.n[0]])),
                optimal_point,
                left.unwrap(),
                right.unwrap(),
            )
        }
    }

    Some(new_point)
}

///  TODO: move this into impl of Halfplane
/// Intersect one halfplane with a set of halfplanes.
fn intersect_halfplane_with_other_halfplanes(
    plane: &Halfplane,
    other_planes: &[Halfplane],
) -> (Option<f64>, Option<f64>) {
    debug!(
        "intersect_halfplane_with_other_halfplanes({:#?}, {:#?})",
        plane, other_planes
    );
    let mut left = -f64::INFINITY;
    let mut right = f64::INFINITY;

    // calculate direction of halfplane
    let direction = arr1(&[plane.n[1], -plane.n[0]]);

    for other_plane in other_planes {
        // perform crazy math to intersect halfplane
        // See https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect for reference
        let other_dir = arr1(&[other_plane.n[1], -other_plane.n[0]]);
        let num = cross(&(&other_plane.u - &plane.u), &other_dir);
        let den = cross(&direction, &other_dir);

        if den < f64::EPSILON {
            if num < f64::EPSILON {
                // lines are parallel
                return (None, None);
            }
            continue;
        }

        let t = num / den;
        if den > 0.0 {
            right = right.min(t);
        } else {
            left = left.max(t);
        }

        if left > right {
            // no intersection
            // TODO: this somehow is buggy...
            return (None, None);
        }
    }

    (Some(left), Some(right))
}
