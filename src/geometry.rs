use ndarray::{arr1, Array1};

use crate::{
    halfplane::Halfplane,
    math::{
        angle2vec, angle_diff, arcsin, closest_point_on_line, cross, norm, normalize, vec2angle,
    },
    obstacle::Obstacle,
    participant::Participant,
};

/// Calculate u and n for a static obstacle.
pub fn obstacle_collision(
    a: &Participant,
    obstacle: &Obstacle,
    tau: f64,
) -> (Array1<f64>, Array1<f64>) {
    let r_vec = normalize(&a.velocity) * (a.radius + a.confidence + obstacle.radius);
    let start = &obstacle.start - &(&a.position + &r_vec);
    let end = &obstacle.end - &(&a.position - &r_vec);

    let dist_vec =
        closest_point_on_line(&obstacle.start, &obstacle.end, &a.position, 0.0, 1.0) - &a.position;

    let u: Array1<f64>;
    let n: Array1<f64>;

    if norm(&dist_vec) < &a.radius + &a.confidence + &obstacle.radius {
        u = -&dist_vec - &a.velocity;
        n = normalize(&-dist_vec)
    } else {
        let start_tau = &start / tau;
        let end_tau = &end / tau;
        let closest = closest_point_on_line(&start_tau, &end_tau, &arr1(&[0.0, 0.0]), 0.0, 1.0);
        let w = closest;
        u = &w - &a.velocity;
        n = -normalize(&w);
    }

    return (u, n);
}

/// Calculate the fastest way out of a disk.
fn out_of_disk(
    disk_center: &Array1<f64>,
    disk_r: f64,
    velocity: &Array1<f64>,
) -> (Array1<f64>, Array1<f64>) {
    let rel_vec = velocity - disk_center;
    let w_length = norm(&rel_vec);
    // rotate vector to outside of disk by 10 degrees
    let n = angle2vec(vec2angle(&rel_vec) + 10.0);
    // calculate length of "u" (i.e., the way out of the disk)
    let u = &n * (disk_r - w_length);
    return (u, n);
}

/// Get adjustment velocities (u and n) with respect to another participant.
pub fn get_adjustment_velocities(
    a: &Participant,
    b: &Participant,
    tau: f64,
) -> (Array1<f64>, Array1<f64>) {
    let x = &b.position - &a.position;
    let r = &a.radius + &a.confidence + &b.radius + &b.confidence;
    let v = &a.velocity - &b.velocity;

    // check, if we are currently colliding
    if norm(&x) < r {
        return out_of_disk(&x, r, &v);
    }

    let disk_center = &x / tau;
    let disk_r = r / tau;
    let adjusted_disk_center = &disk_center * (1.0 - (r / norm(&x)).powi(2));
    // check, if we will collide with front disk
    if norm(&v) < norm(&adjusted_disk_center) && norm(&(&v - &adjusted_disk_center)) < r {
        return out_of_disk(&disk_center, disk_r, &v);
    }

    // get angles for relative positions and velocities
    let position_angle = vec2angle(&x);
    let velocity_angle = vec2angle(&v);

    let difference_angle = angle_diff(position_angle, velocity_angle);

    // calculate the angles of the left and right side of the cone
    let side_angle = arcsin(r, norm(&x));
    let right_side_angle = position_angle - side_angle;
    let left_side_angle = position_angle + side_angle;

    // calculate the vectors for left and right side of the cone
    let right_side = angle2vec(right_side_angle);
    let left_side = angle2vec(left_side_angle);

    // calculate closest points on sides of the cone, respectively
    let left_point = closest_point_on_line(&arr1(&[0.0, 0.0]), &left_side, &v, 0.0, 1.0);
    let right_point = closest_point_on_line(&arr1(&[0.0, 0.0]), &right_side, &v, 0.0, 1.0);

    // calculate vectors from current velocity to closest points in left and right sides
    let left_u = &left_point - &v;
    let right_u = &right_point - &v;

    // calculate length of vectors to closest points
    let left_dist = norm(&left_u);
    let right_dist = norm(&right_u);

    // decide, which side is closer and pick that u
    let mut u: Array1<f64>;
    if left_dist < right_dist {
        u = left_u;
    } else {
        u = right_u;
    }

    // "n" is usually just a normalized version of "u".
    // but if we are outside of the cone, we need to use the inverse
    let mut n = normalize(&u);
    if difference_angle > side_angle {
        n *= -1.0;
        u /= 2.0;
    } else {
        u *= 2.0;
    }

    return (u, n);
}

/// Intersect all given halfplanes and find the best velocity
pub fn halfplane_intersection(
    halfplanes_u: &Vec<Halfplane>,
    current_velocity: &Array1<f64>,
    optimal_point: &Array1<f64>,
) -> Option<Array1<f64>> {
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

    return Some(new_point);
}

/// Intersect one halfplane with a set of halfplanes.
fn intersect_halfplane_with_other_halfplanes(
    plane: &Halfplane,
    other_planes: &[Halfplane],
) -> (Option<f64>, Option<f64>) {
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

        if den == 0.0 {
            if num == 0.0 {
                return (None, None);
            }
            continue;
        }

        let t = &num / &den;
        if den > 0.0 {
            right = right.min(t);
        } else {
            left = left.max(t);
        }

        if left > right {
            return (None, None);
        }
    }

    return (Some(left), Some(right));
}
