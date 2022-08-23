use ndarray::{arr1, Array1};

use crate::{
    geometry::{halfplane_intersection, out_of_disk},
    halfplane::Halfplane,
    math::{
        self, angle2vec, angle_diff, arcsin, closest_point_on_line, dist, max, norm, normalize,
        vec2angle,
    },
    obstacle::Obstacle,
};

#[cfg(feature = "participant_obstacles")]
const PART_OBSTACLES: bool = true;

#[cfg(not(feature = "participant_obstacles"))]
const PART_OBSTACLES: bool = false;

const EPSILON: f64 = 0.0;
/// Struct for representing a participant in the ORCA algorithm.
#[derive(Debug, Clone)]
pub struct Participant {
    pub position: Array1<f64>,
    pub velocity: Array1<f64>,
    pub radius: f64,
    pub confidence: f64,
    pub vmax: f64,
    pub target: Array1<f64>,
    pub in_obstacle: bool,
}

impl Participant {
    pub fn update_position(&mut self, position: &Array1<f64>) {
        self.position = position.clone();
        self.velocity = &self.target - &self.position;
        if math::norm(&self.velocity) > self.vmax {
            self.velocity = math::normalize(&self.velocity) * self.vmax
        }
    }

    pub fn in_obstacle(&mut self) {
        self.in_obstacle = true;
    }

    pub fn is_static(&self) -> bool {
        norm(&self.velocity) <= EPSILON
    }

    /// Perform ORCA for a set of other participants and obstacles.
    ///
    /// The return value is the new velocity for the provided participant ("we").
    pub fn orca(
        &self,
        others: &mut [Participant],
        obstacles: &[Obstacle],
        tau: f64,
    ) -> Array1<f64> {
        let (mut halfplanes, obstacle_planes) = self.generate_halfplanes(others, obstacles, tau);

        let mut new_vel: Option<Array1<f64>> = None;
        while new_vel.is_none() {
            // combine halfplanes
            let mut hp = halfplanes.to_vec();
            let mut op = obstacle_planes.to_vec();
            hp.append(&mut op);
            // calculate new velocity
            new_vel = halfplane_intersection(&hp, &self.velocity, &self.velocity);
            // adjust halfplanes (move them outward)
            let mut new_halfplanes: Vec<Halfplane> = Vec::new();
            for l in halfplanes {
                new_halfplanes.push(Halfplane {
                    u: &l.u - &(&l.n * 0.0001),
                    n: l.n.clone(),
                });
            }
            halfplanes = new_halfplanes
        }
        // shorten our new velocity, so we are not faster than our max velocity
        let mut vel = new_vel.unwrap();
        if norm(&vel) > self.vmax {
            vel = normalize(&vel) * self.vmax;
        }
        vel
    }
    /// Generate the halfplanes for other participants and static obstacles.
    ///
    /// Participants, which are too close to each other will be combined into a static obstacle aswell.
    fn generate_halfplanes(
        &self,
        participants: &mut [Participant],
        obstacles: &[Obstacle],
        tau: f64,
    ) -> (Vec<Halfplane>, Vec<Halfplane>) {
        let mut obstacle_planes = Vec::new();
        let mut halfplanes = Vec::new();
        // some cloning, since rust does not like multiple mutable borrows
        let parts = participants.to_vec();
        for (i, p) in parts.iter().enumerate() {
            let mut in_obstacle = false;
            if PART_OBSTACLES && p.is_static() {
                for (j, other) in participants.iter_mut().enumerate() {
                    // check, if "p" is already part of an obstacle
                    if i == j {
                        in_obstacle = other.in_obstacle
                    }
                    // ignore all participants before "other"
                    if j < i + 1 {
                        continue;
                    }
                    // check, if "other" is static _and_ too close
                    if other.is_static()
                        && dist(&p.position, &other.position)
                            < p.radius
                                + other.radius
                                + p.confidence
                                + other.confidence
                                + (self.confidence + self.radius) * 2.0
                    {
                        let (u, n) = self.obstacle_collision(
                            &Obstacle {
                                start: p.position.clone(),
                                end: other.position.clone(),
                                radius: max(
                                    p.radius + p.confidence,
                                    other.radius + other.confidence,
                                ),
                            },
                            tau,
                        );
                        obstacle_planes.push(Halfplane { u, n });
                        in_obstacle = true;
                        other.in_obstacle();
                    }
                }
                // if this participant is not part of a static obstacle, we treat it like a "normal" participant
                if !in_obstacle {
                    let (u, n) = self.get_adjustment_velocities(p, tau);
                    halfplanes.push(Halfplane { u, n });
                }
            } else {
                let (u, n) = self.get_adjustment_velocities(p, tau);
                halfplanes.push(Halfplane { u, n });
            }
        }

        for obstacle in obstacles {
            let (u, n) = self.obstacle_collision(obstacle, tau);
            obstacle_planes.push(Halfplane { u, n });
        }

        (halfplanes, obstacle_planes)
    }

    /// Get adjustment velocities (u and n) with respect to another participant.
    fn get_adjustment_velocities(&self, b: &Participant, tau: f64) -> (Array1<f64>, Array1<f64>) {
        let x = &b.position - &self.position;
        let r = self.radius + self.confidence + b.radius + b.confidence;
        let v = &self.velocity - &b.velocity;

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

        (u, n)
    }

    /// Calculate u and n for a static obstacle.
    fn obstacle_collision(&self, obstacle: &Obstacle, tau: f64) -> (Array1<f64>, Array1<f64>) {
        let r_vec = normalize(&self.velocity) * (self.radius + self.confidence + obstacle.radius);
        let start = &obstacle.start - &(&self.position + &r_vec);
        let end = &obstacle.end - &(&self.position - &r_vec);

        let dist_vec =
            closest_point_on_line(&obstacle.start, &obstacle.end, &self.position, 0.0, 1.0)
                - &self.position;

        let u: Array1<f64>;
        let n: Array1<f64>;

        if norm(&dist_vec) < self.radius + self.confidence + obstacle.radius {
            u = -&dist_vec - &self.velocity;
            n = normalize(&-dist_vec)
        } else {
            let start_tau = &start / tau;
            let end_tau = &end / tau;
            let closest = closest_point_on_line(&start_tau, &end_tau, &arr1(&[0.0, 0.0]), 0.0, 1.0);
            let w = closest;
            u = &w - &self.velocity;
            n = -normalize(&w);
        }

        (u, n)
    }
}
