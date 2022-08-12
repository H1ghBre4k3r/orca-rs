mod geometry;
mod halfplane;
mod math;
pub mod obstacle;
pub mod participant;

use geometry::{get_adjustment_velocities, halfplane_intersection, obstacle_collision};
use halfplane::Halfplane;
use math::{dist, max, norm, normalize};
use ndarray::Array1;
use obstacle::Obstacle;
use participant::Participant;

#[cfg(feature = "participant_obstacles")]
const PART_OBSTACLES: bool = true;

#[cfg(not(feature = "participant_obstacles"))]
const PART_OBSTACLES: bool = false;

const EPSILON: f64 = 0.0;

/// Perform ORCA for a set of other participants and obstacles.
///
/// The return value is the new velocity for the provided participant ("we").
pub fn orca(
    we: &Participant,
    participants: &mut [Participant],
    obstacles: &[Obstacle],
    tau: f64,
) -> Array1<f64> {
    let (mut halfplanes, obstacle_planes) = generate_halfplanes(we, participants, obstacles, tau);

    let mut new_vel: Option<Array1<f64>> = None;
    while new_vel.is_none() {
        // combine halfplanes
        let mut hp = halfplanes.to_vec();
        let mut op = obstacle_planes.to_vec();
        hp.append(&mut op);
        // calculate new velocity
        new_vel = halfplane_intersection(&hp, &we.velocity, &we.velocity);
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
    if norm(&vel) > we.vmax {
        vel = normalize(&vel) * we.vmax;
    }
    vel
}

fn is_static(p: &Participant) -> bool {
    norm(&p.velocity) <= EPSILON
}

/// Generate the halfplanes for other participants and static obstacles.
///
/// Participants, which are too close to each other will be combined into a static obstacle aswell.
fn generate_halfplanes(
    we: &Participant,
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
        if PART_OBSTACLES && is_static(p) {
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
                if is_static(other)
                    && dist(&p.position, &other.position)
                        < p.radius
                            + other.radius
                            + p.confidence
                            + other.confidence
                            + (we.confidence + we.radius) * 2.0
                {
                    let (u, n) = obstacle_collision(
                        we,
                        &Obstacle {
                            start: p.position.clone(),
                            end: other.position.clone(),
                            radius: max(p.radius + p.confidence, other.radius + other.confidence),
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
                let (u, n) = get_adjustment_velocities(we, p, tau);
                halfplanes.push(Halfplane { u, n });
            }
        } else {
            let (u, n) = get_adjustment_velocities(we, p, tau);
            halfplanes.push(Halfplane { u, n });
        }
    }

    for o in obstacles {
        let (u, n) = obstacle_collision(we, o, tau);
        obstacle_planes.push(Halfplane { u, n });
    }

    (halfplanes, obstacle_planes)
}
