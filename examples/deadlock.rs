use orca_rs::{ndarray::arr1, obstacle::Obstacle, participant::Participant};

fn main() {
    env_logger::init();
    let position = arr1(&[1.171660304069519, 0.22933036088943481]);
    let mut we = Participant::new(position.clone(), arr1(&[0.0, 0.0]), 0.15, 0.0)
        .with_inner_state(0.2, arr1(&[0.2, 0.2]));
    we.update_position(&position);

    let mut others = vec![
        Participant::new(
            arr1(&[0.876857340335846, 0.8371948003768921]),
            arr1(&[0.0, 0.0]),
            0.15,
            0.0,
        ),
        Participant::new(
            arr1(&[1.1547585725784302, 0.3865005373954773]),
            arr1(&[0.0, 0.0]),
            0.15,
            0.0,
        ),
    ];
    println!(
        "{}",
        we.orca(
            &mut others,
            //
            &generate_bounding_obstacles(1.6, 2.0),
            // &[],
            2.0
        )
    );
}

pub fn generate_bounding_obstacles(width: f64, height: f64) -> [Obstacle; 2] {
    [
        // Obstacle {
        //     start: arr1(&[0.0, 0.0]),
        //     end: arr1(&[0.0, height]),
        //     radius: 0.01,
        // },
        // Obstacle {
        //     start: arr1(&[0.0, height]),
        //     end: arr1(&[width, height]),
        //     radius: 0.01,
        // },
        // these two obstacles produce "false" half planes
        Obstacle {
            start: arr1(&[width, height]),
            end: arr1(&[width, 0.0]),
            radius: 0.01,
        },
        Obstacle {
            start: arr1(&[width, 0.0]),
            end: arr1(&[0.0, 0.0]),
            radius: 0.01,
        },
    ]
}
