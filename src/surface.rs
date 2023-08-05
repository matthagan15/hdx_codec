use std::collections::HashMap;

use rand::prelude::*;
use mhgl::PGraph;

struct ErrorTracker {
    pub x: bool,
    pub z: bool,
    error_probability: f64,
}

impl ErrorTracker {
    fn new(error_probability: f64) -> Self {
        ErrorTracker { x: false, z: false, error_probability }
    }
    fn flip_x(&mut self) {
        self.x = true ^ self.x;
    }

    fn flip_z(&mut self) {
        self.z = true ^ self.z;
    }

    fn meas(&self) -> (bool, bool) {
        (self.x, self.z)
    }

    fn sample_errors(&mut self) {
        let mut rng = thread_rng();
        self.x = rng.gen_bool(self.error_probability);
        self.z = rng.gen_bool(self.error_probability);
    }
}

struct SurfaceCode {
    // assume lattice for now
    length: usize,
    height: usize,
    distance: f64, // TODO: might only need usize
    graph: PGraph<u32>,
    lattice_to_node: HashMap<(usize, usize), u32>,
    lattice_index_to_errors: HashMap<(usize, usize), ErrorTracker>,
}

