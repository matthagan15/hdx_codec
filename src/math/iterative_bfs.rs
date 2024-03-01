use std::{
    collections::{HashMap, HashSet, VecDeque}, hash::Hash, io::Write, path::PathBuf, rc::Rc, time::Instant
};

use mhgl::HGraph;

use crate::math::coset_complex::*;


use super::polymatrix::PolyMatrix;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CosetRep {
    rep: PolyMatrix,
    type_ix: usize,
}

impl CosetRep {
    pub fn new(rep: &PolyMatrix, type_ix: usize, subgroups: &Rc<Subgroups>) -> Self {
        let coset = subgroups.get_coset(rep, type_ix);
        Self {
            type_ix,
            rep: coset[0].clone(),
        }
    }
}


#[derive(Debug, Clone)]
pub struct TriangleRep {
    rep_zero: PolyMatrix,
    rep_one: PolyMatrix,
    rep_two: PolyMatrix,
    pub distance_from_origin: usize,
    subgroups: Rc<Subgroups>,
}

impl TriangleRep {
    pub fn new(rep: &PolyMatrix, distance: usize, subgroups: Rc<Subgroups>) -> Self {
        let r0 = subgroups.get_coset_rep(rep, 0);
        let r1 = subgroups.get_coset_rep(rep, 1);
        let r2 = subgroups.get_coset_rep(rep, 2);
        Self {
            rep_zero: r0,
            rep_one: r1,
            rep_two: r2,
            distance_from_origin: distance,
            subgroups: subgroups.clone(),
        }
    }

    pub fn complete_star(&self, subgroups: &Rc<Subgroups>) -> Vec<Self> {
        let c0 = CosetRep::new(&self.rep_zero, 0, subgroups);
        let mut star_zero = star_of_vertex_rep_final(&c0, subgroups, self.distance_from_origin);

        let c1 = CosetRep::new(&self.rep_one, 1, subgroups);
        let mut star_one = star_of_vertex_rep_final(&c1, subgroups, self.distance_from_origin);

        let c2 = CosetRep::new(&self.rep_two, 2, subgroups);
        let mut star_two = star_of_vertex_rep_final(&c2, subgroups, self.distance_from_origin);

        let mut ret = Vec::with_capacity(star_zero.len() + star_one.len() + star_two.len());
        ret.append(&mut star_zero);
        ret.append(&mut star_one);
        ret.append(&mut star_two);
        ret
    }

    /// Returns the cosets (which store the canonical aka sorted representative of the coset)
    ///
    pub fn get_cosets(&self) -> (CosetRep, CosetRep, CosetRep) {
        let c0 = CosetRep::new(&self.rep_zero, 0, &self.subgroups);
        let c1 = CosetRep::new(&self.rep_one, 1, &self.subgroups);
        let c2 = CosetRep::new(&self.rep_two, 2, &self.subgroups);
        (c0, c1, c2)
    }
}

impl PartialEq for TriangleRep {
    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }

    fn eq(&self, other: &Self) -> bool {
        let sub_zero = self.subgroups.get_coset(&self.rep_zero, 0);
        let sub_1 = self.subgroups.get_coset(&self.rep_one, 1);
        let sub_2 = self.subgroups.get_coset(&self.rep_two, 2);
        let rhs_zero = self.subgroups.get_coset(&other.rep_zero, 0);
        let rhs_1 = self.subgroups.get_coset(&other.rep_one, 1);
        let rhs_2 = self.subgroups.get_coset(&other.rep_two, 2);
        sub_zero == rhs_zero && sub_1 == rhs_1 && sub_2 == rhs_2
    }
}

impl Eq for TriangleRep {
}

impl Hash for TriangleRep {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.rep_zero.hash(state);
        self.rep_one.hash(state);
        self.rep_two.hash(state);
    }
}

// pub fn star_of_vertex_rep(coset: &Coset, subgroups: &Rc<Subgroups>, distance_of_vertex: usize) -> Vec<TriangleRep> {
//     let mut t_reps = Vec::new();
//     for rep in subgroups.get_coset(&coset.rep, coset.type_ix) {
//         let t_rep = TriangleRep::new(&rep, distance_of_vertex + 1, subgroups.clone());
//         t_reps.push(t_rep);
//     }
//     t_reps
// }

pub fn star_of_vertex_rep_final(coset: &CosetRep, subgroups: &Rc<Subgroups>, distance_of_vertex: usize) -> Vec<TriangleRep> {
    let mut t_reps = Vec::new();
    for rep in subgroups.get_coset(&coset.rep, coset.type_ix) {
        let t_rep = TriangleRep::new(&rep, distance_of_vertex + 1, subgroups.clone());
        t_reps.push(t_rep);
    }
    t_reps
}


/// Flushes all triangles from visited if their distance is strictly less than 
/// the provided `flushing_upper_limit`
fn flush_visited(visited: &mut HashSet<TriangleRep>, flushing_upper_limit: usize, coset_to_node: &mut HashMap<CosetRep, u32>) {
    let mut to_flush = Vec::new();
    let mut to_save = HashSet::new();
    for triangle in visited.drain() {
        if triangle.distance_from_origin <= flushing_upper_limit {
            to_flush.push(triangle);
        } else {
            to_save.insert(triangle);
        }
    }
    for t in to_save.into_iter() {
        visited.insert(t);
    }
    // let file_path = "/Users/matt/repos/qec/tmp/triangle_toilet.csv";
    // let mut file = std::fs::File::create(file_path).expect("No toilet??");
    for t in to_flush {
        // let s = serde_json::to_string(&t).expect("cannot serialize triangle.");
        // file.write_all(s.as_bytes()).expect("cannot write triangle bytes");
        // file.write_all(",".as_bytes()).expect("cannot write comma");
        let (c0, c1, c2) = t.get_cosets();
        coset_to_node.remove(&c0);
        coset_to_node.remove(&c1);
        coset_to_node.remove(&c2);
    }
}

// TODO: Need to figure out how to check that I'm actually exploring
// the entire group. That is a difficult test to run, because I don't
// think the group can fit on RAM. 

/// Currently computes the entire group using Breadth-First-Search
/// starting at the identity matrix over the generators provided.
fn triangle_based_bfs(subgroups: Subgroups, verbose: bool) -> HGraph {
    if verbose {
        println!("Computing group with the following parameters:");
        println!("quotient: {:}", subgroups.quotient);
        println!("dim: {:}", subgroups.dim);
        for (t, set) in subgroups.type_to_generators.iter() {
            println!("{:}", "*".repeat(50));
            println!("Type: {:}", t);
            for s in set.iter() {
                println!("{:}", s);
            }
        }
    }
    let dim = subgroups.dim as u32;
    let p = subgroups.quotient.field_mod;
    let deg = subgroups.quotient.degree() as u32;
    let total_expected_matrices = (p.pow(deg)).pow((dim * dim) - 1) ;
    let subgroups = Rc::new(subgroups);
    let mut hg = HGraph::new();
    let mut coset_to_node: HashMap<CosetRep, u32> = HashMap::new();
    let e = PolyMatrix::id(subgroups.dim, subgroups.quotient.clone());
    let starting_triangle = TriangleRep::new(&e, 0, subgroups.clone());

    // TODO: Currently a matrix is being stored twice while it is in
    // the frontier as we also put it in visited. Instead just keep track
    // of an extra bit if the matrix is visited or not?
    let mut completed = 0_u32;
    let start_time = Instant::now();
    let mut frontier = VecDeque::from([starting_triangle.clone()]);
    let mut visited = HashSet::from([starting_triangle.clone()]);
    
    
    let mut counter = 0;
    let mut last_flushed_distance = 0;
    while frontier.len() > 0 {
        counter += 1;
        if (counter % 100) == 0 {
            let num_polys = subgroups.quotient.field_mod.pow(subgroups.quotient.degree() as u32);
            println!("{:}", ".".repeat(50));
            println!(
                "frontier length: {:} , {:.4}% of visitied.",
                frontier.len(),
                frontier.len() as f64 / visited.len() as f64
            );
            println!("completed length: {:}", completed);
            println!("visited length: {:}", visited.len());
            let tot_time = start_time.elapsed();
            let ns_per_tri = tot_time.as_nanos() / completed as u128;
            let time_per_triangle = tot_time / completed;
            println!("time per completed triangle: {:} ns", ns_per_tri);
            let predicted_time_left = time_per_triangle * (total_expected_matrices  - completed);
            let percent_rem = 1. - (completed as f64 / total_expected_matrices as f64);
            println!("predicted seconds left: {:}", predicted_time_left.as_secs());
            println!("percent remaining: {:}", percent_rem);
        }
        let x = frontier.pop_front().expect("no frontier?");
        let mut new_edge_nodes = Vec::new();
        let (r0, r1, r2) = x.get_cosets();
        if coset_to_node.contains_key(&r0) == false {
            let new_node = hg.add_nodes(1)[0];
            new_edge_nodes.push(new_node);
            coset_to_node.insert(r0.clone(), new_node);
        } else {
            new_edge_nodes.push(*coset_to_node.get(&r0).unwrap());
        }
        if coset_to_node.contains_key(&r1) == false {
            let new_node = hg.add_nodes(1)[0];
            new_edge_nodes.push(new_node);
            coset_to_node.insert(r1.clone(), new_node);
        } else {
            new_edge_nodes.push(*coset_to_node.get(&r1).unwrap());
        }
        if coset_to_node.contains_key(&r2) == false {
            let new_node = hg.add_nodes(1)[0];
            new_edge_nodes.push(new_node);
            coset_to_node.insert(r2.clone(), new_node);
        } else {
            new_edge_nodes.push(*coset_to_node.get(&r2).unwrap());
        }

        hg.create_edge_no_dups(&new_edge_nodes[..]);

        if x.distance_from_origin > last_flushed_distance + 1 {
            if verbose {
                println!("{:}", "#".repeat(55));
                println!("Flushing shit");
                flush_visited(&mut visited, last_flushed_distance + 1, &mut coset_to_node)
            }
            last_flushed_distance += 1;
        }
        let x_neighbors = x.complete_star(&subgroups);
        for new_neighbor in x_neighbors {
            if visited.contains(&new_neighbor) == false {
                visited.insert(new_neighbor.clone());
                frontier.push_back(new_neighbor);
            }
        }
        completed += 1;
    }
    hg
}


mod tests {
    use std::{io::Write, path::PathBuf, rc::Rc, str::FromStr};

    use crate::math::{polymatrix::PolyMatrix, polynomial::FiniteFieldPolynomial};

    use super::{triangle_based_bfs, Subgroups, TriangleRep};

    fn simple_quotient_and_field() -> (u32, FiniteFieldPolynomial) {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        (p, q)
    }
    #[test]
    fn test_triangle_rep_star() {
        let (p, q) = simple_quotient_and_field();
        let sg = Rc::new(Subgroups::new(3, &q));
        let e = PolyMatrix::id(3, q.clone());
        let t = TriangleRep::new(&e, 0, sg.clone());
        let star = t.complete_star(&sg);
        println!("subgroups:");
        for (t, v) in sg.type_to_generators.iter() {
            println!("type: {:}", t);
            println!("gens:");
            for g in v {
                println!("{:}", g);
            }
        }
        for rep in star {
            println!("type 0 rep");
            println!("{:}", rep.rep_zero);
            println!("type 1 rep");
            println!("{:}", rep.rep_one);
            println!("type 2 rep");
            println!("{:}", rep.rep_two);
        }
    }

    #[test]
    fn test_triangle_bfs() {
        let p =3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let cg = Subgroups::new(3, &q);
        let hg = triangle_based_bfs(cg, true);

        let hgraph_file_path = PathBuf::from_str("/Users/matt/repos/qec/tmp/triangle_bfs.hg").expect("No pathbuf?");
        let mut hgraph_file =
            std::fs::File::create(hgraph_file_path).expect("Could not create hgraph file.");
        let hgraph_string =
            serde_json::to_string(&hg).expect("Could not serialize hgraph.");
        hgraph_file
            .write(hgraph_string.as_bytes())
            .expect("Could not write hgraph to file.");
    }
}