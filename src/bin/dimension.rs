extern crate itertools;

use itertools::Itertools;
use std::io::{self, BufRead};
use temperley_lieb_cat::*;

pub fn gram_matrix(n: usize, m: usize) -> Vec<Vec<i64>> {
    fn ok(tab: &Vec<usize>) -> bool {
        for (c, i) in tab.iter().enumerate() {
            if *i < 2 * (c + 1) {
                return false;
            }
        }
        true
    }
    let mut monic_diagrams = Vec::new();
    let prop = m;
    let rn = (n - prop) / 2;
    for lt in (1..(n + 1)).combinations(rn).filter(ok) {
        monic_diagrams.push(TLDiagram::new(
            n,
            m,
            lt.iter().fold(0, |l, n| l | (1 << n)),
            0,
        ));
    }

    let x = monic_diagrams.len();
    let mut gm: Vec<Vec<i64>> = vec![vec![0; x as usize]; x as usize];
    for i in 0..x {
        for j in 0..x {
            let (a, b) = &(monic_diagrams[i].involute()) * &monic_diagrams[j];
            if b == TLDiagram::id(m) {
                if a == 0 {
                    gm[i][j] = 1
                } else {
                    gm[i][j] = (2 as i64).pow(a as u32);
                }
            } else {
                gm[i][j] = 0;
            }
        }
    }
    return gm;
}

pub fn row_echelon_form(matrix: & Vec<Vec<i64>>, p: i64) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let mut matrix_out: Vec<Vec<i64>> = matrix.to_vec();
    let mut pivot = 0;
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();

    let mut identity_matrix: Vec<Vec<i64>> =
        vec![vec![0; column_count as usize]; row_count as usize];
    for k in 0..row_count {
        //hopefully row_count is less than column_count
        identity_matrix[k][k] = 1;
    }
    //println!("{:?}", identity_matrix);
    for r in 0..row_count {
        matrix_out = reduce_mod_p(&mut matrix_out, p);
        if column_count <= pivot {
            break;
        }
        let mut i = r;
        while matrix_out[i][pivot] == 0 {
            i = i + 1;
            if i == row_count {
                i = r;
                pivot = pivot + 1;
                if column_count == pivot {
                    pivot = pivot - 1;
                    break;
                }
            }
        }
         matrix_out.swap(r, i);
         identity_matrix.swap(r, i);

        let q = matrix_out[r][pivot];
        let a = ((q % p) + p) % p;
        if a != 0 {
            let mod_inverse = mod_inv(((a % p) + p) % p, p);
            for j in 0..column_count {
                matrix_out[r][j] = matrix_out[r][j] * mod_inverse;
                identity_matrix[r][j] = identity_matrix[r][j] * mod_inverse;
            }
        }
        for j in 0..row_count {
            if j != r {
                let hold = matrix_out[j][pivot];
                for k in 0..column_count {
                    matrix_out[j][k] = matrix_out[j][k] - (hold * matrix_out[r][k]);
                    identity_matrix[j][k] = identity_matrix[j][k] - (hold * identity_matrix[r][k]);
                }
            }
        }
        pivot = pivot + 1;
    }
    (matrix_out, identity_matrix)
}


pub fn reduce_mod_p(matrix: &Vec<Vec<i64>>, p: i64) -> Vec<Vec<i64>> {
    let mut matrix_out: Vec<Vec<i64>> = matrix.to_vec();
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();
    for i in 0..row_count {
        for j in 0..column_count {
            matrix_out[i][j] = ((matrix_out[i][j] % p) + p) % p;
        }
    }
    return matrix_out;
}

fn mod_inv(a: i64, module: i64) -> i64 {
    assert!(a % module != 0, "number is 0 mod...");
    let mut mn = (module, a);
    let mut xy = (0, 1);

    while mn.1 != 0 {
        xy = (xy.1, xy.0 - (mn.0 / mn.1) * xy.1);
        mn = (mn.1, mn.0 % mn.1);
    }

    while xy.0 < 0 {
        xy.0 += module;
    }
    xy.0
}

pub fn num_zero_rows(matrix: &Vec<Vec<i64>>) -> i64 {
    let mut num_zeroes = 0; 
    let zero_vec = vec![0; matrix[0].len()]; 
    for x in matrix.iter(){
        if x == &zero_vec {  num_zeroes+=1;}
    }
    return num_zeroes;
}

fn reform_inner_product(new_basis: &Vec<Vec<i64>>, old_gram: &Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    let mut matrix_out: Vec<Vec<i64>> = vec![vec![0; new_basis.len()]; new_basis.len()];
    for i in 0..new_basis.len() {
        for j in 0..new_basis.len() {
            let mut inner = 0;
            for k in 0..old_gram.len() {
                for l in 0..old_gram[0].len() {
                    inner += old_gram[k][l] * new_basis[i][k] * new_basis[j][l];
                }
            }
            matrix_out[i][j] = inner;
        }
    }

    matrix_out
}

fn recursive_ops(m: usize, n: usize, p: i64) {
    let mut g = gram_matrix(m, n);
    //println!("gram Matrix of {0}, {1} is :", m, n);
    //println!("{:?}", &g);

    let mut h = reduce_mod_p(&mut g, p);
    //println!("Reduced mod {}", p);
    //println!("{:?}", &h);

    let (mut reduced_matrix, mut transform_matrix) = row_echelon_form(&mut h, p);

    //println!("Row echelon form: \n{:?}", reduce_mod_p(&reduced_matrix, p));
    //println!(
    //     "transformation matrix:\n{:?}",
    //    reduce_mod_p(&transform_matrix, p)
    // );
    let mut nzr = num_zero_rows(&reduced_matrix);
    let mut rank = (reduced_matrix.len() as i64) - nzr;
    println!("dimension of head is: {}", &rank);

    while nzr != 0 as i64 {
        let basis_of_rad = &transform_matrix[(rank as usize)..];
        let mut bor: Vec<Vec<i64>> = Vec::from(basis_of_rad);
        // println!("new basis: {:?}", &basis_of_rad);
        g = reform_inner_product(&mut bor, &mut g);
        //println!("new gram Matrix: {:?}", &g);
        for i in 0..g.len() {
            for j in 0..g.len() {
                g[i][j] = g[i][j] / p;
            }
        }
        //println!("divide by p: \n{:?}", &g);

        h = reduce_mod_p(&mut g, p);
        //println!("Reduced mod {}", p);
        // println!("{:?}", &h);

        let (temp1, temp2) = row_echelon_form(&mut h, p);
        reduced_matrix = temp1;
        transform_matrix = temp2;

        println!("Row echelon form: \n{:?}", reduce_mod_p(&reduced_matrix, p));
        println!(
            "transformation matrix:\n{:?}",
            reduce_mod_p(&transform_matrix, p)
        );
        nzr = num_zero_rows(&reduced_matrix);
        rank = (reduced_matrix.len() as i64) - nzr;
        println!("dimension of head is: {}", rank);
    }
}

fn main() {
    println!("Enter m n p seperated by spaces");
    let reader = io::stdin();
    let numbers: Vec<usize> = reader
        .lock()
        .lines()
        .next()
        .unwrap()
        .unwrap()
        .split(' ')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();

    recursive_ops(numbers[0], numbers[1], numbers[2] as i64);
}
