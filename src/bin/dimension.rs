extern crate itertools;

use itertools::Itertools;
use temperley_lieb_cat::*;
use tl_jantzen::mod_inv;

use std::env;

/// Calculate the Gram matrix for a standard module
///
/// Calculates explicitly the Gram matrix (in the diagram
/// basis) of the standard module S(n, m) with δ = 2.
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
                gm[i][j] = (2 as i64).pow(a as u32);
            }
        }
    }
    return gm;
}

/// Calculate the row echelon form of a matrix
///
/// Calculates the unreduced row echelon form of a matrix 
/// while performing the same operations on the identity matrix
/// with the same dimensions. Returns both matricies. Division is 
/// performed modulo the argument p
pub fn row_echelon_form(matrix: &Vec<Vec<i64>>, p: i64) -> (Vec<Vec<i64>>, Vec<Vec<i64>>, usize) {
    let mut matrix_out: Vec<Vec<i64>> = matrix.to_vec();
    let mut pivot = 0;
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();
    let mut rank =0;

    let mut identity_matrix: Vec<Vec<i64>> =
        vec![vec![0; column_count as usize]; row_count as usize];
    for k in 0..row_count {
        identity_matrix[k][k] = 1;
    }
    'rowLoop: for r in 0..row_count {
       
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
                    break 'rowLoop;
                }
            }
        }

        matrix_out.swap(r, i);
        identity_matrix.swap(r, i);

        let q = matrix_out[r][pivot];

        let mod_inverse = mod_inv(((q % p) + p) % p, p);

        for j in r + 1..row_count {
            let hold = matrix_out[j][pivot];
            for k in 0..column_count {
                matrix_out[j][k] = matrix_out[j][k] - (hold * matrix_out[r][k] * mod_inverse);
                identity_matrix[j][k] =
                    identity_matrix[j][k] - (hold * identity_matrix[r][k] * mod_inverse);
            }
        }
         rank = r;
        pivot = pivot + 1;

    }
    
    return (matrix_out, identity_matrix, rank+1);
}

/// Reduces a matrix modulo p
///
/// Reduces each entry in a matrix mod p. Returns the corresponding matrix
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

/// This fn needs to be modified...
fn reform_inner_product(new_basis: &[Vec<i64>], old_gram: &Vec<Vec<i64>>) -> Vec<Vec<i64>> {
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

/// This fn needs to be modified...
fn recursive_ops(m: usize, n: usize, p: i64) {
    let mut g = gram_matrix(m, n);

    loop {
        let mut h = reduce_mod_p(&mut g, p);
        let (reduced_matrix, transform_matrix, rank) = row_echelon_form(&mut h, p);

        println!("dimension of head is: {}", &rank);
        if rank == reduced_matrix.len() {
            break;
        }
        let basis_of_rad = &transform_matrix[(rank as usize)..];
        g = reform_inner_product(basis_of_rad, &mut g);
        for i in 0..g.len() {
            for j in 0..g.len() {
                g[i][j] = g[i][j] / p;
            }
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        println!("Enter m n p seperated by spaces");
        return;
    }
    let n: usize = args[1].parse().expect("Please enter an integer for n");
    let m: usize = args[2].parse().expect("Please enter an integer for m");
    let p: i64 = args[3].parse().expect("Please enter an integer for p");

    recursive_ops(n, m, p);
}
#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_num_zero_rows() {
        let m = vec![
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            vec![0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            vec![0, -1, 0, 1, 0, 0, 0, 0, 0, 0],
            vec![-1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, -1, 0, 1, 0, 0, 0],
            vec![0, 1, 0, -1, 0, 0, 0, 1, 0, 0],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            vec![0, -1, 0, 1, 0, 0, 0, -1, 0, 1],
        ];
        assert_eq!(num_zero_rows(&m), 2);
    }

    #[test]
    fn test_reduce_mod_p() {
        let m = vec![
            vec![2, 0, -13, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, -1, 0, 1],
        ];

        let n = vec![[2, 0, 2, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 2, 0, 1]];
        assert_eq!(reduce_mod_p(&m, 3), n);
    }

    #[test]
    fn test_row_echelon_form() {
        let m = vec![
            vec![2, 0, -13, 0],
            vec![4, 1, 0, 0],
            vec![0, 8, 1, 0],
            vec![5, -1, 0, 1],
        ];
        let ans = (
            vec![
                vec![2, 0, 2, 0],
                vec![0, 1, 2, 0],
                vec![0, 0, 0, 1],
                vec![0, 0, 0, 0],
            ],
            vec![
                vec![1, 0, 0, 0],
                vec![-2, 1, 0, 0],
                vec![0, -2, 0, 1],
                vec![4, -2, 1, 0],
            ],
        );
        assert_eq!(row_echelon_form(&m, 3), ans);
    }

    #[test]
    fn test_gram_matrix() {
        let m: Vec<Vec<i64>> = vec![
            vec![4, 2, 0, 0, 2, 1, 0, 0, 2, 0, 0, 0, 0, 0],
            vec![2, 4, 2, 0, 1, 2, 1, 0, 1, 0, 0, 0, 0, 0],
            vec![0, 2, 4, 2, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0],
            vec![0, 0, 2, 4, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0],
            vec![2, 1, 0, 0, 4, 2, 0, 0, 1, 0, 0, 1, 0, 0],
            vec![1, 2, 1, 0, 2, 4, 2, 0, 2, 1, 0, 2, 0, 0],
            vec![0, 1, 2, 1, 0, 2, 4, 2, 1, 2, 1, 1, 0, 0],
            vec![0, 0, 1, 2, 0, 0, 2, 4, 0, 1, 2, 0, 0, 0],
            vec![2, 1, 0, 0, 1, 2, 1, 0, 4, 2, 0, 1, 0, 1],
            vec![0, 0, 0, 0, 0, 1, 2, 1, 2, 4, 2, 2, 1, 2],
            vec![0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 4, 1, 2, 1],
            vec![0, 0, 0, 0, 1, 2, 1, 0, 1, 2, 1, 4, 2, 1],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 4, 2],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 4],
        ];
        assert_eq!(gram_matrix(7, 3), m);
    }
}
