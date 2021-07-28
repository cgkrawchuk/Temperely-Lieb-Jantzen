extern crate itertools;

use itertools::Itertools;
use temperley_lieb_cat::*;
use tl_jantzen::{mod_inv, Matrix};

use std::env;

/// Calculate the Gram matrix for a standard module
///
/// Calculates explicitly the Gram matrix (in the diagram
/// basis) of the standard module S(n, m) with δ = 2.
pub fn gram_matrix(n: usize, m: usize) -> Matrix {
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
    let mut gm = Matrix::new(x, x);
    for i in 0..x {
        for j in 0..x {
            let (a, b) = (monic_diagrams[i].involute()) * monic_diagrams[j];
            if b == TLDiagram::id(m) {
                gm[(i, j)] = (2_i64).pow(a as u32);
            }
        }
    }
    gm
}

/// Calculate the row echelon form of a matrix
///
/// Calculates the unreduced row echelon form of a matrix
/// while performing the same operations on the identity matrix
/// with the same dimensions. Returns both matricies. Division is
/// performed modulo the argument p
pub fn row_echelon_form(matrix: &Matrix, p: i64) -> (Matrix, Matrix, usize) {
    let mut matrix_out = matrix.clone();
    let mut pivot_row = 0;
    let mut rank = 0;

    let mut identity_matrix = Matrix::identity(matrix.rows);

    'col_loop: for column in 0..matrix.cols {
        reduce_mod_p(&mut matrix_out, p);

        let mut i = pivot_row;
        while matrix_out[(i, column)] == 0 {
            i += 1;
            if i == matrix.rows {
                continue 'col_loop;
            }
        }

        matrix_out.swap_rows(pivot_row, i);
        identity_matrix.swap_rows(pivot_row, i);

        let q = matrix_out[(pivot_row, column)];

        let mod_inverse = mod_inv(q, p);

        for j in pivot_row + 1..matrix.rows {
            let hold = matrix_out[(j, column)];
            for k in 0..matrix.cols {
                matrix_out[(j, k)] -= hold * matrix_out[(pivot_row, k)] * mod_inverse;
            }
            for k in 0..identity_matrix.cols {
                identity_matrix[(j, k)] -= hold * identity_matrix[(pivot_row, k)] * mod_inverse;
            }
        }
        pivot_row += 1;
        if pivot_row == matrix.rows {
            break;
        }
    }
    reduce_mod_p(&mut identity_matrix, p);
    (matrix_out, identity_matrix, pivot_row)
}

}

/// Reduces a matrix modulo p
///
/// Reduces each entry in a matrix mod p in-place.
pub fn reduce_mod_p(matrix: &mut Matrix, p: i64) {
    for i in 0..matrix.rows {
        for j in 0..matrix.cols {
            matrix[(i, j)] = ((matrix[(i, j)] % p) + p) % p;
        }
    }
}

/// This fn needs to be modified...
fn reform_inner_product(new_basis: &[&[i64]], old_gram: &Matrix) -> Matrix {
    let mut matrix_out = Matrix::new(old_gram.rows, new_basis.len());
    for i in 0..old_gram.rows {
        for j in 0..new_basis.len() {
            let mut inner = 0;
            for k in 0..old_gram.cols {
                inner += old_gram[(i, k)] * new_basis[j][k];
            }

            matrix_out[(i, j)] = inner;
        }
    }

    matrix_out
}

/// This fn needs to be modified...
fn recursive_ops(m: usize, n: usize, p: i64) {
    let mut g = gram_matrix(m, n);

    loop {
        let mut h = g.clone();
        reduce_mod_p(&mut h, p);
        let (reduced_matrix, transform_matrix, rank) = row_echelon_form(&h, p);

        println!("reduced matrix: {}", &reduced_matrix);
        println!("dimension of head is: {}", &rank);

        if rank == reduced_matrix.cols {
            break;
        }
        let rows = transform_matrix.row_list();
        let basis_of_rad = &rows[rank..];
        g = reform_inner_product(basis_of_rad, &g);
        for i in 0..g.rows {
            for j in 0..g.cols {
                g[(i, j)] /= p;
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
    fn test_reduce_mod_p() {
        let mut m: Matrix = vec![
            vec![2, 0, -13, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, -1, 0, 1],
        ]
        .into();

        let n: Matrix = vec![
            vec![2, 0, 2, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 2, 0, 1],
        ]
        .into();
<<<<<<< HEAD
        assert_eq!(reduce_mod_p(&m, 3), n);
    }


=======
        reduce_mod_p(&mut m, 3);
        assert_eq!(m, n);
    }
>>>>>>> 2aef112bd3e16f9143d0b88b7b8bf14f2b79f662

    #[test]
    fn test_row_echelon_form() {
        let m1: Matrix = vec![
            vec![2, 0, -13, 0],
            vec![4, 1, 0, 0],
            vec![0, 8, 1, 0],
            vec![5, -1, 0, 1],
        ]
        .into();
        let ans1 = (
            vec![
                vec![2, 0, 2, 0],
                vec![0, 1, 2, 0],
                vec![0, 0, 0, 1],
                vec![0, 0, 0, 0],
            ]
            .into(),
            vec![
                vec![1, 0, 0, 0],
                vec![-2, 1, 0, 0],
                vec![0, -2, 0, 1],
                vec![4, -2, 1, 0],
            ]
            .into(),
            3,
        );
        assert_eq!(row_echelon_form(&m1, 3), ans1);
        let m2: Matrix = vec![vec![4, 0, -3, 0], vec![4, 1, 0, 5], vec![0, 8, 1, 0]].into();
        let ans2 = (
            vec![vec![1, 0, 0, 0], vec![0, 1, 0, 2], vec![0, 0, 1, 2]].into(),
            vec![
                vec![1, 0, 0, 0],
                vec![-1, 1, 0, 0],
                vec![2, -2, 1, 0],
                vec![0, 0, 0, 1],
            ]
            .into(),
            3,
        );
        assert_eq!(row_echelon_form(&m2, 3), ans2);
    }

    #[test]
    fn test_gram_matrix() {
        let m: Matrix = vec![
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
        ]
        .into();
        assert_eq!(gram_matrix(7, 3), m);
    }
}
