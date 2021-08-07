use crate::matrix::*;
use crate::util::*;

extern crate itertools;

use crate::{binom, extended_euclid, Matrix};
use itertools::Itertools;
use std::ops::Range;
use temperley_lieb_cat::*;

use std::io::{self, BufRead};

use std::cmp::min;

///Bubble sorts columns by their maximum entry
pub fn sort_cols_max(A: &mut Matrix) -> (Matrix) {
    let mut A = A.clone();
    for i in 0..A.cols {
        for j in 0..A.cols - i - 1 {
            if A.col(j).iter().max() > A.col(j + 1).iter().max() {
                A.swap_cols(j, j + 1);
            }
        }
    }
    return A;
}

/// Calculates the square of the Euclidean norm of a vector
pub fn euclid_norm_squared(v: &[i64]) -> i64 {
    v.iter().map(|x| x * x).sum()
}

///Returns the maximum value of a matrix
pub fn max_norm(A: &Matrix) -> i64 {
    return *A.entries().iter().max().unwrap();
}

///Reduces matrix coefficients using normcol method
pub fn normcol(B: &mut Matrix) -> Matrix {
    let mut A = B.clone();
    sort_cols_max(&mut A);

    let mut norm_s = max_norm(&A);
    let mut norm_a = norm_s;
    while norm_s <= norm_a {
        norm_a = norm_s;
        for k in (1..A.cols).rev() {
            if A[(0, k)].signum() == A[(0, k - 1)].signum() {
                A.add_multiple_coli_to_colj(-1, k - 1, k);
            } else {
                A.add_multiple_coli_to_colj(1, k - 1, k);
            }
        }
        A = sort_cols_max(&mut A);

        norm_s = max_norm(&A);
    }

    return A;
}

///lexicographically sorts the columns of a matrix
pub fn sort_cols_lex(B: &mut Matrix, frm:usize) -> Matrix {
    let mut A = B.clone();
    for i in frm..A.cols {
        for j in frm..(A.cols - i - 1) {
            if A.col(j) < A.col(j + 1) {
                A.swap_cols(j, j + 1);
            }
        }
    }

    return A;
}

///Reduces the entries of a matrix using the rosser method
pub fn rosser(B: &mut Matrix,row: usize) -> Matrix {
    let mut A = B.clone();
    for k in (row..A.cols) {
        if A[(row, k)] < 0 {
            A.add_multiple_coli_to_colj(-2, k, k);
        }
    }
    A = sort_cols_lex(&mut A,row);

    for _ in 0..5 {
        A.add_multiple_coli_to_colj(-A[(row, row)] / A[(row, row+1)], 1, 0);
        A = sort_cols_lex(&mut A,row);
    }

    return A;
}

/// Chooses the pivot for SNF calculations
///
/// The pivot is chosen from the submatrix of A in the given rows
/// and columns.
///
/// The pivot position and the weight of this pivot (smaller is better)
/// are returned.
fn find_pivot(A: &Matrix, rows: Range<usize>, cols: Range<usize>) -> Option<(usize, usize, i64)> {
    let mut indices = None;
    for i in rows {
        for j in cols.clone() {
            if A[(i, j)] == 0 {
                continue;
            }
            let weight = euclid_norm_squared(&A.col(j)) * euclid_norm_squared(&A.row(i));
            match indices {
                None => {
                    indices = Some((i, j, weight));
                }
                Some((_, _, min_weight)) => {
                    if weight < min_weight {
                        indices = Some((i, j, weight));
                    }
                }
            }
        }
    }
    indices
}

/// Calculates the smith normal form of A
///
/// Calculates the snf of A (done in place)
pub fn snf(A: &mut Matrix) -> Matrix {
    let mut A = A.clone();

    //this loops over the rows of A. t increments once A[(t,t)] is the only nonzero element of row t and column t of A.
    //probably only need to iterate up to A.rows-1 ?
    for t in 0..A.rows {
        let pivot = find_pivot(&A, t..A.rows, t..A.cols);
        match pivot {
            None => {
                continue;
            }
            Some((i, j, _)) => {
                A.swap_rows(i, t);
                A.swap_cols(j, t);
            }
        }

        //check that the pivot divides all the elements in column t (lower than row t)
        while !((&A.row(t)[t + 1..] == vec![0; A.cols - t - 1])
            && (&A.col(t)[t + 1..] == vec![0; A.rows - t - 1]))
        {
            for k in t + 1..A.rows {
                let a = A[(t, t)];
                let b = A[(k, t)];

                if b == 0 {
                    continue;
                } else if b % a == 0 {
                    let n = b / a;
                    A.add_multiple_rowi_to_rowj(-1 * n, t, k);
                } else {
                    let (gcd, sigma, tao) = extended_euclid(a, b);

                    let alpha = a / gcd;
                    let gamma = b / gcd;

                    let mut rowt = vec![0; A.cols];
                    let mut rowk = vec![0; A.cols];

                    for j in 0..A.cols {
                        rowt[j] = sigma * A[(t, j)] + tao * A[(k, j)];
                        rowk[j] = -gamma * A[(t, j)] + alpha * A[(k, j)];
                    }

                    for j in 0..A.cols {
                        A[(t, j)] = rowt[j];
                        A[(k, j)] = rowk[j];
                    }
                }
            }

            for k in t + 1..A.cols {
                let a = A[(t, t)];
                let b = A[(t, k)];

                if b == 0 {
                    continue;
                } else if b % a == 0 {
                    let n = b / a;
                    A.add_multiple_coli_to_colj(-1 * n, t, k);
                } else {
                    let (gcd, sigma, tao) = extended_euclid(a, b);

                    let alpha = a / gcd;
                    let gamma = b / gcd;
                    let mut L = Matrix::identity(A.cols);
                    L[(t, t)] = sigma;
                    L[(t, k)] = -1 * gamma;
                    L[(k, t)] = tao;
                    L[(k, k)] = alpha;

                    //the following lines are a horrible hack because i couldn't get matrix multiplication to work here

                    A = A.clone() * &L;
                }
            }
        }
    }
    for i in 0..A.cols {
        A[(i, i)] = A[(i, i)].abs();
    }
    return A;
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]

    fn test_snf() {
        let mut a: Matrix = vec![vec![-4, -6, 7], vec![2, 2, 4], vec![6, 6, 15]].into();
        let b = snf(&mut a);

        let snf_a: Matrix = vec![vec![2, 0, 0], vec![0, 1, 0], vec![0, 0, 6]].into();
        let snf_s: Matrix = vec![vec![-2, -1, 0], vec![1, 0, 0], vec![3, -3, 1]].into();
        let snf_t: Matrix = vec![vec![1, 1, 2], vec![0, 2, -15], vec![0, 1, -7]].into();

        assert_eq!(b, snf_a);

        a = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, 4, 16]].into();
        let b = snf(&mut a);

        let snf_a: Matrix = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 156]].into();
        let snf_s: Matrix = vec![vec![1, 0, 0], vec![-3, 9, -1], vec![5, -8, 1]].into();
        let snf_t: Matrix = vec![vec![1, 2, 2], vec![0, 1, 10], vec![0, 0, 1]].into();

        assert_eq!(b, snf_a);

        let mut a: Matrix = Matrix::identity(4);
        let b = snf(&mut a);

        let snf_a: Matrix = Matrix::identity(4);

        let snf_s: Matrix = Matrix::identity(4);
        let snf_t: Matrix = Matrix::identity(4);
        assert_eq!(b, snf_a);

        let mut a: Matrix = Matrix::new(4, 4);
        let b = snf(&mut a);

        let snf_a: Matrix = Matrix::new(4, 4);

        let snf_s: Matrix = Matrix::identity(4);
        let snf_t: Matrix = Matrix::identity(4);

        assert_eq!(b, snf_a);
    }
}
