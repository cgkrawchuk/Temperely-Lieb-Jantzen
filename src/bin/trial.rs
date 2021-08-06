extern crate itertools;

use itertools::Itertools;
use std::ops::Range;
use temperley_lieb_cat::*;
use tl_jantzen::{binom, extended_euclid, Matrix};

use std::io::{self, BufRead};

use std::cmp::min;
use std::env;

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
pub fn sort_cols_lex(B: &mut Matrix) -> Matrix {
    let mut A = B.clone();
    for i in 0..A.cols {
        for j in 0..(A.cols - i - 1) {
            if A.col(j) < A.col(j + 1) {
                A.swap_cols(j, j + 1);
            }
        }
    }

    return A;
}

///Reduces the entries of a matrix using the rosser method
pub fn rosser(B: &mut Matrix) -> Matrix {
    let mut A = B.clone();
    for k in (0..A.cols) {
        if A[(0, k)] < 0 {
            A.add_multiple_coli_to_colj(-2, k, k);
        }
    }
    A = sort_cols_lex(&mut A);

    while A[(0, 1)] != 0 {
        A.add_multiple_coli_to_colj(-A[(0, 0)] / A[(0, 1)], 1, 0);
        A = sort_cols_lex(&mut A);
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

/// Calculates the snf as before but with new method for choosing pivot
pub fn snf(B: &mut Matrix) -> (Matrix) {
    let mut A = B.clone();
    let t = 0;

    for t in (0..A.rows - 1) {
        let indices = find_pivot(&A, t..A.rows, t..A.cols);
        match indices {
            None => {
                break;
            }
            Some((i, j, _)) => {
                A.swap_rows(t, i);
                A.swap_cols(t, j);
            }
        }

        let zero_row = vec![0; A.cols - t - 1];

        let zero_col = vec![0; A.rows - t - 1];

        while &A.row(t)[t + 1..] != zero_row || &A.col(t)[t + 1..] != zero_col {
            for k in t + 1..A.rows {
                let a = A[(t, t)];
                let b = A[(k, t)];

                if b == 0 {
                    continue;
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
                } else {
                    let (gcd, sigma, tao) = extended_euclid(a, b);

                    let alpha = a / gcd;
                    let gamma = b / gcd;

                    let mut colt = vec![0; A.rows];
                    let mut colk = vec![0; A.rows];

                    for i in 0..A.rows {
                        colt[i] = sigma * A[(i, t)] + tao * A[(i, k)];
                        colk[i] = -gamma * A[(i, t)] + alpha * A[(i, k)];
                    }

                    for i in 0..A.rows {
                        A[(i, t)] = colt[i];
                        A[(i, k)] = colk[i];
                    }
                }
            }

            let indices_row = find_pivot(&A, t..t + 1, t..A.cols).unwrap();
            let indices_col = find_pivot(&A, t..A.rows, t..t + 1).unwrap();

            let indices;
            if indices_row.2 < indices_col.2 {
                indices = (indices_row.0, indices_row.1);
            } else {
                indices = (indices_col.0, indices_col.1);
            }

            A.swap_rows(t, indices.0);
            A.swap_cols(t, indices.1);
        }
    }

    return A;
}

fn main() {
    let mut a: Matrix = vec![vec![-4, 2, 7], vec![2, 7, 4], vec![6, 6, 3]].into();

    println!("{}", snf(&mut a));

    a = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, 4, 16]].into();

    println!("{}", snf(&mut a));

    a = Matrix::new(3, 3);

    println!("{}", snf(&mut a));

    let mut g: Matrix = vec![
        vec![4, 2, 0, 2, 1, 0, 2, 0, 0],
        vec![2, 4, 2, 1, 2, 1, 1, 0, 0],
        vec![0, 2, 4, 0, 1, 2, 0, 0, 0],
        vec![2, 1, 0, 4, 2, 0, 1, 0, 1],
        vec![1, 2, 1, 2, 4, 2, 2, 1, 2],
        vec![0, 1, 2, 0, 2, 4, 1, 2, 1],
        vec![2, 1, 0, 1, 2, 1, 4, 2, 1],
        vec![0, 0, 0, 0, 1, 2, 2, 4, 2],
        vec![0, 0, 0, 1, 2, 1, 1, 2, 4],
    ]
    .into();

    println!("{}", snf(&mut g));
}
