extern crate itertools;

use itertools::Itertools;
use temperley_lieb_cat::*;
use tl_jantzen::{binom, extended_euclid, Matrix};

use std::io::{self, BufRead};

use std::cmp::min;
use std::env;

pub fn sort_cols_max(A: &mut Matrix) -> (Matrix) {
    let mut A = A.clone();
    for i in 0..A.cols {
        for j in 0..(A.cols - i - 1) {
            if A.return_col_j(j).iter().max().unwrap() > A.return_col_j(j + 1).iter().max().unwrap()
            {
                A.swap_cols(j, j + 1);
            }
        }
    }
    return A;
}
pub fn euclid_norm_squared(v: Vec<i64>) -> i64 {
    let mut ans = 0;
    for x in v.iter() {
        ans += x * x;
    }
    ans
}

pub fn max_norm(A: &Matrix) -> i64 {
    let data = A.entries();
    return *data.iter().max().unwrap();
}

pub fn sign(x: i64) -> (i64) {
    if x < 0 {
        return -1;
    } else if x > 0 {
        return 1;
    } else {
        return 0;
    }
}

pub fn normcol(B: &mut Matrix) -> Matrix {
    let mut A = B.clone();
    sort_cols_max(&mut A);

    let mut norm_s = max_norm(&A);
    let mut norm_a = norm_s;
    while norm_s <= norm_a {
        norm_a = norm_s;
        for k in (1..A.cols).rev() {
            if sign(A[(0, k)]) == sign(A[(0, k - 1)]) {
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

pub fn sort_cols_lex(B: &mut Matrix) -> Matrix {
    let mut A = B.clone();
    for i in 0..A.cols {
        for j in 0..(A.cols - i - 1) {
            if A.return_col_j(j) < A.return_col_j(j + 1) {
                A.swap_cols(j, j + 1);
            }
        }
    }

    return A;
}

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

pub fn snf(B: &mut Matrix) -> (Matrix) {
    let mut A = B.clone();
    let t = 0;

    for t in (0..A.rows - 1) {
        let mut count = 0;
        let mut minValue = i64::MAX;
        let mut indicies = (t, t);
        for i in t..A.rows {
            for j in t..A.cols {
                if A[(i, j)] == 0 {
                    count += 1;
                    continue;
                } else if euclid_norm_squared(A.return_col_j(j))
                    * euclid_norm_squared(A.return_row_i(i))
                    < minValue
                {
                    minValue = euclid_norm_squared(A.return_col_j(j))
                        * euclid_norm_squared(A.return_row_i(i));
                    indicies.0 = i;
                    indicies.1 = j;
                }
            }
        }

        if count == (A.rows - t) * (A.cols - t) {
            break;
        }

        A.swap_rows(t, indicies.0);
        A.swap_cols(t, indicies.1);

        let zero_row = vec![0; A.cols - t - 1];

        let zero_col = vec![0; A.rows - t - 1];

        while &A.return_row_i(t)[t + 1..] != zero_row || &A.return_col_j(t)[t + 1..] != zero_col {
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
                    let mut L = Matrix::identity(A.cols);
                    L[(t, t)] = sigma;
                    L[(t, k)] = -1 * gamma;
                    L[(k, t)] = tao;
                    L[(k, k)] = alpha;

                    //the following lines are a horrible hack because i couldn't get matrix multiplication to work here

                    A = A.clone() * &L;
                }
            }

            indicies = (t, t);
            minValue = i64::MAX;
            let i = t;
            for j in (t..A.cols) {
                if A[(i, j)] == 0 {
                    count += 1;
                    continue;
                } else if euclid_norm_squared(A.return_col_j(j))
                    * euclid_norm_squared(A.return_row_i(i))
                    < minValue
                {
                    minValue = euclid_norm_squared(A.return_col_j(j))
                        * euclid_norm_squared(A.return_row_i(i));
                    indicies.0 = i;
                    indicies.1 = j;
                }
            }

            let j = t;
            for i in (t..A.rows) {
                if A[(i, j)] == 0 {
                    count += 1;
                    continue;
                } else if euclid_norm_squared(A.return_col_j(j))
                    * euclid_norm_squared(A.return_row_i(i))
                    < minValue
                {
                    minValue = euclid_norm_squared(A.return_col_j(j))
                        * euclid_norm_squared(A.return_row_i(i));
                    indicies.0 = i;
                    indicies.1 = j;
                }
            }

            A.swap_rows(t, indicies.0);
            A.swap_cols(t, indicies.1);
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
