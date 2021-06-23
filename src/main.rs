extern crate itertools;

use itertools::Itertools;
use std::io::{stdout, BufWriter};
use temperley_lieb_cat::*;

pub fn gramMatrix(n: usize, m: usize) -> Vec<Vec<i64>> {
    fn ok(tab: &Vec<usize>) -> bool {
        for (c, i) in tab.iter().enumerate() {
            if *i < 2 * (c + 1) {
                return false;
            }
        }
        true
    }
    let mut monicDiagrams = Vec::new();
    let prop = m;
    let rn = (n - prop) / 2;
    for lt in (1..(n + 1)).combinations(rn).filter(ok) {
        monicDiagrams.push(TLDiagram::new(
            n,
            m,
            lt.iter().fold(0, |l, n| l | (1 << n)),
            0,
        ));
    }

    let x = monicDiagrams.len();
    let mut gm: Vec<Vec<i64>> = vec![vec![0; x as usize]; x as usize];
    for i in 0..x {
        for j in 0..x {
            let (a, b) = &(monicDiagrams[i].involute()) * &monicDiagrams[j];
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

pub fn rowEchelonForm(matrix: &mut Vec<Vec<i64>>, p: i64) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let mut matrix_out: Vec<Vec<i64>> = matrix.to_vec();
    let mut pivot = 0;
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();

    let mut identityMatrix: Vec<Vec<i64>> =
        vec![vec![0; column_count as usize]; row_count as usize];
    for k in 0..row_count {
        //hopefully row_count is less than column_count
        identityMatrix[k][k] = 1;
    }

    for r in 0..row_count {
        matrix_out = reduceModP(&mut matrix_out, p);
        println!("{:?}", matrix_out);
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
        matrix_out = swapRows(&mut matrix_out, r, i, row_count);
        identityMatrix = swapRows(&mut identityMatrix, r, i, row_count);

        let a = matrix_out[r][pivot];
        let modInv = mod_inv((((a % p) + p) % p), p);
        for j in 0..column_count {
            matrix_out[r][j] = matrix_out[r][j] * modInv;
            identityMatrix[r][j] = identityMatrix[r][j] * modInv;
        }
        for j in 0..row_count {
            if j != r {
                let hold = matrix_out[j][pivot];
                for k in 0..column_count {
                    matrix_out[j][k] = matrix_out[j][k] - (hold * matrix_out[r][k]);
                    identityMatrix[j][k] = identityMatrix[j][k] - (hold * identityMatrix[r][k]);
                }
            }
        }
        pivot = pivot + 1;
    }
    (matrix_out, identityMatrix)
}

pub fn swapRows(
    matrixIn: &mut Vec<Vec<i64>>,
    r: usize,
    i: usize,
    row_count: usize,
) -> Vec<Vec<i64>> {
    let mut matrix: Vec<Vec<i64>> = matrixIn.to_vec();
    for j in 0..row_count {
        let temp = matrix[r][j];
        matrix[r][j] = matrix[i][j];
        matrix[i][j] = temp;
    }
    matrix
}

pub fn reduceModP(matrix: &Vec<Vec<i64>>, p: i64) -> Vec<Vec<i64>> {
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

pub fn numZeroRows(matrix: &Vec<Vec<i64>>) -> i64 {
    let mut numZeroes = 0;
    let row_count = matrix.len();
    let column_count = matrix[0].len();
    let zero_vec = vec![0; column_count];
    for x in 0..row_count {
        if matrix[row_count - x - 1] == zero_vec {
            numZeroes = numZeroes + 1;
        } else {
            break;
        } //remove for more general matrix
    }
    return numZeroes;
}

fn reform_inner_product(
    new_basis: &mut Vec<Vec<i64>>,
    old_gram: &mut Vec<Vec<i64>>,
) -> Vec<Vec<i64>> {
    let mut matrix_out: Vec<Vec<i64>> =
        vec![vec![0; new_basis[0].len() as usize]; new_basis.len() as usize];
    for i in 0..new_basis.len() {
        for j in 0..new_basis.len() {
            let basis_i = new_basis[i].clone();
            let basis_j = new_basis[j].clone();
            let mut inner = 0;
            for k in 0..old_gram.len() {
                for l in 0..old_gram[0].len() {
                    inner += old_gram[k][l] * basis_i[k] * basis_j[l];
                }
            }
            matrix_out[i][j] = inner;
        }
    }

    matrix_out
}

fn main() {
    let m = 6;
    let n = 4;
    let p = 3;
    let mut G = gramMatrix(m, n);
    println!("Gram Matrix of {0}, {1} is :", m, n);
    println!("{:?}", &G);

    let mut G = reduceModP(&mut G, p);
    println!("Reduced mod {}", p);
    println!("{:?}", &G);

    let (mut reducedMatrix, mut transformMatrix) = rowEchelonForm(&mut G, p);

    println!("Row echelon form: \n{:?}", reduceModP(&reducedMatrix, p));
    println!(
        "transformation matrix:\n{:?}",
        reduceModP(&transformMatrix, p)
    );

    let rank = (reducedMatrix.len() as i64) - numZeroRows(&reducedMatrix.clone());
    println!("dimension of head is: {}", &rank);

    let basisOfRad = &transformMatrix[(rank as usize)..];
    let mut bor: Vec<Vec<i64>> = Vec::from(basisOfRad);
    println!("{:?}", &basisOfRad);

    let newG = reform_inner_product(&mut bor, &mut G);
    println!("{:?}", &newG);
}
