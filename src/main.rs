extern crate itertools;

use itertools::Itertools;
use std::io::{stdout, BufWriter};
use temperley_lieb_cat::*;

pub fn gramMatrix(n: usize, m: usize) -> Vec<Vec<i32>> {
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
    let mut gm: Vec<Vec<i32>> = vec![vec![0; x as usize]; x as usize];
    for i in 0..x {
        for j in 0..x {
            let (a, b) = &(monicDiagrams[i].involute()) * &monicDiagrams[j];
            if a == 0 {
                if b == TLDiagram::id(m) {
                    gm[i][j] = 1;
                } else {
                    gm[i][j] = 0;
                }
            } else {
                gm[i][j] = (2 as i32).pow(a as u32); //hopefully power of delta doesn't exceed u32!
            }
        }
    }
    return gm;
}

pub fn rowEchelonForm(matrix: &mut Vec<Vec<i32>>, p: i32) -> (Vec<Vec<i32>>, Vec<Vec<i32>>) {
    let mut matrix_out: Vec<Vec<i32>> = matrix.to_vec();
    let mut pivot = 0;
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();

    let mut identityMatrix: Vec<Vec<i32>> =
        vec![vec![0; column_count as usize]; row_count as usize];
    for k in 0..row_count {
        //hopefully row_count is less than column_count
        identityMatrix[k][k] = 1;
    }

    for r in 0..row_count {
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
        for j in 0..row_count {
            let temp1 = matrix_out[r][j];
            let temp2 = identityMatrix[r][j];
            matrix_out[r][j] = matrix_out[i][j];
            matrix_out[i][j] = temp1;
            identityMatrix[r][j] = identityMatrix[i][j];
            identityMatrix[i][j] = temp2;
        }
        let a = matrix_out[r][pivot];
        let modInv = mod_inv((((a % p) + p) % p) as isize, p as isize) as i32;
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

pub fn reduceModP(matrix: &mut Vec<Vec<i32>>, p: i32) -> Vec<Vec<i32>> {
    let mut matrix_out: Vec<Vec<i32>> = matrix.to_vec();
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();
    for i in 0..row_count {
        for j in 0..column_count {
            matrix_out[i][j] = ((matrix_out[i][j] % p) + p) % p;
        }
    }
    return matrix_out;
}

fn mod_inv(a: isize, module: isize) -> isize {
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

pub fn numZeroRows(matrix: &mut Vec<Vec<i32>>) -> i32 {
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

fn reform_inner_product(new_basis: &mut Vec<Vec<i32>>,old_gram: &mut Vec<Vec<i32>>) -> Vec<Vec<i32>> {
    let mut matrix_out: Vec<Vec<i32>> =
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
    let m = 5;
    let n = 3;
    let p = 5;
    let mut G = gramMatrix(m, n);
    println!("Gram Matrix of {0}, {1} is :", m, n);
    println!("{:?}", &G);

    let mut G = reduceModP(&mut G, p);
    println!("Reduced mod {}", p);
    println!("{:?}", &G);

    let (mut reducedMatrix,mut transformMatrix) = rowEchelonForm(&mut G, p);
    
    println!("Row echelon form: \n{:?}", reduceModP(&mut reducedMatrix,p));
    println!("transformation matrix:\n{:?}", reduceModP(&mut transformMatrix,p));

    let rank = (reducedMatrix.len() as i32) - numZeroRows(&mut reducedMatrix.clone());
    println!("dimension of head is: {}", &rank);

    let basisOfRad = &transformMatrix[(rank as usize)..];
    let mut bor: Vec<Vec<i32>> = Vec::from(basisOfRad);
    println!("{:?}", &basisOfRad);

    let newG = reform_inner_product(&mut bor, &mut G);
    println!("{:?}", &newG);
}
