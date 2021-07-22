extern crate itertools;

mod modInverse;

use itertools::Itertools;
use std::io::{self, BufRead};
use temperley_lieb_cat::*;
use crate::modInverse::mod_inv;

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

pub fn row_echelon_form(matrix: &Vec<Vec<i64>>, p: i64) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let mut matrix_out: Vec<Vec<i64>> = matrix.to_vec();
    let mut pivot = 0;
    let row_count = matrix_out.len();
    let column_count = matrix_out[0].len();

    let mut identity_matrix: Vec<Vec<i64>> =
        vec![vec![0; column_count as usize]; row_count as usize];
    for k in 0..row_count {
        
        identity_matrix[k][k] = 1;
    }
   
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
        for j in r+1..row_count {
           
                let hold = matrix_out[j][pivot];
                for k in 0..column_count {
                    matrix_out[j][k] = matrix_out[j][k] - (hold * matrix_out[r][k]);
                    identity_matrix[j][k] = identity_matrix[j][k] - (hold * identity_matrix[r][k]);
                
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



pub fn num_zero_rows(matrix: &Vec<Vec<i64>>) -> i64 {
    let mut num_zeroes = 0;
    let zero_vec = vec![0; matrix[0].len()];
    for x in matrix.iter() {
        if x == &zero_vec {
            num_zeroes += 1;
        }
    }
    return num_zeroes;
}

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

fn recursive_ops(m: usize, n: usize, p: i64) {
    let mut g = gram_matrix(m, n);
    
    let mut h = reduce_mod_p(&mut g, p);
    
    let (mut reduced_matrix, mut transform_matrix) = row_echelon_form(&mut h, p);

    loop{
        let  nzr = num_zero_rows(&reduced_matrix);
        let  rank = (reduced_matrix.len() as i64) - nzr;
        println!("dimension of head is: {}", &rank);
        if nzr ==0 { break;}
        let basis_of_rad = &transform_matrix[(rank as usize)..]; 
        g = reform_inner_product(basis_of_rad, &mut g);
        for i in 0..g.len() {
            for j in 0..g.len() {
                g[i][j] = g[i][j] / p;
            }
        }
        h = reduce_mod_p(&mut g, p);
        let (temp1, temp2) = row_echelon_form(&mut h, p);
        reduced_matrix = temp1;
        transform_matrix = temp2;
        
        
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
#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_mod_inv() {
        assert_eq!(mod_inv(3, 5), 2);
    }

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
    fn test_row_echelon_form(){
        let m = vec![
            vec![2, 0, -13, 0],
            vec![4, 1, 0, 0],
            vec![0, 8, 1, 0],
            vec![5, -1, 0, 1],
        ];
        let ans = (vec![vec![1, 0, 1, 0], vec![0, 1, 2, 0], vec![0, 0, 0, 1], vec![0, 0, 0, 0]], vec![vec![2, 0, 0, 0], vec![-2, 1, 0, 0], vec![0, -2, 0, 1], vec![4, -2, 1, 0]]);
    assert_eq!(row_echelon_form(&m,3),ans);

    }

    #[test]
    fn test_gram_matrix(){
        let m: Vec<Vec<i64>>= vec![vec![4, 2, 0, 0, 2, 1, 0, 0, 2, 0, 0, 0, 0, 0], vec![2, 4, 2, 0, 1, 2, 1, 0, 1, 0, 0, 0, 0, 0], vec![0, 2, 4, 2, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0], vec![0, 0, 2, 4, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0], vec![2, 1, 0, 0, 4, 2, 0, 0, 1, 0, 0, 1, 0, 0], vec![1, 2, 1, 0, 2, 4, 2, 0, 2, 1, 0, 2, 0, 0], vec![0, 1, 2, 1, 0, 2, 4, 2, 1, 2, 1, 1, 0, 0], vec![0, 0, 1, 2, 0, 0, 2, 4, 0, 1, 2, 0, 0, 0], vec![2, 1, 0, 0, 1, 2, 1, 0, 4, 2, 0, 1, 0, 1], vec![0, 0, 0, 0, 0, 1, 2, 1, 2, 4, 2, 2, 1, 2], vec![0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 4, 1, 2, 1], vec![0, 0, 0, 0, 1, 2, 1, 0, 1, 2, 1, 4, 2, 1], vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 4, 2], vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 4]];
        assert_eq!(gram_matrix(7,3), m);
    }

}
