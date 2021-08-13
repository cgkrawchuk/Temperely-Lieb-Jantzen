extern crate itertools;

use crate::matrix::*;

use temperley_lieb_cat::*;

use itertools::Itertools;

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

#[cfg(test)]
mod tests {

    use super::*;
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
