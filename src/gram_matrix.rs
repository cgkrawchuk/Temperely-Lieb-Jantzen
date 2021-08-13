extern crate itertools;

use crate::matrix::*;

use crate::util::*;

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