extern crate itertools;

use crate::matrix::*;

use crate::util::*;

///Returns a list of the p-parts of elementary divisors of A
///
///Returns the list 'dimensions' where dimensions[i] is the number of elementary
///divisiors divisible by p^i subtract the number of elementary divisible by p^(i-1)
pub fn elem_div(A: &Matrix, p: i64) -> Vec<i64> {
    //let r = rank(&A,71) as i64;
    let r = A.cols as i64;
    //println!("rank is: {:?}", r);

    let n = A.cols;
    let mut Tp = Vec::new();
    for i in 0..n {
        Tp.push(i);
    }
    let mut A2 = A.clone();
    let mut A1: Matrix = Matrix::new(0, 0);
    let mut inv = Vec::new();
    let mut res = Vec::new();

    while A1.rows < r as usize {
        let i0 = A2.rows;
        let mut B1: Matrix = Matrix::new(0, 0);
        let mut A21 = 0;

        for ii in 0..i0 {
            let mut vv = A2.row(ii);

            for i in 0..A1.rows {
                let c = (vv[Tp[i]] * inv[i]) % p;
                if c != 0 {
                    for y in 0..vv.len() {
                        vv[y] = vv[y] - c * A1.row(i)[y];
                    }
                }
            }
            let mut pos = A1.rows;

            while (pos < n) && (vv[Tp[pos]] % p == 0) {
                pos = pos + 1;
            }

            if pos >= n {
                A21 = A21 + 1;
                for y in 0..vv.len() {
                    vv[y] = vv[y] / p;
                }
                B1 = B1.add_row(vv);
            } else {
                let i1 = A1.rows;
                A1 = A1.add_row(vv);

                if pos != i1 {
                    let x = Tp[pos];
                    Tp[pos] = Tp[i1];
                    Tp[i1] = x;
                }

                inv.push(mod_inv(A1.row(i1)[Tp[i1]], p));
            }
        }

        A2 = B1;

        res.push(r - A1.rows as i64);
    }
    let mut ans = vec![r];
    let mut dimensions = Vec::new();
    ans.append(&mut res);
    for i in 1..ans.len() {
        dimensions.push(ans[i - 1] - ans[i]);
    }
    return dimensions;
}
