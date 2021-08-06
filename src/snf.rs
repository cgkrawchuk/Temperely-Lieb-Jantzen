use crate::matrix::*;
use crate::util::*;

///Calculates the smith normal form of A
///
///Calculates the snf of A (done in place) and returns the invertible matricies
/// S and T such that S*(snf(A))*T=A
pub fn snf(A: &mut Matrix) -> (Matrix, Matrix, Matrix) {
    let mut A = A.clone();
    let mut S = Matrix::identity(A.rows);
    let mut T = Matrix::identity(A.cols);

    //this loops over the rows of A. t increments once A[(t,t)] is the only nonzero element of row t and column t of A.
    //probably only need to iterate up to A.rows-1 ?
    for t in 0..A.rows {
        //finding a pivot:
        //first search down column t starting at row t
        for i in t..A.rows {
            if A[(i, t)] != 0 && i == t {
                break;
            } else if A[(i, t)] != 0 && i > t {
                A.swap_rows(i, t);
            }
        }

        //if we haven't found a pivot, next we check row t starting at column t
        if A[(t, t)] == 0 {
            for i in t..A.cols {
                if A[(t, i)] != 0 {
                    A.swap_cols(t, i);
                }
            }
        }
        //this must mean row t and column t are all 0, so we go to the next row and try again
        if A[(t, t)] == 0 {
            continue;
        }

        //check that the pivot divides all the elements in column t (lower than row t)
        loop {
            let mut made_changes = false;

            for k in t + 1..A.rows {
                let a = A[(t, t)];
                let b = A[(k, t)];

                if b == 0 {
                    continue;
                } else if b % a == 0 {
                    let n = b / a;
                    A.add_multiple_rowi_to_rowj(-1 * n, t, k);
                    S.add_multiple_coli_to_colj(n, k, t);
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

                    let mut colt = vec![0; S.rows];
                    let mut colk = vec![0; S.rows];

                    for j in 0..S.rows {
                        colt[j] = alpha * S[(j, t)] + gamma * S[(j, k)];
                        colk[j] = -tao * S[(j, t)] + sigma * S[(j, k)];
                    }

                    for j in 0..S.rows {
                        S[(j, t)] = colt[j];
                        S[(j, k)] = colk[j];
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
                    T.add_multiple_rowi_to_rowj(n, k, t);
                    made_changes = true;
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

                    let mut rowt = vec![0; T.cols];
                    let mut rowk = vec![0; T.cols];

                    for j in 0..T.cols {
                        rowt[j] = alpha * T[(t, j)] + gamma * T[(k, j)];
                        rowk[j] = -tao * T[(t, j)] + sigma * T[(k, j)];
                    }

                    for j in 0..T.cols {
                        T[(t, j)] = rowt[j];
                        T[(k, j)] = rowk[j];
                    }

                    made_changes = true;
                }
            }

            if !made_changes {
                break;
            }
        }
    }
    return (S, A, T);
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]

    //put in zero matrix, put in identity
    fn test_snf() {
        let mut a: Matrix = vec![vec![-4, -6, 7], vec![2, 2, 4], vec![6, 6, 15]].into();
        let (s, b, t) = snf(&mut a);

        let snf_a: Matrix = vec![vec![2, 0, 0], vec![0, 1, 0], vec![0, 0, 6]].into();
        let snf_s: Matrix = vec![vec![-2, -1, 0], vec![1, 0, 0], vec![3, -3, 1]].into();
        let snf_t: Matrix = vec![vec![1, 1, 2], vec![0, 2, -15], vec![0, 1, -7]].into();

        assert_eq!(b, snf_a);
        assert_eq!(s, snf_s);
        assert_eq!(t, snf_t);

        a = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, 4, 16]].into();
        let (s, b, t) = snf(&mut a);

        let snf_a: Matrix = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 156]].into();
        let snf_s: Matrix = vec![vec![1, 0, 0], vec![-3, 9, -1], vec![5, -8, 1]].into();
        let snf_t: Matrix = vec![vec![1, 2, 2], vec![0, 1, 10], vec![0, 0, 1]].into();

        assert_eq!(b, snf_a);
        assert_eq!(s, snf_s);
        assert_eq!(t, snf_t);

        let mut a: Matrix = Matrix::identity(4);
        let (s, b, t) = snf(&mut a);

        let snf_a: Matrix = Matrix::identity(4);
        
        let snf_s: Matrix=Matrix::identity(4);
        let snf_t: Matrix=Matrix::identity(4);
        assert_eq!(b, snf_a);
        assert_eq!(s, snf_s);
        assert_eq!(t, snf_t);

        let mut a: Matrix = Matrix::new(4, 4);
        let (s, b, t) = snf(&mut a);

        let snf_a: Matrix = Matrix::new(4,4);
        
        let snf_s:Matrix= Matrix::identity(4);
        let snf_t: Matrix=Matrix::identity(4);

        assert_eq!(b, snf_a);
        assert_eq!(s, snf_s);
        assert_eq!(t, snf_t);
    }
}
