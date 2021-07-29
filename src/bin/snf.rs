use tl_jantzen::{Matrix,extended_euclid};



///Calculates the smith normal form of A 
///
///Calculates the snf of A (done in place) and returns the invertible matricies
/// S and T such that S*(snf(A))*T=A
pub fn snf(A: &mut Matrix) -> (Matrix, Matrix) {
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
                    let mut L = Matrix::identity(A.rows);
                    L[(t, t)] = sigma;
                    L[(t, k)] = tao;
                    L[(k, t)] = -1 * gamma;
                    L[(k, k)] = alpha;
                    let mut L_inverse = Matrix::identity(A.rows);
                    L_inverse[(t, t)] = alpha;
                    L_inverse[(t, k)] = -1 * tao;
                    L_inverse[(k, t)] = gamma;
                    L_inverse[(k, k)] = sigma;
                    *A = &L * A;

                    S = S * L_inverse;
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
                    L[(k, t)] = tao;
                    L[(t, k)] = -1 * gamma;
                    L[(k, k)] = alpha;

                    let mut L_inverse = Matrix::identity(A.cols);
                    L_inverse[(t, t)] = alpha;
                    L_inverse[(k, t)] = -1 * tao;
                    L_inverse[(t, k)] = gamma;
                    L_inverse[(k, k)] = sigma;
                    
                    //the following lines are a horrible hack because i couldn't get matrix multiplication to work here
                    let mut m: Matrix = Matrix::new(A.cols,A.cols);
                    
                    for i in 0..A.rows {
                        for j in 0..A.cols {
                            for k in 0..A.cols {
                                m[(i, j)] += A[(i, k)] * L[(k, j)];
                            }
                        }
                    }
                    *A = m;
                    

                    T = L_inverse * T;

                    made_changes = true;
                }
            }

            if !made_changes {
                break;
            }
        }
    }
    return (S, T);
}

fn main() {
    let mut a: Matrix = vec![vec![-4, -6, 7], vec![2, 2, 4], vec![6, 6, 15]].into();
    let (s, t) = snf(&mut a);
    println!("a is: {}", a);
    println!("s is: {}", s);
    println!("t is: {}", t);


}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_snf(){
        let mut a: Matrix = vec![vec![-4, -6, 7], vec![2, 2, 4], vec![6, 6, 15]].into();
        let (s, t) = snf(&mut a);

        let snf_a: Matrix = vec![vec![2, 0, 0], vec![0, 1, 0], vec![0, 0, 6]].into();
        let snf_s: Matrix = vec![vec![-2, -1, 0], vec![1, 0, 0], vec![3, -3, 1]].into();
        let snf_t: Matrix = vec![vec![1, 1, 2], vec![0, 2, -15], vec![0, 1, -7]].into();

        assert_eq!(a, snf_a);
        assert_eq!(s, snf_s);
        assert_eq!(t, snf_t);

        a = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, 4, 16]].into();
        let (s, t) = snf(&mut a);

        let snf_a: Matrix = vec![vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 156]].into();
        let snf_s: Matrix = vec![vec![1, 0, 0], vec![-3, 9, -1], vec![5, -8, 1]].into();
        let snf_t: Matrix = vec![vec![1, 2, 2], vec![0, 1, 10], vec![0, 0, 1]].into();

        assert_eq!(a, snf_a);
        assert_eq!(s, snf_s);
        assert_eq!(t, snf_t);
    }

    
}
