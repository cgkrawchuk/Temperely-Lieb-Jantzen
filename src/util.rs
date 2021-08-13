use crate::matrix::*;

/// Calculate the binomial coefficients
///
/// Calculates $n$ choose $k$ where this is understood
/// to be zero if $k<0$, $n<0$ or $n>k$.
pub fn binom(n: i64, k: i64) -> i64 {
    if k < 0 || n < 0 {
        return 0;
    }
    let mut res = 1;
    for i in 0..k {
        res = (res * (n - i)) / (i + 1);
    }
    res
}
/// Calculates the mod-inverse of 'a'
///
/// Calculates the multiplicative inverse of 'a' modulo 'module'
/// Asserts that 'a' is not divisible by 'module'
pub fn mod_inv(c: i64, module: i64) -> i64 {
    assert!(c % module != 0, "number is 0 mod...");
    let a = ((c % module) + module) % module;
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

///Function for extended euclidian algorithm
///
/// Returns the gcd of x and y and the corresponding Bezout coefficients
pub fn extended_euclid(x: i64, y: i64) -> (i64, i64, i64) {
    let mut s = 0;
    let mut old_s = 1;
    let mut t = 1;
    let mut old_t = 0;
    let mut r = y;
    let mut old_r = x;
    let mut temp;
    while r != 0 {
        let quotient = old_r / r;
        temp = r;
        r = old_r - quotient * r;
        old_r = temp;

        temp = s;
        s = old_s - quotient * s;
        old_s = temp;

        temp = t;
        t = old_t - quotient * t;
        old_t = temp;
    }
    if old_r < 0 {
        old_r = old_r * (-1);
        old_s = old_s * (-1);
        old_t = old_t * (-1);
    }
    return (old_r, old_s, old_t);
}

pub fn rank(matrix: &Matrix, p: i64) -> usize {
    let id = Matrix::identity(matrix.rows);
    row_echelon_form(matrix, &id, p).2
}

/// Reduces a matrix modulo p
///
/// Reduces each entry in a matrix mod p in-place.
pub fn reduce_mod_p(matrix: &mut Matrix, p: i64) {
    for i in 0..matrix.rows {
        for j in 0..matrix.cols {
            matrix[(i, j)] = ((matrix[(i, j)] % p) + p) % p;
        }
    }
}

/// Calculate the row echelon form of a matrix
///
/// Calculates the unreduced row echelon form of a matrix
/// while performing the same operations on the identity matrix
/// with the same dimensions. Returns both matricies. Division is
/// performed modulo the argument p
pub fn row_echelon_form(mx: &Matrix, id: &Matrix, p: i64) -> (Matrix, Matrix, usize) {
    let mut identity_matrix = id.clone();
    let mut C: Matrix = mx.clone();
    let mut pivot_row = 0;

    'col_loop: for column in 0..mx.cols {
        let mut i = pivot_row;
        while C[(i, column)] % p == 0 {
            i += 1;
            if i == mx.rows {
                continue 'col_loop;
            }
        }

        C.swap_rows(pivot_row, i);

        identity_matrix.swap_rows(pivot_row, i);

        let q = C[(pivot_row, column)];

        let mod_inverse = mod_inv(q, p);

        for j in pivot_row + 1..mx.rows {
            let hold = ((C[(j, column)] % p) + p) % p;
            for k in 0..mx.cols {
                C[(j, k)] -= hold * C[(pivot_row, k)] * mod_inverse;
            }
            for k in 0..identity_matrix.cols {
                identity_matrix[(j, k)] -= hold * identity_matrix[(pivot_row, k)] * mod_inverse;
            }
        }
        pivot_row += 1;
        if pivot_row == mx.rows {
            break;
        }
    }

    (C, identity_matrix, pivot_row)
}

/// Converts an i64 into the vector of its digits in base p
///
/// Returns the list of values of 'r' in base 'p'
/// mutiplied by the corresponding power of p. The values are
/// given in descending order  
pub fn convert_base_p(mut r: i64, p: i64) -> Vec<i64> {
    assert!(2 <= p);
    let mut digits: Vec<i64> = Vec::new();
    let mut index = 0;
    while r > 0 {
        digits.push((r % p) * (p.pow(index as u32)));
        r /= p;
        index += 1;
    }
    digits.reverse();
    digits
}

/// Calculates the 'p'-adic valuation of 'a'
///
/// Returns the highest power of 'p' that divides
/// 'a'. Asserts that 'a' is not 0.
pub fn p_adic_val(mut a: i64, p: i64) -> i64 {
    assert!(a != 0, "vp(0)=infinity");
    let mut index = 0;
    while a % p == 0 {
        a /= p;
        index += 1;
    }
    index
}

/// Checks for a solution to the knapsack problem
///
/// Accepts a vector of pairs and checks if the sum of
/// the second entries of a subset of these pairs totals
/// to 'sum'. Returns TRUE if a solutions exists and
/// FALSE otherwise
pub fn knapsack_sols_exists(v: &[(i64, i64)], sum: i64) -> bool {
    if v.is_empty() {
        sum == 0
    } else {
        let with_first = knapsack_sols_exists(&v[1..], sum - v[0].1);
        let without_first = knapsack_sols_exists(&v[1..], sum);
        if sum != 0 {
            assert!(!(with_first & without_first), "multiple solutions");
        }
        with_first | without_first
    }
}

/// Returns a solution to the knapsack problem if it exists
///
/// Accepts a vector of pairs and returns a subset of these pairs if
/// the sum of their second entry totals to 'sum'
pub fn knapsack_sols(v: &[(i64, i64)], mut sum: i64) -> Vec<(i64, i64)> {
    let mut ans: Vec<(i64, i64)> = Vec::new();
    knapsack_sols_exists(&v, sum);
    for (n, x) in v.iter().enumerate() {
        if knapsack_sols_exists(&v[(n + 1)..], sum - x.1) {
            ans.push(*x);
            sum -= x.1;
        }
    }
    ans
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_mod_inv() {
        assert_eq!(mod_inv(3, 5), 2);
        assert_eq!(mod_inv(-58, 3), 2);
    }

    #[test]
    fn test_rank() {
        let a: Matrix = vec![vec![1, 2, 3], vec![-1, -1, 0], vec![3, 8, 9]].into();
        assert_eq!(rank(&a, 71), 3);
        let a: Matrix = Matrix::identity(3);
        assert_eq!(rank(&a, 71), 3);

        let a = Matrix::new(3, 3);
        assert_eq!(rank(&a, 71), 0);

        let a: Matrix = vec![vec![1, 2, 1], vec![-2, -3, 1], vec![3, 5, 0]].into();
        assert_eq!(rank(&a, 71), 2);
    }

    #[test]
    fn test_reduce_mod_p() {
        let mut a: Matrix = vec![vec![1, 2, 1], vec![-2, -3, 1], vec![3, 5, 0]].into();
        let ans: Matrix = vec![vec![1, 2, 1], vec![1, 0, 1], vec![0, 2, 0]].into();
        reduce_mod_p(&mut a, 3);
        assert_eq!(a, ans);

        let mut m: Matrix = vec![
            vec![2, 0, -13, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, -1, 0, 1],
        ]
        .into();

        let n: Matrix = vec![
            vec![2, 0, 2, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 2, 0, 1],
        ]
        .into();
        reduce_mod_p(&mut m, 3);
        assert_eq!(m, n);
    }

    #[test]
    fn test_convert_base_p() {
        let v = vec![4, 2, 1];
        assert_eq!(convert_base_p(7, 2), v);
        assert_eq!(convert_base_p(0, 2), vec![]);
        let v = vec![243, 0, 54, 0, 0, 2];
        assert_eq!(convert_base_p(299, 3), v);
    }

    #[test]
    #[should_panic]
    fn test_p_adic_val_infinity() {
        p_adic_val(0, 2);
    }
    #[test]
    fn test_p_adic_val() {
        assert_eq!(p_adic_val(54, 3), 3);
        assert_eq!(p_adic_val(1, 2), 0);
    }

    #[test]
    fn test_row_echelon_form() {
        let m1: Matrix = vec![
            vec![2, 0, -13, 0],
            vec![4, 1, 0, 0],
            vec![0, 8, 1, 0],
            vec![5, -1, 0, 1],
        ]
        .into();
        let ans1 = vec![
            vec![2, 0, -13, 0],
            vec![0, 1, 26, 0],
            vec![-3, -3, 0, 1],
            vec![0, 6, -51, 0],
        ]
        .into();
        let id = Matrix::identity(m1.rows);
        assert_eq!(row_echelon_form(&m1, &id, 3).0, ans1);
        let m2: Matrix = vec![vec![4, 0, -3, 0], vec![4, 1, 0, 5], vec![0, 8, 1, 0]].into();
        let ans2: Matrix = vec![vec![4, 0, -3, 0], vec![0, 1, 3, 5], vec![0, 6, -5, -10]].into();
        assert_eq!(row_echelon_form(&m2, &id, 3).0, ans2);
    }

    #[test]
    fn test_binom() {
        assert_eq!(binom(6, 4), 15);
        assert_eq!(binom(-1, 3), 0);
        assert_eq!(binom(4, 0), 1);
        assert_eq!(binom(4, 7), 0);
    }

    #[test]
    fn binomial_symmetry() {
        for n in 0..20 {
            for m in 0..=n {
                assert_eq!(binom(n, m), binom(n, n - m));
            }
        }
    }

    #[test]
    fn test_extended_euclid() {
        assert_eq!(extended_euclid(24, 18), (6, 1, -1));
        assert_eq!(extended_euclid(54, 36), (18, 1, -1));
        assert_eq!(extended_euclid(120, 428860), (20, 3574, -1));
        assert_eq!(extended_euclid(95642, 1681), (1, 682, -38803));
    }

    #[test]
    #[should_panic]
    fn test_multiple_sols() {
        let v = vec![(0, 1), (0, 2), (0, 3), (0, 10)];
        knapsack_sols_exists(&v, 13);
    }

    #[test]
    fn test_knapsack_sols_exists() {
        let v = vec![(0, 1), (0, 2), (0, 3), (0, 10)];
        let w = vec![(0, 10), (1, 3), (2, 2), (3, 2)];
        let sum = 13;
        assert_eq!(knapsack_sols_exists(&w, sum), true);

        let w: Vec<(i64, i64)> = vec![];

        assert_eq!(knapsack_sols_exists(&w, sum), false);

        let sum = 0;

        assert_eq!(knapsack_sols_exists(&v, sum), true);

        assert_eq!(knapsack_sols_exists(&w, sum), true);
    }

    #[test]
    fn test_knapsack_sols() {
        let v = vec![(0, 1), (0, 3), (0, 10)];
        let sum = 13;
        assert_eq!(knapsack_sols(&v, sum), vec![(0, 3), (0, 10)]);

        let w: Vec<(i64, i64)> = vec![];

        assert_eq!(knapsack_sols(&w, sum), w);
        let sum = 0;

        assert_eq!(knapsack_sols(&v, sum), w);

        assert_eq!(knapsack_sols(&w, sum), w);
    }
}
