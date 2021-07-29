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
pub fn mod_inv(a: i64, module: i64) -> i64 {
    assert!(a % module != 0, "number is 0 mod...");
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

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_mod_inv() {
        assert_eq!(mod_inv(3, 5), 2);
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
}
