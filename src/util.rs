pub fn binom(n: i64, k: i64) -> i64 {
    if k < 0 || n < 0 {
        return 0;
    }
    let mut res = 1;
    for i in 0..k {
        res = (res * (n - i)) / (i + 1);
    }
    return res;
}

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
    }
}
