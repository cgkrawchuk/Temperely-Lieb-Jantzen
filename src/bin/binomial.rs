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