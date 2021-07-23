use std::io::{self, BufRead};
use tl_jantzen::binom;

/// Converts an i64 into the vector of its digits in base p
///
/// Returns the list of values of 'r' in base 'p'
/// mutiplied by the corresponding power of p. The values are 
/// given in descending order  
fn convert_base_p(mut r: i64, p: i64) -> Vec<i64> {
    assert!(2 <= p);
    let mut digits: Vec<i64> = Vec::new();
    let mut index = 0;
    while r > 0 {
        digits.push((r % p) * (p.pow(index as u32)));
        r = r / p;
        index += 1;
    }
    digits.reverse();
    return digits;
}

/// Calculates the 'p'-adic valuation of 'a' 
///
/// Returns the highest power of 'p' that divides
/// 'a'. Asserts that 'a' is not 0.
fn p_adic_val(mut a: i64, p: i64) -> i64 {
    assert!(a != 0, "vp(0)=infinity");
    let mut index = 0;
    while a % p == 0 {
        a = a / p;
        index += 1;
    }
    return index;
}

/// Calculates the value of Supp('r') with provided prime 'p'
///
/// Returns the vector of all values in the set Supp('r')
fn supp(r: i64, p: i64) -> Vec<i64> {
    let digits: Vec<i64> = convert_base_p(r + 1, p);
    let mut set: Vec<i64> = Vec::new();
    set.push(digits[0]);
    let mut new_set: Vec<i64>;
    for i in 1..digits.len() {
        new_set = Vec::new();
        for x in set.iter() {
            new_set.push(x + digits[i]);
            new_set.push(x - digits[i]);
        }
        set = new_set.clone();
    }
    for i in 0..set.len() {
        set[i] -= 1;
    }
    set.sort_unstable();
    set.dedup();
    return set;
}

/// Determines whether x is wedge-greater-than y
///
/// Returns TRUE if x is greater than y under the wedge
// ordering and FALSE otherwise
fn wedge_greater_than(x: i64, y: i64, p: i64) -> bool {
    let mut g: Vec<i64> = convert_base_p(x, p);
    let mut h: Vec<i64> = convert_base_p(y, p);
    if g.len() > h.len() {
        let mut r = vec![0; g.len() - h.len()];
        r.append(&mut h);
        h = r;
    }
    if g.len() < h.len() {
        let mut r = vec![0; h.len() - g.len()];
        r.append(&mut g);
        g = r;
    }
    for i in 0..g.len() {
        if h[i] > g[i] {
            return false;
        }
    }
    return true;
}

/// Calculates the value e-tilde(n,m)
///
/// Returns -1, 0, or 1 according to the formula
/// for e-tilde
fn e_tilde(n: i64, m: i64, p: i64) -> i64 {
    let ans: i64;
    let mut g = convert_base_p((n + m) / 2, p);
    let mut h = convert_base_p(m, p);
    if g.len() > h.len() {
        let mut r = vec![0; g.len() - h.len()];
        r.append(&mut h);
        h = r;
    }
    if g.len() < h.len() {
        let mut r = vec![0; h.len() - g.len()];
        r.append(&mut g);
        g = r;
    }

    g.reverse();
    h.reverse();
    if (n - m) % 2 != 0 {
        ans = 0;
    } else if p_adic_val((n + m) / 2, p) > p_adic_val(m, p)
        && wedge_greater_than(((n + m) / 2) - 1, m, p)
    {
        ans = -1;
    } else if p_adic_val((n + m) / 2, p) == p_adic_val(m, p)
        && wedge_greater_than((n + m) / 2, m, p)
        && g[(p_adic_val(m, p)) as usize] == h[(p_adic_val(m, p)) as usize]
    {
        ans = 1;
    } else {
        ans = 0;
    }
    return ans;
}

/// Calculates the dimension of the simple module D(n,m)
///
/// Returns the integer dimension of the simple module
/// indexed by 'n' and 'm' in characteristic 'p'
fn dimension(n: i64, m: i64, p: i64) -> i64 {
    let mut ans = 0;
    for r in 0..((n - m) / 2) + 1 {
        ans += e_tilde(n - 2 * r + 1, m + 1, p) * (binom(n, r) - binom(n, r - 1));
    }
    return ans;
}

/// Checks for a solution to the knapsack problem
///
/// Accepts a vector of pairs and checks if the sum of 
/// the second entries of a subset of these pairs totals 
/// to 'sum'. Returns TRUE if a solutions exists and 
/// FALSE otherwise
fn knapsack_sols_exists(v: &[(i64, i64)], sum: i64) -> bool {
    if v.is_empty() {
        sum == 0
    } else {
        knapsack_sols_exists(&v[1..], sum) | knapsack_sols_exists(&v[1..], sum - v[0].1)
    }
}

/// Returns a solution to the knapsack problem if it exists
///
/// Accepts a vector of pairs and returns a subset of these pairs if 
/// the sum of their second entry totals to 'sum'
fn knapsack_sols(v: &Vec<(i64, i64)>, mut sum: i64) -> Vec<(i64, i64)> {
    let mut ans: Vec<(i64, i64)> = Vec::new();
    for (n, x) in v.iter().enumerate() {
        if knapsack_sols_exists(&v[(n + 1)..], sum - x.1) {
            ans.push(*x);
            sum -= x.1;
        }
    }
    return ans;
}

/// Finds the simple modules in a Jantzen layer of a cell module of TLn
///
/// Returns the indicies corresponding to simple modules in the Jantzen 
/// layer of dimension 'dim'
fn find_simples(n: i64, m: i64, p: i64, dim: i64) -> Vec<(i64, i64)> {
    let mut indicies: Vec<i64> = Vec::new();
    for r in m..n + 1 {
        if (n - r) % 2 == 0 {
            let support: Vec<i64> = supp(r, p);
            if support.contains(&m) {
                indicies.push(r);
            }
        }
    }

    let mut simples: Vec<(i64, i64)> = Vec::new();
    for k in &indicies {
        simples.push((*k, dimension(n, *k, p)));
    }

    let ans: Vec<(i64, i64)> = knapsack_sols(&simples, dim);

    return ans;
}

fn main() {
    println!("Enter m n p dim seperated by spaces");
    let reader = io::stdin();
    let numbers: Vec<i64> = reader
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

    println!(
        "simples \n{:?}",
        find_simples(numbers[0], numbers[1], numbers[2], numbers[3])
    );
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_convert_base_p() {
        let v = vec![4, 2, 1];
        assert_eq!(convert_base_p(7, 2), v);
    }

    #[test]
    fn test_p_adic_val() {
        assert_eq!(p_adic_val(54, 3), 3);
    }

    #[test]
    fn test_supp() {
        let v = vec![0, 2, 4, 6];
        assert_eq!(supp(6, 2), v);
    }

    #[test]
    fn test_wedge_greater_than() {
        assert_eq!(wedge_greater_than(8, 5, 3), true);
    }

    #[test]
    fn test_e_tilde() {
        let e = vec![
            vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            vec![0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            vec![0, -1, 0, 1, 0, 0, 0, 0, 0, 0],
            vec![-1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            vec![0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            vec![1, 0, 0, 0, -1, 0, 1, 0, 0, 0],
            vec![0, 1, 0, -1, 0, 0, 0, 1, 0, 0],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            vec![0, -1, 0, 1, 0, 0, 0, -1, 0, 1],
        ];
        for i in 1..e.len() {
            for j in 1..e.len() {
                assert_eq!(e[i - 1][j - 1], e_tilde(i as i64, j as i64, 3));
            }
        }
    }

    #[test]
    fn test_dimension() {
        assert_eq!(dimension(6, 2, 2), 4);
    }

    #[test]
    fn test_knapsack_sols_exists() {
        let v = vec![(0, 1), (0, 2), (0, 3), (0, 10)];
        let w = vec![(0, 10), (1, 3), (2, 2), (3, 2)];
        let sum = 13;
        assert_eq!(knapsack_sols_exists(&w, sum), true);
    }

    #[test]
    fn test_knapsack_sols() {
        let v = vec![(0, 1), (0, 3), (0, 10)];
        let sum = 13;
        assert_eq!(knapsack_sols(&v, sum), vec![(0, 3), (0, 10)]);
    }
}
