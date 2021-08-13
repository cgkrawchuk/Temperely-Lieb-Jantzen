use tl_jantzen::{binom, convert_base_p, elem_div, gram_matrix, knapsack_sols, p_adic_val, Matrix};

use std::env;

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
    set
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
    true
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
    ans
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
    ans
}

pub fn find_layers(m: i64, n: i64, p: i64) {
    let mut g = gram_matrix(m as usize, n as usize);

    let dimensions = elem_div(&mut g, p);

    let mut indicies: Vec<i64> = Vec::new();
    for r in n..m + 1 {
        if (m - r) % 2 == 0 {
            let support: Vec<i64> = supp(r, p);
            if support.contains(&n) {
                indicies.push(r);
            }
        }
    }

    let mut simples: Vec<(i64, i64)> = Vec::new();
    for k in &indicies {
        simples.push((*k, dimension(m, *k, p)));
    }

    println!("simples: {:?}", simples);
    for y in dimensions {
        println!("layer has dimension: {:?}", y);

        let ans: Vec<(i64, i64)> = knapsack_sols(&simples, y);

        println!("simples are: {:?}", ans);

        for z in ans {
            let index = simples.iter().position(|x| *x == z).unwrap();
            simples.remove(index);
        }
    }
}

fn main() {
    /*
    for i in 16..20 {
        for j in (0..i).rev() {
            if (i - j) % 2 == 0 {
                println!("n: {0}, m: {1}, p: 2", i, j);
                find_layers(i, j, 2);
                println!("");
            }
        }
    }
    */

    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        println!("Enter m n p seperated by spaces");
        return;
    }
    let n: i64 = args[1].parse().expect("Please enter an integer for n");
    let m: i64 = args[2].parse().expect("Please enter an integer for m");
    let p: i64 = args[3].parse().expect("Please enter an integer for p");

    find_layers(n, m, p);
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_dimension() {
        assert_eq!(dimension(6, 2, 2), 4);
        assert_eq!(dimension(8, 4, 3), 13);
        assert_eq!(dimension(10, 4, 5), 75);
        assert_eq!(dimension(0, 0, 2), 1);
    }

    #[test]
    fn test_supp() {
        let v = vec![0, 2, 4, 6];
        assert_eq!(supp(6, 2), v);

        let v = vec![6, 8];

        assert_eq!(supp(8, 8), v);

        let v = vec![0];

        assert_eq!(supp(0, 4), v);

        for n in 0..100 {
            assert!(supp(n, 2).contains(&n));
        }
    }

    #[test]
    fn test_wedge_greater_than() {
        assert_eq!(wedge_greater_than(8, 5, 3), true);
        assert_eq!(wedge_greater_than(0, 0, 3), true);
        assert_eq!(wedge_greater_than(15, 8, 2), true);
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
    fn test_for_overflow() {
        find_layers(8, 2, 2);
        find_layers(10, 2, 2);
        find_layers(10, 6, 3);
    }
}
