use std::io::{self,BufRead};


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

fn p_adic_val(mut a: i64, p: i64) -> i64 {
    assert!(a != 0, "vp(0)=infinity");
    let mut index = 0;
    while a % p == 0 {
        a = a / p;
        index += 1;
    }
    return index;
}

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

fn wedge_greater_than(x: i64, y: i64, p: i64) -> bool {
    let mut ans = true;
    let n: Vec<i64> = convert_base_p(x, p);
    let r: Vec<i64> = convert_base_p(y, p);
    for i in 0..n.len() {
        if r[i] > n[i] {
            ans = false;
            break;
        }
    }
    return ans;
}

fn e_tilde(n: i64, m: i64, p: i64) -> i64 {
    let ans: i64;
    let mut g = convert_base_p((n + m) / 2, p);
    let mut h = convert_base_p(m, p);
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
        && g[(p_adic_val(m, p) ) as usize] == h[(p_adic_val(m, p) ) as usize]
    {
        ans = 1;
    } else {
        ans = 0;
    }
    return ans;
}

fn binom(n: i64, k: i64) -> i64 {
    if k<0{ return 0;}
    let mut res = 1;
    for i in 0..k {
        res = (res * (n - i)) / (i + 1);
    }
    return res;
}

fn dimension(n: i64, m: i64, p: i64) -> i64 {
    let mut ans = 0;
    for r in 0..((n - m) / 2)+1 {
        ans += e_tilde(n - 2 * r + 1, m + 1, p) * (binom(n, r) - binom(n, r - 1));
    }
    return ans;
}

fn knapsack_count_sols(v: &Vec<(i64, i64)>, i: usize, sum: i64) -> i64 {
    if i >= v.len() {
        if sum == 0 {
            return 1;
        } else {
            return 0;
        }
    }
  return knapsack_count_sols(&v, i + 1, sum) + knapsack_count_sols(v, i + 1, sum - v[i].1);

}

fn knapsack_sols(v: &Vec<(i64, i64)>, mut sum: i64) -> Vec<(i64, i64)> {
    let mut ans: Vec<(i64, i64)> = Vec::new();
    for i in v {
        for x in v {
            if knapsack_count_sols(v, (i.1 + 1) as usize, sum - x.1) > 0 {
                ans.push(*x);
                sum -= x.1;
            }
        }
    }
    return ans;
}

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
    let numbers: Vec<i64> = 
        reader.lock()                           
              .lines().next().unwrap().unwrap() 
              .split(' ').map(|s| s.trim())     
              .filter(|s| !s.is_empty())        
              .map(|s| s.parse().unwrap())      
              .collect();                       
    
    


    println!("simples \n{:?}", find_simples(numbers[0] , numbers[1], numbers[2], numbers[3] ));

    //println!("e \n{:?}", e_tilde(5, 3, 2));
    //println!("e \n{:?}", (p_adic_val(3,2)-1)as usize);
    

    
   // println!("dim \n{:?}", dimension(12, i, 2));}}
}
