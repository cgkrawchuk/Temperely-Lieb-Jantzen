use tl_jantzen::{gram_matrix, row_echelon_form, Matrix};

use std::env;

/// This fn needs to be modified...
fn recursive_ops(m: usize, n: usize, p: i64) {
    let mut g = gram_matrix(m, n);
    let mut id = Matrix::identity(g.rows);

    let mut v = Vec::new();
    v.push(0);
    let mut r = 0;

    while r < g.rows {
        println!("r is: {}", r);
        let values = row_echelon_form(&g, &id, p);
        g = values.0;
        id = values.1;
        let rank = values.2;

        v.push(rank);

        println!("basis of layer is: ");
        for i in r..rank {
            println!("{:?}", id.row(i));
        }

        for i in rank..g.rows {
            for j in 0..g.cols {
                g[(i, j)] /= p;
            }
        }
        for i in 0..rank {
            for j in 0..id.cols {
                id[(i, j)] *= p;
            }
        }

        r += v[v.len() - 1] - v[v.len() - 2];
    }

    for i in 1..v.len() {
        println!("{:?}", v[i] - v[i - 1]);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        println!("Enter m n p seperated by spaces");
        return;
    }
    let n: usize = args[1].parse().expect("Please enter an integer for n");
    let m: usize = args[2].parse().expect("Please enter an integer for m");
    let p: i64 = args[3].parse().expect("Please enter an integer for p");

    recursive_ops(n, m, p);
}
#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_gram_matrix() {
        let m: Matrix = vec![
            vec![4, 2, 0, 0, 2, 1, 0, 0, 2, 0, 0, 0, 0, 0],
            vec![2, 4, 2, 0, 1, 2, 1, 0, 1, 0, 0, 0, 0, 0],
            vec![0, 2, 4, 2, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0],
            vec![0, 0, 2, 4, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0],
            vec![2, 1, 0, 0, 4, 2, 0, 0, 1, 0, 0, 1, 0, 0],
            vec![1, 2, 1, 0, 2, 4, 2, 0, 2, 1, 0, 2, 0, 0],
            vec![0, 1, 2, 1, 0, 2, 4, 2, 1, 2, 1, 1, 0, 0],
            vec![0, 0, 1, 2, 0, 0, 2, 4, 0, 1, 2, 0, 0, 0],
            vec![2, 1, 0, 0, 1, 2, 1, 0, 4, 2, 0, 1, 0, 1],
            vec![0, 0, 0, 0, 0, 1, 2, 1, 2, 4, 2, 2, 1, 2],
            vec![0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 4, 1, 2, 1],
            vec![0, 0, 0, 0, 1, 2, 1, 0, 1, 2, 1, 4, 2, 1],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 4, 2],
            vec![0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 4],
        ]
        .into();
        assert_eq!(gram_matrix(7, 3), m);
    }
}
