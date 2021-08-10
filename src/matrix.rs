/// A matrix type that holds i64
///
/// This matrix is not required to be square.
///
/// This type is preferred to alternatives such as `Vec<Vec<i64>>`
/// because it is assured each row is the same length.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    data: Vec<i64>,
}

impl Matrix {
    /// Construct a new, zeroed matrix of the given shape
    pub fn new(rows: usize, cols: usize) -> Matrix {
        Matrix {
            rows,
            cols,
            data: vec![0; rows * cols],
        }
    }

    /// Construct an n x n identity matrix
    pub fn identity(n: usize) -> Matrix {
        let mut ans = Matrix::new(n, n);
        for i in 0..n {
            ans[(i, i)] = 1;
        }
        ans
    }

    /// Swap two rows of this matrix in-place
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        // TODO (robert) : this is inefficient again.  See about using
        // compiler built in swap.
        for k in 0..self.cols {
            self.data.swap(i * self.cols + k, j * self.cols + k);
        }
    }

    /// Swap two columns of this matrix in-place
    pub fn swap_cols(&mut self, i: usize, j: usize) {
        // TODO (robert) : this is inefficient again.  See about using
        // compiler built in swap.
        for k in 0..self.rows {
            self.data.swap(k * self.cols + i, k * self.cols + j);
        }
    }

    /// Adds n times row i to row j to acheive a new row j
    pub fn add_multiple_rowi_to_rowj(&mut self, n: i64, i: usize, j: usize) {
        for k in 0..self.cols {
            self.data[j * self.cols + k] += n * self.data[i * self.cols + k];
        }
    }

    /// Adds n times column i to col jumn to acheive a new col j
    pub fn add_multiple_coli_to_colj(&mut self, n: i64, i: usize, j: usize) {
        for k in 0..self.rows {
            self.data[k * self.cols + j] += n * self.data[k * self.cols + i];
        }
    }

    /// Get a list of slices into the rows
    pub fn row_list(&self) -> Vec<&[i64]> {
        let mut ans = Vec::new();
        for i in 0..self.rows {
            ans.push(&self.data[i * self.cols..(i + 1) * self.cols]);
        }
        ans
    }

    pub fn col(&self, j: usize) -> Vec<i64> {
        let mut col = Vec::new();
        for i in 0..self.rows {
            col.push(self[(i, j)]);
        }
        col
    }

    pub fn row(&self, i: usize) -> Vec<i64> {
        let mut row = Vec::new();
        for j in 0..self.cols {
            row.push(self[(i, j)]);
        }
        row
    }

    pub fn entries(&self) -> &[i64] {
        &self.data
    }

    pub fn add_row(&self, v: Vec<i64>) -> Matrix {
        let mut ans: Matrix;
        if self.data.len() == 0 {
            ans = Matrix::new(1, v.len());
            for i in 0..v.len() {
                ans[(0, i)] = v[i];
            }
        } else {
            assert!(v.len() == self.cols);
            ans = Matrix::new(self.rows + 1, self.cols);
            for i in 0..self.rows + 1 {
                for j in 0..self.cols {
                    if i < self.rows {
                        ans[(i, j)] = self.data[i * self.cols + j];
                    } else {
                        ans[(i, j)] = v[j];
                    }
                }
            }
        }
        return ans;
    }
}

impl core::ops::Index<(usize, usize)> for Matrix {
    type Output = i64;

    fn index(&self, index: (usize, usize)) -> &i64 {
        assert!(index.0 < self.rows);
        assert!(index.1 < self.cols);
        &self.data[index.0 * self.cols + index.1]
    }
}

impl core::ops::IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut i64 {
        assert!(index.0 < self.rows);
        assert!(index.1 < self.cols);
        &mut self.data[index.0 * self.cols + index.1]
    }
}

impl core::convert::From<Vec<Vec<i64>>> for Matrix {
    fn from(other: Vec<Vec<i64>>) -> Matrix {
        let rows = other.len();
        assert!(rows > 0);
        let cols = other[0].len();
        let mut ans = Matrix::new(rows, cols);
        for row in 0..rows {
            assert_eq!(cols, other[row].len());
            for col in 0..cols {
                ans[(row, col)] = other[row][col];
            }
        }
        ans
    }
}

impl core::convert::From<&[&[i64]]> for Matrix {
    fn from(other: &[&[i64]]) -> Matrix {
        let rows = other.len();
        assert!(rows > 0);
        let cols = other[0].len();
        let mut ans = Matrix::new(rows, cols);
        for row in 0..rows {
            assert_eq!(cols, other[row].len());
            for col in 0..cols {
                ans[(row, col)] = other[row][col];
            }
        }
        ans
    }
}

impl core::ops::Index<usize> for Matrix {
    type Output = [i64];

    fn index(&self, index: usize) -> &[i64] {
        assert!(index < self.rows);
        &self.data[index * self.cols..(index + 1) * self.cols]
    }
}

impl core::ops::Add<Matrix> for Matrix {
    type Output = Matrix;

    fn add(self, rhs: Matrix) -> Matrix {
        assert_eq!(
            self.rows, rhs.rows,
            "Cannot add matrices of different sizes"
        );
        assert_eq!(
            self.cols, rhs.cols,
            "Cannot add matrices of different sizes"
        );
        Matrix {
            data: self.data.iter().zip(rhs.data).map(|(x, y)| x + y).collect(),
            ..self
        }
    }
}

impl core::ops::Sub<Matrix> for Matrix {
    type Output = Matrix;

    fn sub(self, rhs: Matrix) -> Matrix {
        assert_eq!(
            self.rows, rhs.rows,
            "Cannot subtract matrices of different sizes"
        );
        assert_eq!(
            self.cols, rhs.cols,
            "Cannot subtract matrices of different sizes"
        );
        Matrix {
            data: self.data.iter().zip(rhs.data).map(|(x, y)| x - y).collect(),
            ..self
        }
    }
}

impl core::ops::Mul<&Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: &Matrix) -> Matrix {
        &self * rhs
    }
}

impl core::ops::Mul<Matrix> for &Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Matrix) -> Matrix {
        self * &rhs
    }
}

impl core::ops::Mul<&Matrix> for &Matrix {
    type Output = Matrix;

    fn mul(self, rhs: &Matrix) -> Matrix {
        assert_eq!(self.cols, rhs.rows, "Cannot multiply mismatched matrices");
        let mut ans = Matrix::new(self.rows, rhs.cols);
        for i in 0..self.rows {
            for j in 0..rhs.cols {
                for k in 0..self.cols {
                    ans[(i, j)] += self[(i, k)] * rhs[(k, j)]
                }
            }
        }
        ans
    }
}

impl core::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for row in 0..self.rows {
            write!(f, "[")?;
            for col in 0..self.cols {
                write!(f, " {:2}", self[(row, col)])?;
            }
            write!(f, "]")?;
            if row < self.rows - 1 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn interesting_matrix(rows: usize, cols: usize) -> Matrix {
        let mut matrix = Matrix::new(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                matrix[(r, c)] = (r * cols + c) as i64;
            }
        }
        matrix
    }
    #[test]
    fn new_matrix() {
        let rows = 10;
        let cols = 12;
        let matrix = Matrix::new(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                assert_eq!(matrix[(r, c)], 0);
            }
        }
    }

    #[test]
    fn identity_matrix() {
        let rows = 10;
        let matrix = Matrix::identity(rows);
        for r in 0..rows {
            for c in 0..rows {
                if r == c {
                    assert_eq!(matrix[(r, c)], 1);
                } else {
                    assert_eq!(matrix[(r, c)], 0);
                }
            }
        }
    }

    #[test]
    fn test_add_row() {
        let mut matrix = Matrix::identity(3);
        let v = vec![9, 1, -3];
        matrix = matrix.add_row(v);
        let ans: Matrix = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1], vec![9, 1, -3]].into();
        assert_eq!(matrix, ans);
        let v = vec![9, 1, -3];
        let mut matrix = Matrix::new(0, 0);
        let ans: Matrix = vec![vec![9, 1, -3]].into();
        matrix = matrix.add_row(v);
        assert_eq!(matrix, ans);
    }

    #[test]
    fn swap_matrix_rows() {
        let rows = 10;
        let cols = 8;
        let mut matrix = interesting_matrix(rows, cols);

        matrix.swap_rows(4, 7);

        for r in 0..rows {
            for c in 0..cols {
                match r {
                    4 => assert_eq!(matrix[(r, c)] as usize, 7 * cols + c),
                    7 => assert_eq!(matrix[(r, c)] as usize, 4 * cols + c),
                    r => assert_eq!(matrix[(r, c)] as usize, r * cols + c),
                }
            }
        }
    }

    #[test]
    fn test_col() {
        let mut mat: Matrix = vec![vec![1, 0, -1], vec![0, 1, 0], vec![-3, 0, 1]].into();

        assert_eq!(mat.col(0), vec![1, 0, -3]);

        assert_eq!(mat.col(2), vec![-1, 0, 1]);
    }

    #[test]
    fn test_row() {
        let mut mat: Matrix = vec![vec![1, 0, -1], vec![0, 1, 0], vec![-3, 0, 1]].into();

        assert_eq!(mat.row(0), vec![1, 0, -1]);

        assert_eq!(mat.row(2), vec![-3, 0, 1]);
    }

    #[test]
    #[should_panic]
    fn return_out_of_bounds_cols() {
        let mut mat: Matrix = vec![vec![1, 0, -1], vec![0, 1, 0], vec![-3, 0, 1]].into();
        mat.col(3);
        mat.row(4);
    }

    #[test]
    fn swap_matrix_cols() {
        let mut mat: Matrix = vec![vec![1, 0, -1], vec![0, 1, 0], vec![-3, 0, 1]].into();
        let ans: Matrix = vec![vec![-1, 0, 1], vec![0, 1, 0], vec![1, 0, -3]].into();
        mat.swap_cols(0, 2);
        assert_eq!(mat, ans);
    }

    #[test]
    fn test_add_rows() {
        let mut mat: Matrix = Matrix::identity(3);
        let ans: Matrix = vec![vec![1, 0, 0], vec![3, 1, 0], vec![0, 0, 1]].into();
        mat.add_multiple_rowi_to_rowj(3, 0, 1);
        assert_eq!(mat, ans);
    }

    #[test]
    fn test_add_cols() {
        let mut mat: Matrix = Matrix::identity(3);
        let ans: Matrix = vec![vec![1, 0, 0], vec![-2, 1, 0], vec![0, 0, 1]].into();
        mat.add_multiple_coli_to_colj(-2, 1, 0);
        assert_eq!(mat, ans);
    }

    #[test]
    fn matrix_row_list() {
        let rows = 10;
        let cols = 8;
        let matrix = interesting_matrix(rows, cols);

        let list = matrix.row_list();

        for (r, row) in list.iter().enumerate() {
            for c in 0..cols {
                assert_eq!(row[c] as usize, r * cols + c);
            }
        }
    }

    #[test]
    fn matrix_indexing() {
        let rows = 10;
        let cols = 8;
        let matrix = interesting_matrix(rows, cols);

        for r in 0..rows {
            for c in 0..cols {
                assert_eq!(matrix[r][c], matrix[(r, c)]);
            }
        }
    }

    #[test]
    #[rustfmt::skip]
    fn matrix_add() {
        let a : Matrix = vec![
            vec![1,2,3],
            vec![-1,-1,0],
            vec![3,8,9],

        ].into();
        let b : Matrix = vec![
            vec![4,2,4],
            vec![1,5,1],
            vec![6,2,9],
        ].into();
        let c : Matrix = vec![
            vec![5,4,7],
            vec![0,4,1],
            vec![9,10,18],

        ].into();
            assert_eq!(a + b, c);
    }

    #[test]
    #[rustfmt::skip]
    fn matrix_sub() {
        let a : Matrix = vec![
            vec![1,2,3],
            vec![-1,-1,0],
            vec![3,8,9],

        ].into();
        let b : Matrix = vec![
            vec![4,2,4],
            vec![1,5,1],
            vec![6,2,9],
        ].into();
        let c : Matrix = vec![
            vec![-3,0,-1],
            vec![-2,-6,-1],
            vec![-3,6,0],

        ].into();
            assert_eq!(a - b, c);
    }

    #[test]
    #[rustfmt::skip]
    fn matrix_mul() {
        let a : Matrix = vec![
            vec![1,2,3],
            vec![-1,-1,0],
            vec![3,8,9],

        ].into();
        let b : Matrix = vec![
            vec![4,2,4],
            vec![1,5,1],
            vec![6,2,9],
        ].into();
        let c : Matrix = vec![
            vec![24,18,33],
            vec![-5,-7,-5],
            vec![74,64,101],

        ].into();
        assert_eq!(&a * &b, c);
        assert_eq!(&a * b.clone(), c);
        assert_eq!(a.clone() * &b, c);
        assert_eq!(a * &b, c);
    }
}
