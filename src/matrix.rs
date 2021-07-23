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

    /// Get a list of slices into the rows
    pub fn row_list(&self) -> Vec<&[i64]> {
        let mut ans = Vec::new();
        for i in 0..self.rows {
            ans.push(&self.data[i * self.cols..(i + 1) * self.cols]);
        }
        ans
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

impl core::ops::Index<usize> for Matrix {
    type Output = [i64];

    fn index(&self, index: usize) -> &[i64] {
        assert!(index < self.rows);
        &self.data[index * self.cols..(index + 1) * self.cols]
    }
}

impl core::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for row in 0..self.rows {
            write!(f, "[")?;
            for col in 0..self.cols {
                write!(f, " {:3}", self[(row, col)])?;
            }
            write!(f, "]\n")?;
        }
        Ok(())
    }
}