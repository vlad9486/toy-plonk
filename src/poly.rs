use std::{
    fmt,
    ops::{Add, Mul, Sub},
};

use super::{field::Field, matrix, srs::Srs, curve::Point};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Poly<const N: usize>([Field<17>; N]);

impl<const N: usize> From<Field<17>> for Poly<N> {
    fn from(value: Field<17>) -> Self {
        Poly::new([value]).resize()
    }
}

impl<const N: usize> Poly<N> {
    // length of v must be N
    pub fn new<T>(a: [T; N]) -> Self
    where
        T: Into<Field<17>>,
    {
        let mut s = Poly([Field::new(0); N]);
        for (x, a) in s.0.iter_mut().zip(a) {
            *x = a.into();
        }
        s
    }

    // length of v must be N
    pub fn interpolate_vector<T>(v: [T; N]) -> Self
    where
        T: Into<Field<17>> + Clone,
    {
        match N {
            4 => {
                let h_inv = matrix::h_inv_matrix();
                let mut s = [Field::new(0); N];
                for (row, res) in h_inv.into_iter().zip(&mut s) {
                    *res = row
                        .iter()
                        .zip(v.clone())
                        .fold(Field::new(0), |a, (m, input)| a + *m * input.into());
                }
                Self(s)
            }
            _ => unimplemented!(),
        }
    }

    pub fn evaluate(&self, at: Field<17>) -> Field<17> {
        let mut t = Field::new(1);
        let mut r = Field::new(0);
        for x in self.0 {
            r = r + t * x;
            t = t * at;
        }
        r
    }

    pub fn evaluate_srs(&self) -> Point {
        let srs = Srs::new(2.into());
        srs.zip(&self.0).fold(Point::Inf, |a, (sr, coefficient)| {
            a + sr.scalar_mul(coefficient.cast())
        })
    }

    pub fn substitute(mut self, coefficient: Field<17>) -> Self {
        let mut t = Field::new(1);
        for x in &mut self.0 {
            *x = *x * t;
            t = t * coefficient;
        }
        self
    }

    pub fn divide_by_h4(self) -> Self {
        // h4 = x^4 - 1
        let mut t = Self::new([0; N]);
        for (i, s) in self.0.iter().enumerate().rev() {
            if i >= N - 4 {
                t.0[i - 4] = *s;
            } else if i >= 4 {
                t.0[i - 4] = *s + t.0[i];
            }
        }
        t
    }

    pub fn split<const K: usize>(self) -> (Poly<K>, Poly<K>, Poly<K>) {
        let mut lo = Poly([0.into(); K]);
        let mut mid = Poly([0.into(); K]);
        let mut hi = Poly([0.into(); K]);
        for i in 0..K {
            lo.0[i] = self.0[i];
            mid.0[i] = self.0[i + K];
            hi.0[i] = self.0[i + 2 * K];
        }
        (lo, mid, hi)
    }

    // self / (x - x0)
    pub fn divide_by_first(self, x0: Field<17>) -> Self {
        let mut t = Self::new([0; N]);
        for (i, s) in self.0.iter().enumerate().rev() {
            if i == N - 1 {
                t.0[i - 1] = *s;
            } else if i >= 1 {
                t.0[i - 1] = *s + x0 * t.0[i];
            } else {
                assert_eq!(*s + x0 * t.0[i], Field::new(0));
            }
        }
        t
    }
}

impl<const N: usize> Poly<N> {
    pub fn resize<const M: usize>(self) -> Poly<M> {
        let mut s = Poly([Field::new(0); M]);
        for i in 0..N.min(M) {
            s.0[i] = self.0[i];
        }
        s
    }
}

impl<const N: usize> Add<Poly<N>> for Poly<N> {
    type Output = Poly<N>;

    fn add(self, rhs: Poly<N>) -> Self::Output {
        let mut s = self;
        s.0.iter_mut()
            .zip(self.0.iter().zip(rhs.0.iter()))
            .for_each(|(c, (a, b))| {
                *c = *a + *b;
            });
        s
    }
}

impl<const N: usize> Sub<Poly<N>> for Poly<N> {
    type Output = Poly<N>;

    fn sub(self, rhs: Poly<N>) -> Self::Output {
        let mut s = self;
        s.0.iter_mut()
            .zip(self.0.iter().zip(rhs.0.iter()))
            .for_each(|(c, (a, b))| {
                *c = *a - *b;
            });
        s
    }
}

impl<const N: usize, const M: usize> Mul<Poly<M>> for Poly<N> {
    type Output = Poly<N>;

    fn mul(self, rhs: Poly<M>) -> Self::Output {
        let mut out = Poly([Field::new(0); N]);
        for (i, x) in out.0.iter_mut().enumerate() {
            *x = (0..(i + 1))
                .map(|j| self.0[j] * rhs.0[i - j])
                .sum::<Field<17>>();
        }
        out
    }
}

impl<const N: usize> Mul<Field<17>> for Poly<N> {
    type Output = Poly<N>;

    fn mul(self, rhs: Field<17>) -> Self::Output {
        let mut out = Poly([Field::new(0); N]);
        for (x, s) in out.0.iter_mut().zip(&self.0) {
            *x = *s * rhs;
        }
        out
    }
}

impl<const N: usize> fmt::Display for Poly<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for (n, x) in self.0.iter().enumerate().filter(|(_, x)| !x.is_zero()) {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            if n == 0 {
                write!(f, "{x}")?;
            } else if n == 1 {
                if Field::from(1) == *x {
                    write!(f, "x")?;
                } else {
                    write!(f, "{x} * x")?;
                }
            } else {
                if Field::from(1) == *x {
                    write!(f, "x ^ {n}")?;
                } else {
                    write!(f, "{x} * x ^ {n}")?;
                }
            }
        }
        if first {
            write!(f, "0")?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[test]
fn interpolate() {
    let p = Poly::<4>::interpolate_vector([3, 4, 5, 9]);
    assert_eq!(p.0, [1.into(), 13.into(), 3.into(), 3.into()]);
    println!("{p}");
}

#[cfg(test)]
#[test]
fn multiply() {
    let h = Poly::new([16, 0, 0, 0, 1]);
    let fa = Poly::new([1, 13, 3, 3]);
    let b = Poly::new([4, 7]);

    let x = (h.resize::<6>() * b.resize::<6>()) + fa.resize();
    assert_eq!(x, Poly::new([14, 6, 3, 3, 4, 7]));
    println!("{x}");
}
