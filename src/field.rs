//! The Finite Field F_101

use std::{
    ops::{Add, Sub, Neg, Mul, Div, AddAssign},
    fmt,
    iter::Sum,
};

use rand::distributions::{Distribution, Standard};

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Field<const M: u8>(u8);

impl<const M: u8> Distribution<Field<M>> for Standard {
    fn sample<R>(&self, rng: &mut R) -> Field<M>
    where
        R: rand::Rng + ?Sized,
    {
        loop {
            let x = rng.gen::<u8>();
            if x < M {
                break Field::new(x);
            }
        }
    }
}

impl<const M: u8> fmt::Display for Field<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const M: u8> Field<M> {
    pub const fn new(x: u8) -> Self {
        if x < M {
            Self(x)
        } else {
            panic!()
        }
    }

    pub const fn is_zero(&self) -> bool {
        self.0 == 0
    }

    pub const fn check_bit(self, other: u8) -> bool {
        (self.0 & other) != 0
    }

    pub const fn cast<const K: u8>(self) -> Field<K> {
        Field::new(self.0)
    }

    pub fn power(self, n: usize) -> Self {
        let mut r = Self::new(1);
        for _ in 0..n {
            r = r * self;
        }
        r
    }
}

impl<const M: u8> From<u8> for Field<M> {
    fn from(value: u8) -> Self {
        Field(value % M)
    }
}

impl<const M: u8> Neg for Field<M> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.0 == 0 {
            self
        } else {
            Field(M - self.0)
        }
    }
}

impl<const M: u8> Add for Field<M> {
    type Output = Self;

    fn add(self, rhs: Field<M>) -> Self::Output {
        Field((self.0 + rhs.0) % M)
    }
}

impl<const M: u8> AddAssign for Field<M> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 = (self.0 + rhs.0) % M
    }
}

impl<const M: u8> Sum for Field<M> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut s = Self::default();
        for item in iter {
            s += item;
        }
        s
    }
}

impl<const M: u8> Sub for Field<M> {
    type Output = Self;

    fn sub(self, rhs: Field<M>) -> Self::Output {
        self + (-rhs)
    }
}

impl<const M: u8> Mul for Field<M> {
    type Output = Self;

    fn mul(self, rhs: Field<M>) -> Self::Output {
        let x = (self.0 as u16 * rhs.0 as u16) % (M as u16);
        Field(x as u8)
    }
}

impl<const M: u8> Div for Field<M> {
    type Output = Self;

    fn div(self, rhs: Field<M>) -> Self::Output {
        for x in 0..M {
            let x = Self::from(x);
            if rhs * x == self {
                return x;
            }
        }

        unreachable!()
    }
}
