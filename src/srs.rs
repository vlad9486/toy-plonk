//! Structured reference string

use super::{field::Field, curve::Point};

pub struct Srs {
    this: Point,
    s: Field<101>,
}

impl Srs {
    pub const fn new(s: Field<101>) -> Self {
        Srs { this: Point::G, s }
    }
}

impl Iterator for Srs {
    type Item = Point;

    fn next(&mut self) -> Option<Self::Item> {
        let r = self.this;
        self.this = self.this.scalar_mul(self.s);
        Some(r)
    }
}

#[cfg(test)]
#[test]
fn srs_2() {
    Srs::new(2.into()).take(10).for_each(|p| println!("{p}"));
}
