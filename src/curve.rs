//! Curve y^2 = x^3 + 3

use std::{ops::Add, fmt};

use super::field::Field;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Point {
    Inf,
    Reg { x: Field<101>, y: Field<101> },
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Point::Inf => write!(f, "inf"),
            Point::Reg { x, y } => write!(f, "{{ {x}, {y} }}"),
        }
    }
}

impl Point {
    pub const G: Self = Point::Reg {
        x: Field::new(1),
        y: Field::new(2),
    };

    /// Validating constructor
    pub fn new(x: Field<101>, y: Field<101>) -> Option<Point> {
        if y * y == x * x * x + 3.into() {
            Some(Point::Reg { x, y })
        } else {
            None
        }
    }

    pub fn inv(self) -> Self {
        match self {
            Point::Inf => Point::Inf,
            Point::Reg { x, y } => Point::Reg { x, y: -y },
        }
    }

    pub fn double(self) -> Self {
        let Point::Reg { x, y } = self else {
            return Point::Inf
        };
        let m = (Field::from(3) * x * x) / (Field::from(2) * y);
        Point::Reg {
            x: m * m - Field::from(2) * x,
            y: m * (Field::from(3) * x - m * m) - y,
        }
    }

    pub fn scalar_mul(self, s: Field<101>) -> Self {
        let mut x = self;
        let mut r = Self::Inf;
        for i in 0..7 {
            if s.check_bit(1 << i) {
                r = r + x;
            }
            x = x.double();
        }
        r
    }

    pub const CHEAT_TABLE: [Self; 17] = [
        Self::Inf,
        Self::Reg {
            x: Field::new(1),
            y: Field::new(2),
        },
        Self::Reg {
            x: Field::new(68),
            y: Field::new(74),
        },
        Self::Reg {
            x: Field::new(26),
            y: Field::new(45),
        },
        Self::Reg {
            x: Field::new(65),
            y: Field::new(98),
        },
        Self::Reg {
            x: Field::new(12),
            y: Field::new(32),
        },
        Self::Reg {
            x: Field::new(32),
            y: Field::new(42),
        },
        Self::Reg {
            x: Field::new(91),
            y: Field::new(35),
        },
        Self::Reg {
            x: Field::new(18),
            y: Field::new(49),
        },
        Self::Reg {
            x: Field::new(18),
            y: Field::new(52),
        },
        Self::Reg {
            x: Field::new(91),
            y: Field::new(66),
        },
        Self::Reg {
            x: Field::new(32),
            y: Field::new(59),
        },
        Self::Reg {
            x: Field::new(12),
            y: Field::new(69),
        },
        Self::Reg {
            x: Field::new(65),
            y: Field::new(3),
        },
        Self::Reg {
            x: Field::new(26),
            y: Field::new(56),
        },
        Self::Reg {
            x: Field::new(68),
            y: Field::new(27),
        },
        Self::Reg {
            x: Field::new(1),
            y: Field::new(99),
        },
    ];
}

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, rhs: Point) -> Self::Output {
        if self == rhs {
            return self.double();
        }
        let Point::Reg { x: rx, y: ry } = rhs else {
            return self;
        };
        let Point::Reg { x: lx, y: ly } = self else {
            return rhs;
        };
        if lx == rx {
            return Point::Inf;
        }

        let l = (ry - ly) / (rx - lx);
        let x = l * l - lx - rx;
        Point::Reg {
            x,
            y: l * (lx - x) - ly,
        }
    }
}

#[cfg(test)]
#[test]
fn doubling() {
    let mut x = Point::G;
    for _ in 0..5 {
        println!("{x}");
        x = x.double();
    }
}

#[cfg(test)]
#[test]
fn table() {
    let mut x = Point::Inf;
    for i in 0..17 {
        println!("{x}");
        assert_eq!(x, Point::G.scalar_mul(Field::new(i)));
        x = x + Point::G;
    }
    assert_eq!(x, Point::Inf);
    assert_eq!(Point::G.scalar_mul(Field::new(17)), Point::Inf);
}
