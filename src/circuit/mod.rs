use std::fmt;

use crate::field::ExtendedField101;

use super::{field::Field, curve::Point, poly::Poly, matrix};

pub struct Assignment<const N: usize> {
    a: [u8; N],
    b: [u8; N],
    c: [u8; N],
}

pub struct Selector<const N: usize> {
    l: Poly<N>,
    r: Poly<N>,
    zero: Poly<N>,
    m: Poly<N>,
    c: Poly<N>,
}

pub struct CopyConstraints<const N: usize> {
    s1: Poly<N>,
    s2: Poly<N>,
    s3: Poly<N>,
}

pub struct Circuit<const N: usize> {
    assignments: Assignment<N>,
    selectors: Selector<N>,
    copy_constraints: CopyConstraints<N>,
}

pub struct Proof {
    a: Point,
    b: Point,
    c: Point,
    zet: Point,
    t_lo: Point,
    t_mid: Point,
    t_hi: Point,
    w_zeta: Point,
    w_zeta_omega: Point,
    a_val: Field<17>,
    b_val: Field<17>,
    c_val: Field<17>,
    s1_val: Field<17>,
    s2_val: Field<17>,
    r_val: Field<17>,
    z_omega_val: Field<17>,
}

impl fmt::Display for Proof {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "a(s) = {}", self.a)?;
        writeln!(f, "b(s) = {}", self.b)?;
        writeln!(f, "c(s) = {}", self.c)?;
        writeln!(f, "zet(s) = {}", self.zet)?;
        writeln!(f, "t_lo(s) = {}", self.t_lo)?;
        writeln!(f, "t_mid(s) = {}", self.t_mid)?;
        writeln!(f, "t_hi(s) = {}", self.t_hi)?;
        writeln!(f, "w_zeta(s) = {}", self.w_zeta)?;
        writeln!(f, "w_zeta_omega(s) = {}", self.w_zeta_omega)?;
        writeln!(f, "a(zeta) = {}", self.a_val)?;
        writeln!(f, "b(zeta) = {}", self.b_val)?;
        writeln!(f, "c(zeta) = {}", self.c_val)?;
        writeln!(f, "sigma1(zeta) = {}", self.s1_val)?;
        writeln!(f, "sigma2(zeta) = {}", self.s2_val)?;
        writeln!(f, "r(zeta) = {}", self.r_val)?;
        writeln!(f, "zet(zeta * omega) = {}", self.z_omega_val)?;

        Ok(())
    }
}

impl Circuit<4> {
    pub fn pythagorean() -> Self {
        // Pythagorean triple
        // 4 equations:
        // 3 * 3 - 9 = 0
        // 4 * 4 - 16 = 0
        // 5 * 5 - 25 = 0
        // 9 + 16 - 25 = 0
        Circuit {
            assignments: Assignment {
                a: [3, 4, 5, 9],
                b: [3, 4, 5, 16],
                c: [9, 16, 25, 25],
            },
            selectors: Selector {
                l: Poly::interpolate_vector([0, 0, 0, 1]),
                r: Poly::interpolate_vector([0, 0, 0, 1]),
                zero: Poly::interpolate_vector([16, 16, 16, 16]), // 16 aka -1
                m: Poly::interpolate_vector([1, 1, 1, 0]),
                c: Poly::interpolate_vector([0, 0, 0, 0]),
            },
            // TODO:
            // I did not really understand where those numbers are come from
            copy_constraints: CopyConstraints {
                s1: Poly::interpolate_vector([2, 8, 15, 3]),
                s2: Poly::interpolate_vector([1, 4, 16, 12]),
                s3: Poly::interpolate_vector([13, 9, 5, 14]),
            },
        }
    }

    // `beta` and `gamma` is challenge
    #[rustfmt::skip]
    pub fn prove(
        &self,
        alpha: Field<17>,
        beta: Field<17>,
        gamma: Field<17>,
        zeta: Field<17>,
        nu: Field<17>,
    ) -> Proof {
        const SIZE: usize = 4;

        // let r = rand::random::<[Field<17>; 9]>();
        let random = [7, 4, 11, 12, 16, 2, 14, 11, 7];

        let k1 = Field::<17>::new(2);
        let k2 = Field::<17>::new(3);

        // h = x^4 - 1
        let h = Poly::new([16, 0, 0, 0, 1]).resize::<{ SIZE + 2 }>();

        // round 1
        // hide assignments using random `r`
        let a = (h * Poly::new([random[1], random[0]]).resize::<{ SIZE + 2 }>())
            + Poly::interpolate_vector(self.assignments.a).resize();
        // println!("a = {a}");
        let b = (h * Poly::new([random[3], random[2]]).resize::<{ SIZE + 2 }>())
            + Poly::interpolate_vector(self.assignments.b).resize();
        // println!("b = {b}");
        let c = (h * Poly::new([random[5], random[4]]).resize::<{ SIZE + 2 }>())
            + Poly::interpolate_vector(self.assignments.c).resize();
        // println!("c = {c}");

        // round 2
        let asg = &self.assignments;
        let cc = &self.copy_constraints;
        let factor = |i: usize| {
            Field::<17>::new(1) * (Field::<17>::from(asg.a[i]) + matrix::OMEGA[i] * beta + gamma)
                / (Field::<17>::from(asg.a[i]) + cc.s1.evaluate(matrix::OMEGA[i]) * beta + gamma)
                * (Field::<17>::from(asg.b[i]) + k1 * matrix::OMEGA[i] * beta + gamma)
                / (Field::<17>::from(asg.b[i]) + cc.s2.evaluate(matrix::OMEGA[i]) * beta + gamma)
                * (Field::<17>::from(asg.c[i]) + k2 * matrix::OMEGA[i] * beta + gamma)
                / (Field::<17>::from(asg.c[i]) + cc.s3.evaluate(matrix::OMEGA[i]) * beta + gamma)
        };

        // SIZE = 4
        let acc0 = Field::<17>::new(1);
        let acc1 = acc0 * factor(0);
        let acc2 = acc1 * factor(1);
        let acc3 = acc2 * factor(2);
        let acc = Poly::interpolate_vector([acc0, acc1, acc2, acc3]);
        // println!("acc = {acc}");

        let h = h.resize::<7>();
        let zet = (Poly::new([random[8], random[7], random[6]]).resize::<{ SIZE + 3 }>() * h)
            + acc.resize::<7>();
        // println!("zet = {zet}");

        // Lagrange basis polynomial
        let l1 = Poly::interpolate_vector([1, 0, 0, 0]);

        const D: usize = 3 * SIZE + 10;
        let se = &self.selectors;
        let t_zh = a.resize::<D>() * b.resize::<D>() * se.m.resize::<D>()
            + a.resize::<D>() * se.l.resize::<D>() + b.resize::<D>() * se.r.resize::<D>()
            + c.resize::<D>() * se.zero.resize::<D>() + se.c.resize::<D>()
            +
            (a.resize::<D>() + Poly::new([gamma, beta]).resize::<D>()) *
            (b.resize::<D>() + Poly::new([gamma, beta * k1]).resize::<D>()) *
            (c.resize::<D>() + Poly::new([gamma, beta * k2]).resize::<D>()) * zet.resize::<D>() * alpha
            -
            (a.resize::<D>() + cc.s1.resize::<D>() * beta + Poly::new([gamma]).resize::<D>()) *
            (b.resize::<D>() + cc.s2.resize::<D>() * beta + Poly::new([gamma]).resize::<D>()) *
            (c.resize::<D>() + cc.s3.resize::<D>() * beta + Poly::new([gamma]).resize::<D>()) * zet.substitute(matrix::OMEGA[1]).resize::<D>() * alpha  // zet(omega * x)
            +
            (zet.resize::<D>() - Poly::new([1]).resize::<D>()) * l1.resize::<D>() * (alpha * alpha);
        // println!("t_zh = {t_zh}");
        let t = t_zh.divide_by_h4();
        // println!("t = {t}");
        let (t_lo, t_mid, t_hi) = t.split::<{ SIZE + 2 }>();
        // println!("t_lo = {t_lo}");
        // println!("t_mid = {t_mid}");
        // println!("t_hi = {t_hi}");

        // round 4
        let a_val = a.evaluate(zeta);
        let b_val = b.evaluate(zeta);
        let c_val = c.evaluate(zeta);
        let s1_val = cc.s1.evaluate(zeta);
        let s2_val = cc.s2.evaluate(zeta);
        let t_val = t.evaluate(zeta);
        // println!("a(zeta) = {a_val}, b(zeta) = {b_val}, c(zeta) = {c_val}");
        // println!("s1(zeta) = {s1_val}, s2(zeta) = {s2_val}, t(zeta) = {t_val}");

        let z_omega_val = zet.evaluate(zeta * matrix::OMEGA[1]);

        let r = {
            let _1 = (se.m * (a_val * b_val) + se.l * a_val + se.r * b_val + se.zero * c_val).resize();
            let _2 = zet * ((a_val + beta * zeta + gamma) * (b_val + k1 * beta * zeta + gamma) * (c_val + k2 * beta * zeta + gamma) * alpha);
            let _3 = cc.s3.resize() * ((a_val + beta * s1_val + gamma) * (b_val + beta * s2_val + gamma) * beta * z_omega_val * alpha);
            let _4 = zet * (l1.evaluate(zeta) * alpha * alpha);
            _1 + _2 - _3 + _4
        };
        // println!("r = {r}");
        let r_val = r.evaluate(zeta);
        // println!("r(zeta) = {r_val}");

        // round 5
        // w_zeta * (x - zeta)
        let w_zeta_prod = {
            let q = zeta.power(SIZE + 2);
            (t_lo + t_mid * q + t_hi * (q * q) - t_val.into()).resize()
                + (r - r_val.into()) * nu
                + (a - a_val.into()).resize() * nu.power(2)
                + (b - b_val.into()).resize() * nu.power(3)
                + (c - c_val.into()).resize() * nu.power(4)
                + (cc.s1 - s1_val.into()).resize() * nu.power(5)
                + (cc.s2 - s2_val.into()).resize() * nu.power(6)
        };
        // println!("w_zeta * (x - zeta) = {w_zeta_prod}");
        let w_zeta = w_zeta_prod.divide_by_first(zeta);
        // println!("w_zeta = {w_zeta}");

        let w_zeta_omega = (zet - z_omega_val.into()).divide_by_first(zeta * matrix::OMEGA[1]);
        // println!("w_zeta_omega = {w_zeta_omega}");

        Proof {
            a: a.evaluate_srs(),
            b: b.evaluate_srs(),
            c: c.evaluate_srs(),
            zet: zet.evaluate_srs(),
            t_lo: t_lo.evaluate_srs(),
            t_mid: t_mid.evaluate_srs(),
            t_hi: t_hi.evaluate_srs(),
            w_zeta: w_zeta.evaluate_srs(),
            w_zeta_omega: w_zeta_omega.evaluate_srs(),
            a_val,
            b_val,
            c_val,
            s1_val,
            s2_val,
            r_val,
            z_omega_val,
        }
    }

    #[rustfmt::skip]
    pub fn verify(
        &self,
        proof: &Proof,
        alpha: Field<17>,
        beta: Field<17>,
        gamma: Field<17>,
        zeta: Field<17>,
        nu: Field<17>,
        u: Field<17>,
    ) -> bool {
        const SIZE: usize = 4;

        let k1 = Field::<17>::new(2);
        let k2 = Field::<17>::new(3);

        // step 4, compute zero polynomial evaluation `zeta ^ n - 1`
        let zet_eval = zeta.power(SIZE) - Field::new(1);
        println!("zeta ^ n - 1 = {zet_eval}");

        // step 5, compute Lagrange polynomial
        let l1_zeta = zet_eval / (Field::new(SIZE as u8) * (zeta - Field::new(1)));
        println!("L1(zeta) = {l1_zeta}");

        // step 6 compute public input polynomial evaluation
        let pi_zeta = Field::<17>::new(0);

        // step 7, compute quotient polynomial evaluation
        let t_val = {
            let _a = proof.a_val + beta * proof.s1_val + gamma;
            let _b = proof.b_val + beta * proof.s2_val + gamma;
            let _c = proof.c_val + gamma;
            let middle = (_a * _b * _c * proof.z_omega_val) * alpha;
            (proof.r_val + pi_zeta - middle - l1_zeta * alpha * alpha) / zet_eval
        };
        println!("t(zeta) = {t_val}");

        // step 8, compute first part of batched polynomial commitment
        let se = &self.selectors;
        let cc = &self.copy_constraints;
        let coefficient_at_zet = {
            (proof.a_val + beta * zeta + gamma)
                * (proof.b_val + beta * k1 * zeta + gamma)
                * (proof.c_val + beta * k2 * zeta + gamma)
                * alpha * nu
            + l1_zeta * alpha * alpha * nu + u
        };
        let coefficient_at_s3 = {
            (proof.a_val + beta * proof.s1_val + gamma) * (proof.b_val + beta * proof.s2_val + gamma) * alpha * nu * beta * proof.z_omega_val
        };
        let d =
            // first line
            se.m.evaluate_srs().scalar_mul((proof.a_val * proof.b_val * nu).cast())
            + se.l.evaluate_srs().scalar_mul((proof.a_val * nu).cast())
            + se.r.evaluate_srs().scalar_mul((proof.b_val * nu).cast())
            + se.zero.evaluate_srs().scalar_mul((proof.c_val * nu).cast())
            + se.c.evaluate_srs().scalar_mul(nu.cast())
            // second line
            + proof.zet.scalar_mul(coefficient_at_zet.cast())
            // third line
            + cc.s3.evaluate_srs().scalar_mul((-coefficient_at_s3).cast());
        println!("d = {d}");

        // step 9, compute full batched polynomial commitment
        let t_zeta = {
            let q = zeta.power(SIZE + 2);
            proof.t_lo + proof.t_mid.scalar_mul(q.cast()) + proof.t_hi.scalar_mul((q * q).cast())
        };
        let f = t_zeta
            + d
            + proof.a.scalar_mul(nu.power(2).cast())
            + proof.b.scalar_mul(nu.power(3).cast())
            + proof.c.scalar_mul(nu.power(4).cast())
            + cc.s1.evaluate_srs().scalar_mul(nu.power(5).cast())
            + cc.s2.evaluate_srs().scalar_mul(nu.power(6).cast());
        println!("f = {f}");

        // step 10, compute group encoded batch evaluation
        let e_val = t_val
            + proof.r_val * nu
            + proof.a_val * nu.power(2)
            + proof.b_val * nu.power(3)
            + proof.c_val * nu.power(4)
            + proof.s1_val * nu.power(5)
            + proof.s2_val * nu.power(6)
            + u * proof.z_omega_val;
        let e = Point::G.scalar_mul(e_val.cast());
        println!("e = {e}");

        // step 11, batch validate all evaluations
        let lhs = proof.w_zeta + proof.w_zeta_omega.scalar_mul(u.cast());
        println!("lhs = {lhs}");
        let rhs = proof.w_zeta.scalar_mul(zeta.cast())
            + proof.w_zeta_omega.scalar_mul((u * zeta * matrix::OMEGA[1]).cast())
            + f + e.scalar_mul(Field::new(16));
        println!("rhs = {rhs}");

        // pairing equation: e(lhs, 2 * G_2) = e(rhs, G_2) where G_2 = { 36, 31u }, 2 * G_2 = { 90, 82u }
        // cheating, use curve logarithm table

        let it = || Point::CHEAT_TABLE.iter().enumerate();
        let lhs_scalar = it().find(|(_, &x)| x.eq(&lhs)).unwrap().0;
        let rhs_scalar = it().find(|(_, &x)| x.eq(&rhs)).unwrap().0;
        // need to verify e(G, 2 * G_2) ^ lhs_scalar = e(G, G_2) ^ rhs_scalar
        // e({ 1, 2 }, { 90, 82u }) ^ lhs_scalar = e({ 1, 2 }, { 36, 31u }) ^ rhs_scalar

        // TODO: perform this computation by hands
        // precomputed f_17(90, 82u)
        let f_17_90_82u = ExtendedField101(68.into(), 47.into());
        // precomputed f_17(36, 31u)
        let f_17_36_31u = ExtendedField101(8.into(), 61.into());

        // 101 is order of prime field, 17 is order of elliptic curve subgroup
        // 101 ^ 2 - 1 / 17 = 600
        // need to check (f_17_90_82u ^ 600) ^ lhs_scalar == (f_17_36_31u ^ 600) ^ rhs_scalar
        // u ^ 2 = -2 (mod 101)

        // naive exponent, not optimal
        let mut l = ExtendedField101::E;
        for _ in 0..(lhs_scalar * 600) {
            l = l * f_17_90_82u;
        }
        let mut r = ExtendedField101::E;
        for _ in 0..(rhs_scalar * 600) {
            r = r * f_17_36_31u;
        }

        println!("l = {l}, r = {r}");

        l == r
    }
}

impl<const N: usize> fmt::Display for Circuit<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "f_a = {}", Poly::interpolate_vector(self.assignments.a))?;
        writeln!(f, "f_b = {}", Poly::interpolate_vector(self.assignments.b))?;
        writeln!(f, "f_c = {}", Poly::interpolate_vector(self.assignments.c))?;
        writeln!(f, "q_l = {}", self.selectors.l)?;
        writeln!(f, "q_r = {}", self.selectors.r)?;
        writeln!(f, "q_0 = {}", self.selectors.zero)?;
        writeln!(f, "q_m = {}", self.selectors.m)?;
        writeln!(f, "q_c = {}", self.selectors.c)?;
        writeln!(f, "s_1 = {}", self.copy_constraints.s1)?;
        writeln!(f, "s_2 = {}", self.copy_constraints.s2)?;
        writeln!(f, "s_3 = {}", self.copy_constraints.s3)?;

        Ok(())
    }
}

#[cfg(test)]
#[test]
fn circuit() {
    let c = Circuit::pythagorean();
    println!("circuit:");
    println!("{c}");

    let alpha = Field::new(15);
    let beta = Field::new(12);
    let gamma = Field::new(13);
    let zeta = Field::new(5);
    let nu = Field::new(12);
    let proof = c.prove(alpha, beta, gamma, zeta, nu);
    println!("proof:");
    println!("{proof}");

    let u = Field::new(4);
    assert!(c.verify(&proof, alpha, beta, gamma, zeta, nu, u));
}
