use super::field::Field;

pub const OMEGA: [Field<17>; 4] = [Field::new(1), Field::new(4), Field::new(16), Field::new(13)];

// this is the matrix for polynomial interpolation
// it consists of roots of 1, and their powers
// it is very easy to find inverse matrix
pub fn h_matrix() -> Vec<[Field<17>; 4]> {
    let width = 4;
    let mut t = [Field::<17>::from(1), 1.into(), 1.into(), 1.into()];
    (0..width)
        .map(|_| {
            let r = t;
            for (t, root) in t.iter_mut().zip(OMEGA) {
                *t = *t * root;
            }
            r
        })
        .collect()
}

pub fn h_inv_matrix() -> Vec<[Field<17>; 4]> {
    vec![
        [13.into(), 13.into(), 13.into(), 13.into()],
        [13.into(), 16.into(), 4.into(), 1.into()],
        [13.into(), 4.into(), 13.into(), 4.into()],
        [13.into(), 1.into(), 4.into(), 16.into()],
    ]
}

#[cfg(test)]
#[test]
fn matrix() {
    let h = h_matrix();
    let h_inv = h_inv_matrix();
    for i in 0..4 {
        for j in 0..4 {
            let x = (0..4).fold(Field::<17>::new(0), |a, k| a + h[i][k] * h_inv[k][j]);
            if i == j {
                assert_eq!(x, Field::new(1));
            } else {
                assert_eq!(x, Field::new(0));
            }
        }
    }
}
