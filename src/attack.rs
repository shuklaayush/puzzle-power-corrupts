use ark_bls12_cheon::{Bls12Cheon, G1Projective as G1, G2Projective as G2, Fr};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::{PrimeField, Field, MontFp};

use num_integer::Roots; 
use std::collections::HashMap;

use crate::utils::{pow_sp, bigInt_to_u128};

fn baby_step_giant_step(eta: Fr, max: u128, gu0: &PairingOutput<Bls12Cheon>, gv0: &PairingOutput<Bls12Cheon>) -> u128 {
    // Create lookup table
    let eta_inv = eta.inverse().unwrap();
    let m = max.sqrt() + 1;
    let mut u: u128 = 0;
    let mut gu = gu0.clone();
    let mut map = HashMap::new();
    while u < m {
        // println!("{u}/{m}");
        map.insert(gu, u);
        gu *= eta_inv;
        u += 1;
    }

    // Check for collisions
    let eta_pow_m = pow_sp(eta, m, 80);
    let m_hat = max / m;
    let mut v: u128 = 0;
    let mut gv = gv0.clone();
    while v <= m_hat {
        // println!("{v}/{m_hat}");
        if map.contains_key(&gv) {
            u = map.get(&gv).copied().unwrap_or(0);
            break;
        }
        gv *= eta_pow_m;
        v += 1;
    }

    return u + m * v;

}

pub fn attack(P: G1, tau_P: G1, tau_d1_P: G1, Q: G2, tau_d2_Q: G2) -> i128 {
    /*

        IMPLEMENT YOUR ATTACK HERE

    */
    // Order of the subgroup
    let q = bigInt_to_u128(Fr::MODULUS);

    // Compute pairings to get three elliptic curve points of the form g, g^a, g^a^d
    let g = Bls12Cheon::pairing(P, Q);
    let g1 = Bls12Cheon::pairing(tau_P, Q);
    let gd = Bls12Cheon::pairing(tau_d1_P, tau_d2_Q);

    let two: Fr = MontFp!("2");
    let two_inv = two.inverse().unwrap();

    let d1: u128 = 11726539;
    let d2: u128 = 690320833;

    let d = d1 + d2;

    // Compute k0
    let eta = pow_sp(two, d, 80);
    let eta_inv = eta.inverse().unwrap();

    // Given A < k0 < B
    let A: u128 = 1089478584172543;
    let B: u128 = 1089547303649280;

    // Rewrite inequality as 0 < (k0 - A) < B - A
    // Solve for (k0 - A) using baby-step-giant-step
    let gu0 = gd * pow_sp(eta_inv, A, 80);
    let gv0 = g;

    let k0 = A + baby_step_giant_step(eta, B - A, &gu0, &gv0);

    println!("k0 : {k0}");

    // Compute k1
    let eta = pow_sp(two, (q - 1) / d, 80);

    let gu0 = g1 * pow_sp(two_inv, k0, 80);
    let gv0 = g;

    let k1 = baby_step_giant_step(eta, d, &gu0, &gv0);

    println!("k1 : {k1}");

    // Get k
    let k = k0 + k1 * (q - 1) / d;

    println!("k  : {k}");

    // Calculate tau
    let tau = pow_sp(two, k, 80);
    let tau = bigInt_to_u128(tau.into()) as i128;

    println!("tau: {tau}");

    return tau;
}
