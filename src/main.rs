use std::ops::{Add, Div, Mul};
use wasm_bindgen::prelude::*;

const PI: f64 = std::f64::consts::PI;
const G: f64 = 9.81;
const L1: f64 = 1.0;
const L2: f64 = 1.0;
const M1: f64 = 1.0;
const M2: f64 = 1.0;
const STEP: f64 = 0.01;

#[derive(Clone, Copy)]
struct RK4 {
    alpha: f64,
    beta: f64,
    alpha_dot: f64,
    beta_dot: f64,
}

impl RK4 {
    fn new(alpha: f64, beta: f64, alpha_dot: f64, beta_dot: f64) -> Self {
        Self {
            alpha,
            beta,
            alpha_dot,
            beta_dot,
        }
    }

    fn alpha_ddot(&self) -> f64 {
        let a = M2 * G * L1 * self.beta.sin() * (self.alpha - self.beta).cos();
        let b = - M2 * L1.powi(2) * self.alpha_dot.powi(2) * (self.alpha - self.beta).sin() * (self.alpha - self.beta).cos();
        let c = - (M1 + M2) * G * L1 * self.alpha.sin();
        let d = - M2 * L1 * L2 * self.beta_dot.powi(2) * (self.alpha - self.beta).sin();
        let e = (M1 + M2) * L1.powi(2);
        let f = - M2 * L1.powi(2) * (self.alpha - self.beta).cos().powi(2);
        return (a + b + c + d) / (e + f)
    }


    fn beta_ddot(&self) -> f64 {
        let a = L1 * self.alpha_dot.powi(2) * (self.alpha - self.beta).sin();
        let b = - L1 * self.alpha_ddot() * (self.alpha - self.beta).cos();
        let c = - G * self.alpha.sin();
        return (a + b + c) / L2
    }


    fn xi_dot(&self) -> Self {
        return Self {
            alpha: self.alpha_dot,
            beta: self.beta_dot,
            alpha_dot: self.alpha_ddot(),
            beta_dot: self.beta_ddot(),
        }
    }


    fn step(&mut self) {
        let k1 = self.xi_dot();
        let k2 = (*self + k1 * STEP / 2.0).xi_dot();
        let k3 = (*self + k2 * STEP / 2.0).xi_dot();
        let k4 = (*self + k3 * STEP).xi_dot();
        let k = (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) * STEP / 6.0;
        
        self.alpha = self.alpha + k.alpha;
        self.beta = self.beta + k.beta;
        self.alpha_dot = self.alpha_dot + k.alpha_dot;
        self.beta_dot = self.beta_dot + k.beta_dot;
    }
}

impl Add<RK4> for RK4 {
    type Output = RK4;
    fn add(self, rhs: RK4) -> Self::Output {
        RK4 {
            alpha: self.alpha + rhs.alpha,
            beta: self.beta + rhs.beta,
            alpha_dot: self.alpha_dot + rhs.alpha_dot,
            beta_dot: self.beta_dot + rhs.beta_dot,
        }
    }
}

impl<T> Mul<T> for RK4
where
    f64: From<T>,
    T: Copy, 
{
    type Output = RK4;
    fn mul(self, rhs: T) -> Self::Output {
        RK4 {
            alpha: self.alpha * f64::from(rhs),
            beta: self.beta * f64::from(rhs),
            alpha_dot: self.alpha_dot * f64::from(rhs),
            beta_dot: self.beta_dot * f64::from(rhs),
        }
    }
}

impl<T> Div<T> for RK4
where
    f64: From<T>,
    T: Copy, 
{
    type Output = RK4;
    fn div(self, rhs: T) -> Self::Output {
        RK4 {
            alpha: self.alpha / f64::from(rhs),
            beta: self.beta / f64::from(rhs),
            alpha_dot: self.alpha_dot / f64::from(rhs),
            beta_dot: self.beta_dot / f64::from(rhs),
        }
    }
}

#[wasm_bindgen]
pub fn run(alpha: f64, beta: f64, alpha_dot: f64, beta_dot: f64) {
    let mut rk4 = RK4::new(alpha, beta, alpha_dot, beta_dot);
    for i in 0..10000 {
        println!("{} {} {} {} {}", i as f64 * STEP, rk4.alpha, rk4.beta, rk4.alpha_dot, rk4.beta_dot);
        rk4.step();
    }
}


fn main() {
    run(PI / 2.0, PI / 2.0, 0.0, 0.0);
}