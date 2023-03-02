use std::ops::{Add, Div, Mul};
use wasm_bindgen::prelude::*;

const PI: f64 = std::f64::consts::PI;
const G: f64 = 9.81;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub struct DoublePendulum {
    alpha: f64,
    beta: f64,
    alpha_dot: f64,
    beta_dot: f64,
    lengths: (f64, f64),
    masses: (f64, f64),
    step: f64,
}

#[wasm_bindgen]
impl DoublePendulum {
    pub fn new(
        alpha: f64,
        beta: f64,
        alpha_dot: f64,
        beta_dot: f64,
        l1: f64,
        l2: f64,
        m1: f64,
        m2: f64,
        step: f64,
    ) -> Self {
        if alpha > PI || alpha < -PI {
            panic!("alpha must be between -PI and PI");
        }
        if beta > PI || beta < -PI {
            panic!("beta must be between -PI and PI");
        }
        if l1 <= 0.0 {
            panic!("l1 must be greater than 0");
        }
        if l2 <= 0.0 {
            panic!("l2 must be greater than 0");
        }
        if m1 <= 0.0 {
            panic!("m1 must be greater than 0");
        }
        if m2 <= 0.0 {
            panic!("m2 must be greater than 0");
        }
        if step <= 0.0 {
            panic!("step must be greater than 0");
        }
        Self {
            alpha,
            beta,
            alpha_dot,
            beta_dot,
            lengths: (l1, l2),
            masses: (m1, m2),
            step,
        }
    }

    fn alpha_ddot(&self) -> f64 {
        let (l1, l2) = self.lengths;
        let (m1, m2) = self.masses;
        let a = m2 * G * l1 * self.beta.sin() * (self.alpha - self.beta).cos();
        let b = - m2 * l1.powi(2) * self.alpha_dot.powi(2) * (self.alpha - self.beta).sin() * (self.alpha - self.beta).cos();
        let c = - (m1 + m2) * G * l1 * self.alpha.sin();
        let d = - m2 * l1 * l2 * self.beta_dot.powi(2) * (self.alpha - self.beta).sin();
        let e = (m1 + m2) * l1.powi(2);
        let f = - m2 * l1.powi(2) * (self.alpha - self.beta).cos().powi(2);
        return (a + b + c + d) / (e + f)
    }


    fn beta_ddot(&self) -> f64 {
        let (l1, l2) = self.lengths;
        let a = l1 * self.alpha_dot.powi(2) * (self.alpha - self.beta).sin();
        let b = - l1 * self.alpha_ddot() * (self.alpha - self.beta).cos();
        let c = - G * self.alpha.sin();
        return (a + b + c) / l2
    }


    fn xi_dot(&self) -> Self {
        return Self {
            alpha: self.alpha_dot,
            beta: self.beta_dot,
            alpha_dot: self.alpha_ddot(),
            beta_dot: self.beta_ddot(),
            lengths: self.lengths,
            masses: self.masses,
            step: self.step,
        }
    }

    pub fn step(&mut self) {
        let k1 = self.xi_dot();
        let k2 = (*self + k1 * self.step / 2.0).xi_dot();
        let k3 = (*self + k2 * self.step / 2.0).xi_dot();
        let k4 = (*self + k3 * self.step).xi_dot();
        let k = (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) * self.step / 6.0;
        
        self.alpha = self.alpha + k.alpha;
        self.beta = self.beta + k.beta;
        self.alpha_dot = self.alpha_dot + k.alpha_dot;
        self.beta_dot = self.beta_dot + k.beta_dot;
    }

    pub fn x1(&self) -> f64 {
        let (l1, _) = self.lengths;
        l1 * self.alpha.sin()
    }

    pub fn y1(&self) -> f64 {
        let (l1, _) = self.lengths;
        l1 * self.alpha.cos()
    }

    pub fn x2(&self) -> f64 {
        let (l1, l2) = self.lengths;
        l1 * self.alpha.sin() + l2 * self.beta.sin()
    }

    pub fn y2(&self) -> f64 {
        let (l1, l2) = self.lengths;
        l1 * self.alpha.cos() + l2 * self.beta.cos()
    }
}

// WASM BINDGEN HAS NOT IMPLEMENTED TRAITS YET

// pub trait RK4 {
//     fn xi_dot(&self) -> Self;
//     fn step(&mut self);
// }


// impl RK4 for DoublePendulum {
//     fn xi_dot(&self) -> Self {
//         return Self {
//             alpha: self.alpha_dot,
//             beta: self.beta_dot,
//             alpha_dot: self.alpha_ddot(),
//             beta_dot: self.beta_ddot(),
//             lengths: self.lengths,
//             masses: self.masses,
//             step: self.step,
//         }
//     }

//     fn step(&mut self) {
//         let k1 = self.xi_dot();
//         let k2 = (*self + k1 * self.step / 2.0).xi_dot();
//         let k3 = (*self + k2 * self.step / 2.0).xi_dot();
//         let k4 = (*self + k3 * self.step).xi_dot();
//         let k = (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) * self.step / 6.0;
        
//         self.alpha = self.alpha + k.alpha;
//         self.beta = self.beta + k.beta;
//         self.alpha_dot = self.alpha_dot + k.alpha_dot;
//         self.beta_dot = self.beta_dot + k.beta_dot;
//     }
// }

impl Add<DoublePendulum> for DoublePendulum {
    type Output = DoublePendulum;
    fn add(self, rhs: DoublePendulum) -> Self::Output {
        DoublePendulum {
            alpha: self.alpha + rhs.alpha,
            beta: self.beta + rhs.beta,
            alpha_dot: self.alpha_dot + rhs.alpha_dot,
            beta_dot: self.beta_dot + rhs.beta_dot,
            lengths: self.lengths,
            masses: self.masses,
            step: self.step,
        }
    }
}

impl<T> Mul<T> for DoublePendulum
where
    f64: From<T>,
    T: Copy, 
{
    type Output = DoublePendulum;
    fn mul(self, rhs: T) -> Self::Output {
        DoublePendulum {
            alpha: self.alpha * f64::from(rhs),
            beta: self.beta * f64::from(rhs),
            alpha_dot: self.alpha_dot * f64::from(rhs),
            beta_dot: self.beta_dot * f64::from(rhs),
            lengths: self.lengths,
            masses: self.masses,
            step: self.step,
        }
    }
}

impl<T> Div<T> for DoublePendulum
where
    f64: From<T>,
    T: Copy, 
{
    type Output = DoublePendulum;
    fn div(self, rhs: T) -> Self::Output {
        DoublePendulum {
            alpha: self.alpha / f64::from(rhs),
            beta: self.beta / f64::from(rhs),
            alpha_dot: self.alpha_dot / f64::from(rhs),
            beta_dot: self.beta_dot / f64::from(rhs),
            lengths: self.lengths,
            masses: self.masses,
            step: self.step,
        }
    }
}

pub fn run(
    alpha: f64,
    beta: f64,
    alpha_dot: f64,
    beta_dot: f64,
    l1: f64,
    l2: f64,
    m1: f64,
    m2: f64,
    step: f64,
) {
    let mut double_pendulum = DoublePendulum::new(alpha, beta, alpha_dot, beta_dot, l1, l2, m1, m2, step);
    let mut count = 0;
    loop {
        println!("{} {} {} {} {}", count as f64 * step, double_pendulum.alpha, double_pendulum.beta, double_pendulum.alpha_dot, double_pendulum.beta_dot);
        double_pendulum.step();
        count += 1;
    }
}


fn main() {
    run(PI / 2.0, PI / 2.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.001);
}