use serde::{Deserialize, Serialize};

/// Deterministic, portable `XorShift64*` PRNG so saved `.popworkset.json` states
/// ("save games") reproduce identical stochastic trajectories across platforms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimRng {
    pub state: u64,
}

impl Default for SimRng {
    fn default() -> Self {
        Self::new(0x9E37_79B9_7F4A_7C15)
    }
}

impl SimRng {
    pub fn new(seed: u64) -> Self {
        let s = if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        };
        Self { state: s }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform float in `[0.0, 1.0)`.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        let bits = self.next_u64() >> 11;
        (bits as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    /// Uniform float in `[low, high)`.
    #[inline]
    pub fn range_f64(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.next_f64()
    }

    /// Uniform integer in `0..upper` (`upper > 0`).
    #[inline]
    pub fn gen_index(&mut self, upper: usize) -> usize {
        if upper <= 1 {
            0
        } else {
            (self.next_u64() % (upper as u64)) as usize
        }
    }

    /// Standard normal deviate via Box-Muller transform.
    pub fn standard_normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Approximate Binomial(n, p) draw for Wright-Fisher genetic drift and survival steps.
    pub fn binomial(&mut self, n: usize, p: f64) -> usize {
        let p_clamped = p.clamp(0.0, 1.0);
        if n == 0 || p_clamped <= 0.0 {
            return 0;
        }
        if p_clamped >= 1.0 {
            return n;
        }
        if n <= 64 {
            let mut successes = 0usize;
            for _ in 0..n {
                if self.next_f64() < p_clamped {
                    successes += 1;
                }
            }
            successes
        } else {
            let mean = (n as f64) * p_clamped;
            let std_dev = (mean * (1.0 - p_clamped)).sqrt();
            let sample = (mean + std_dev * self.standard_normal()).round();
            sample.clamp(0.0, n as f64) as usize
        }
    }

    /// Exponential waiting time with rate `lambda > 0`.
    pub fn exponential(&mut self, lambda: f64) -> f64 {
        if lambda <= 0.0 {
            return 0.0;
        }
        let u = self.next_f64().max(1e-12);
        -u.ln() / lambda
    }
}
