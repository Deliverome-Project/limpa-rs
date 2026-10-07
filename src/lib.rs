//! LIMPA's Gaussian observed-data likelihood, logistic-normal missingness likelihood,
//! and protein priors. Peptide effects sum to zero. The 16-point quadrature matches
//! statmod::gauss.quad.prob; no data imputation enters the likelihood.
use nalgebra::{DMatrix, DVector};

// Preserve the exact decimal exports from statmod for auditability.
#[allow(clippy::excessive_precision)]
const Z: [f64; 16] = [
    -6.6308781983931286,
    -5.472225705949346,
    -4.4929553025200137,
    -3.600873624171546,
    -2.7602450476307001,
    -1.9519803457163334,
    -1.1638291005549648,
    -0.38676060450055721,
    0.38676060450055744,
    1.1638291005549644,
    1.9519803457163321,
    2.7602450476307006,
    3.60087362417155,
    4.4929553025200129,
    5.472225705949346,
    6.6308781983931322,
];
#[allow(clippy::excessive_precision)]
const W: [f64; 16] = [
    1.4978147231618505e-10,
    1.3094732162868192e-07,
    1.5300032162487211e-05,
    0.00052598492657390751,
    0.0072669376011847168,
    0.047284752354014234,
    0.15833837275094964,
    0.28656852123801252,
    0.28656852123801241,
    0.15833837275094959,
    0.047284752354014123,
    0.0072669376011847333,
    0.00052598492657390588,
    1.5300032162487232e-05,
    1.3094732162868129e-07,
    1.4978147231618433e-10,
];

#[derive(Clone, Copy, Debug)]
pub struct Model {
    pub intercept: f64,
    pub slope: f64,
    pub prior_mean: f64,
    pub prior_sd: f64,
    pub prior_logfc: f64,
}
#[derive(Debug)]
pub struct Protein {
    pub samples: usize,
    pub peptides: usize,
    pub sigma: f64,
    /// Row-major peptide x sample values, NaN denotes missing.
    pub y: Vec<f64>,
    pub start: Vec<f64>,
}
#[derive(Debug)]
pub struct Fit {
    pub expression: Vec<f64>,
    pub se: Vec<f64>,
    pub iterations: usize,
    pub gradient_max: f64,
}
struct Evaluation {
    value: f64,
    gradient: DVector<f64>,
    weights: Vec<f64>,
}

/// Stable log probability and first two derivatives of -log P(missing).
fn missing(mu: f64, sigma: f64, m: Model) -> (f64, f64, f64) {
    // Ordinary intensity range: direct integration avoids 16 logs per cell.
    // Keep log-sum-exp below as a fallback for extreme detection probabilities.
    let mut mass = 0.0;
    let mut first = 0.0;
    let mut second = 0.0;
    for k in 0..16 {
        let eta = m.intercept + m.slope * (mu + sigma * Z[k]);
        let q = 1.0 / (1.0 + eta.exp());
        let d = 1.0 - q;
        let a = W[k] * q;
        mass += a;
        first += a * d;
        second += a * d * (2.0 * d - 1.0);
    }
    if mass > 1e-280 {
        let r = first / mass;
        return (
            -mass.ln(),
            m.slope * first / mass,
            m.slope * m.slope * (r * r - second / mass),
        );
    }
    let mut logs = [0.0_f64; 16];
    let mut detected = [0.0_f64; 16];
    let mut top = f64::NEG_INFINITY;
    for k in 0..16 {
        let eta = m.intercept + m.slope * (mu + sigma * Z[k]);
        let softplus = eta.max(0.0) + (-eta.abs()).exp().ln_1p();
        logs[k] = W[k].ln() - softplus;
        detected[k] = if eta >= 0.0 {
            1.0 / (1.0 + (-eta).exp())
        } else {
            let e = eta.exp();
            e / (1.0 + e)
        };
        top = top.max(logs[k]);
    }
    let mut mass = 0.0;
    let mut first = 0.0;
    let mut second = 0.0;
    for k in 0..16 {
        let a = (logs[k] - top).exp();
        let d = detected[k];
        mass += a;
        first += a * d;
        second += a * (2.0 * d * d - d);
    }
    let r = first / mass;
    (
        -(top + mass.ln()),
        m.slope * first / mass,
        m.slope * m.slope * (r * r - second / mass),
    )
}
fn evaluate(p: &Protein, m: Model, beta: &DVector<f64>) -> Evaluation {
    let n = p.samples;
    let q = p.peptides - 1;
    let mut gradient = DVector::zeros(n + q);
    let mut weights = vec![0.0; p.y.len()];
    let mut observed = 0.0;
    let mut missing_log = 0.0;
    // R flattens the peptide matrix column-major. Preserve its reduction order
    // and arithmetic grouping: tiny roundoff can alter BFGS's stopping iteration.
    for j in 0..n {
        for i in 0..p.peptides {
            let k = i * n + j;
            let mut mu = beta[j];
            if i < q {
                mu += beta[n + i];
            } else {
                for a in 0..q {
                    mu -= beta[n + a];
                }
            }
            let (g, h) = if p.y[k].is_nan() {
                let (v, g, h) = missing(mu, p.sigma, m);
                missing_log -= v;
                (g, h)
            } else {
                let r = mu - p.y[k];
                observed += (r / p.sigma).powi(2);
                (r / p.sigma.powi(2), 1.0 / p.sigma.powi(2))
            };
            gradient[j] += g;
            weights[k] = h;
            if i < q {
                gradient[n + i] += g;
            } else {
                for a in 0..q {
                    gradient[n + a] -= g;
                }
            }
        }
    }
    let mut avg = (0..n).map(|j| beta[j]).sum::<f64>() / n as f64;
    avg += (0..n).map(|j| beta[j] - avg).sum::<f64>() / n as f64;
    let mut centered_sum = 0.0;
    for j in 0..n {
        let centered = beta[j] - avg;
        centered_sum += centered * centered;
        gradient[j] +=
            (avg - m.prior_mean) / n as f64 / m.prior_sd.powi(2) + centered / m.prior_logfc.powi(2);
    }
    let prior =
        (avg - m.prior_mean).powi(2) / m.prior_sd.powi(2) + centered_sum / m.prior_logfc.powi(2);
    let value = 0.5 * (observed - 2.0 * missing_log + prior);
    Evaluation {
        value,
        gradient,
        weights,
    }
}

/// H = [A B; B' C], A diagonal + rank one. Eliminate the samples;
/// only the (peptides-1)^2 Schur complement is factorized.
struct Factor {
    invd: Vec<f64>,
    rank: f64,
    u: DMatrix<f64>,
    sinv: DMatrix<f64>,
}
impl Factor {
    fn new(p: &Protein, m: Model, w: &[f64]) -> Result<Self, String> {
        let n = p.samples;
        let q = p.peptides - 1;
        let lambda = 1.0 / m.prior_logfc.powi(2);
        let c = 1.0 / (n as f64).powi(2) / m.prior_sd.powi(2) - lambda / n as f64;
        let invd: Vec<f64> = (0..n)
            .map(|j| 1.0 / (lambda + (0..p.peptides).map(|i| w[i * n + j]).sum::<f64>()))
            .collect();
        let denom = 1.0 + c * invd.iter().sum::<f64>();
        if denom <= 0.0 || !denom.is_finite() {
            return Err("sample Hessian is not positive definite".into());
        }
        let rank = c / denom;
        let b = DMatrix::from_fn(n, q, |j, k| w[k * n + j] - w[q * n + j]);
        let sums: Vec<f64> = (0..q)
            .map(|k| (0..n).map(|l| invd[l] * b[(l, k)]).sum())
            .collect();
        let u = DMatrix::from_fn(n, q, |j, k| invd[j] * (b[(j, k)] - rank * sums[k]));
        let lastsum = (0..n).map(|j| w[q * n + j]).sum::<f64>();
        let mut s = DMatrix::from_fn(q, q, |a, b| {
            lastsum
                + if a == b {
                    (0..n).map(|j| w[a * n + j]).sum::<f64>()
                } else {
                    0.0
                }
        });
        s -= b.transpose() * &u;
        let sinv = if q == 0 {
            DMatrix::zeros(0, 0)
        } else {
            s.cholesky()
                .ok_or("peptide Hessian is not positive definite; unidentified peptide effects")?
                .inverse()
        };
        Ok(Self {
            invd,
            rank,
            u,
            sinv,
        })
    }
    fn solve(&self, g: &DVector<f64>) -> DVector<f64> {
        let n = self.invd.len();
        let q = self.u.ncols();
        let sum = (0..n).map(|j| self.invd[j] * g[j]).sum::<f64>();
        let ag = DVector::from_fn(n, |j, _| self.invd[j] * (g[j] - self.rank * sum));
        let pg = &self.sinv * (g.rows(n, q) - self.u.transpose() * g.rows(0, n));
        let sg = ag - &self.u * &pg;
        DVector::from_iterator(n + q, sg.iter().chain(pg.iter()).copied())
    }
    fn se(&self) -> Vec<f64> {
        let us = &self.u * &self.sinv;
        (0..self.invd.len())
            .map(|j| {
                (self.invd[j] - self.rank * self.invd[j].powi(2)
                    + (0..self.u.ncols())
                        .map(|k| us[(j, k)] * self.u[(j, k)])
                        .sum::<f64>())
                .sqrt()
            })
            .collect()
    }
}
pub fn fit_newton(p: &Protein, m: Model) -> Result<Fit, String> {
    if p.samples < 2
        || p.peptides == 0
        || p.y.len() != p.samples * p.peptides
        || p.start.len() != p.samples + p.peptides - 1
    {
        return Err("invalid protein dimensions".into());
    }
    if ![p.sigma, m.prior_sd, m.prior_logfc]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
        || ![m.intercept, m.slope, m.prior_mean]
            .iter()
            .all(|x| x.is_finite())
        || p.y.iter().any(|x| x.is_infinite())
        || p.start.iter().any(|x| !x.is_finite())
    {
        return Err("invalid numerical input".into());
    }
    // Unobserved precursor effects are weakly identified; the R bridge uses
    // the reference fitter for this input rather than claiming equivalence.
    if p.y
        .chunks(p.samples)
        .filter(|row| row.iter().all(|x| x.is_nan()))
        .count()
        > 0
    {
        return Err("all-missing precursors require the reference fitter".into());
    }
    let mut beta = DVector::from_vec(p.start.clone());
    for iteration in 0..100 {
        let e = evaluate(p, m, &beta);
        let gradient_max = e.gradient.amax();
        let factor = Factor::new(p, m, &e.weights)?;
        if gradient_max < 1e-7 {
            let se = factor.se();
            if se.iter().any(|x| !x.is_finite() || *x <= 0.0) {
                return Err("invalid standard error".into());
            }
            return Ok(Fit {
                expression: beta.rows(0, p.samples).iter().copied().collect(),
                se,
                iterations: iteration,
                gradient_max,
            });
        }
        let step = factor.solve(&e.gradient);
        let descent = e.gradient.dot(&step);
        if !descent.is_finite() || descent <= 0.0 {
            return Err("non-descent Newton direction".into());
        }
        let mut scale = 1.0;
        let mut accepted = false;
        for _ in 0..40 {
            let next = &beta - scale * &step;
            let nextval = evaluate(p, m, &next).value;
            if nextval.is_finite()
                && nextval <= e.value - 1e-4 * scale * descent + 1e-12 * (1.0 + e.value.abs())
            {
                beta = next;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            return Err("Newton line search failed".into());
        }
    }
    Err("Newton solver did not converge in 100 iterations".into())
}

/// Compatibility optimizer: R's optim(method="BFGS") defaults and restart rules.
/// Adapted from R src/appl/optim.c, vmmin, R Core Team (GPL >= 2).
/// Retaining stopping rules matters for weakly identified missing-data fits.
pub fn fit(p: &Protein, m: Model) -> Result<Fit, String> {
    if p.samples < 2
        || p.peptides == 0
        || p.y.len() != p.samples * p.peptides
        || p.start.len() != p.samples + p.peptides - 1
        || ![p.sigma, m.prior_sd, m.prior_logfc]
            .iter()
            .all(|x| x.is_finite() && *x > 0.0)
        || ![m.intercept, m.slope, m.prior_mean]
            .iter()
            .all(|x| x.is_finite())
        || p.y.iter().any(|x| x.is_infinite())
        || p.start.iter().any(|x| !x.is_finite())
    {
        return Err("invalid numerical input or dimensions".into());
    }
    if p.y
        .chunks(p.samples)
        .filter(|r| r.iter().all(|x| x.is_nan()))
        .count()
        > 0
    {
        return Err("all-missing precursors require the reference fitter".into());
    }
    let n = p.start.len();
    let mut beta = DVector::from_vec(p.start.clone());
    let e = evaluate(p, m, &beta);
    let mut fmin = 2.0 * e.value;
    let mut f = fmin;
    let mut g = 2.0 * e.gradient;
    let mut b = vec![vec![0.0; n]; n];
    let mut t = vec![0.0; n];
    let mut x = vec![0.0; n];
    let mut c = vec![0.0; n];
    let mut gradcount = 1;
    let mut iter = 1;
    let mut ilast = 1;
    let reltol = f64::EPSILON.sqrt();
    loop {
        if ilast == gradcount {
            for (i, row) in b.iter_mut().enumerate() {
                row[..i].fill(0.0);
                row[i] = 1.0;
            }
        }
        x.copy_from_slice(beta.as_slice());
        c.copy_from_slice(g.as_slice());
        let mut proj = 0.0;
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..=i {
                s -= b[i][j] * g[j];
            }
            for j in i + 1..n {
                s -= b[j][i] * g[j];
            }
            t[i] = s;
            proj += s * g[i];
        }
        let mut count;
        if proj < 0.0 {
            let mut scale = 1.0;
            loop {
                count = 0;
                for i in 0..n {
                    beta[i] = x[i] + scale * t[i];
                    if 10.0 + x[i] == 10.0 + beta[i] {
                        count += 1;
                    }
                }
                if count == n {
                    break;
                }
                f = 2.0 * evaluate(p, m, &beta).value;
                if f.is_finite() && f <= fmin + proj * scale * 0.0001 {
                    break;
                }
                scale *= 0.2;
            }
            if (f - fmin).abs() <= reltol * (fmin.abs() + reltol) {
                count = n;
                fmin = f;
            }
            if count < n {
                fmin = f;
                g = 2.0 * evaluate(p, m, &beta).gradient;
                gradcount += 1;
                iter += 1;
                let mut d1 = 0.0;
                for i in 0..n {
                    t[i] *= scale;
                    c[i] = g[i] - c[i];
                    d1 += t[i] * c[i];
                }
                if d1 > 0.0 {
                    let mut d2 = 0.0;
                    for i in 0..n {
                        let mut s = 0.0;
                        for j in 0..=i {
                            s += b[i][j] * c[j];
                        }
                        for j in i + 1..n {
                            s += b[j][i] * c[j];
                        }
                        x[i] = s;
                        d2 += s * c[i];
                    }
                    d2 = 1.0 + d2 / d1;
                    for i in 0..n {
                        for j in 0..=i {
                            b[i][j] += (d2 * t[i] * t[j] - x[i] * t[j] - t[i] * x[j]) / d1;
                        }
                    }
                } else {
                    ilast = gradcount;
                }
            } else if ilast < gradcount {
                count = 0;
                ilast = gradcount;
            }
        } else {
            count = 0;
            if ilast == gradcount {
                count = n;
            } else {
                ilast = gradcount;
            }
        }
        if iter >= 100 {
            return Err(
                "R-compatible BFGS reached maxit=100; reference would return convergence=1".into(),
            );
        }
        if gradcount - ilast > 2 * n {
            ilast = gradcount;
        }
        if count == n && ilast == gradcount {
            break;
        }
    }
    let e = evaluate(p, m, &beta);
    let factor = Factor::new(p, m, &e.weights)?;
    let se = factor.se();
    if se.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err("invalid standard error".into());
    }
    Ok(Fit {
        expression: beta.rows(0, p.samples).iter().copied().collect(),
        se,
        iterations: iter,
        gradient_max: 2.0 * e.gradient.amax(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn model() -> Model {
        Model {
            intercept: -12.0,
            slope: 0.75,
            prior_mean: 18.0,
            prior_sd: 3.0,
            prior_logfc: 2.0,
        }
    }
    #[test]
    fn compatibility_fit_matches_limpa_1_4_2_oracle() {
        // Generated by limpa::peptides2Proteins on this synthetic matrix with
        // global imputeByExpTilt(..., dpc.slope=.75, prior.logfc=2) starts.
        let p = Protein {
            samples: 3,
            peptides: 2,
            sigma: 0.4,
            y: vec![18.0, 19.0, f64::NAN, 20.0, f64::NAN, 21.0],
            start: vec![19.0, 18.25, 19.7269230769231, -0.507692307692306],
        };
        let mut m = model();
        m.intercept = -11.0;
        let f = fit(&p, m).unwrap();
        let expected = [19.0067389354125, 19.8567893391462, 19.8685423644212];
        let se = [0.280856430943388, 0.475136210471305, 0.473790957266767];
        for j in 0..3 {
            assert!((f.expression[j] - expected[j]).abs() < 1e-6);
            assert!((f.se[j] - se[j]).abs() < 1e-7);
        }
    }
    #[test]
    fn invalid_and_all_missing_rows_fail_closed() {
        let mut p = Protein {
            samples: 3,
            peptides: 1,
            sigma: 0.4,
            y: vec![f64::NAN; 3],
            start: vec![18.0; 3],
        };
        assert!(fit(&p, model()).is_err());
        p.y = vec![18.0; 3];
        p.sigma = 0.0;
        assert!(fit(&p, model()).is_err());
    }
    #[test]
    fn derivatives_match_finite_differences() {
        for mu in [-1000.0, 10.0, 18.0, 30.0, 1000.0] {
            let (v, g, h) = missing(mu, 0.5, model());
            let eps = 1e-3;
            let (vp, gp, _) = missing(mu + eps, 0.5, model());
            let (vm, gm, _) = missing(mu - eps, 0.5, model());
            assert!(v.is_finite());
            assert!((g - (vp - vm) / (2.0 * eps)).abs() < 1e-7);
            assert!((h - (gp - gm) / (2.0 * eps)).abs() < 1e-7);
        }
    }
    #[test]
    fn complete_constant_data_has_analytic_solution() {
        let p = Protein {
            samples: 4,
            peptides: 3,
            sigma: 0.5,
            y: vec![20.0; 12],
            start: vec![0.0; 6],
        };
        let f = fit_newton(&p, model()).unwrap();
        let expected = (12.0 * 20.0 + 18.0 / 36.0) / (12.0 + 1.0 / 36.0);
        for x in f.expression {
            assert!((x - expected).abs() < 1e-8);
        }
        assert!(f.se.iter().all(|x| x.is_finite() && *x > 0.0));
    }
    #[test]
    fn structured_inverse_matches_dense_hessian() {
        let p = Protein {
            samples: 4,
            peptides: 3,
            sigma: 0.5,
            y: vec![
                20.0,
                19.0,
                f64::NAN,
                21.0,
                18.0,
                19.0,
                20.0,
                21.0,
                22.0,
                23.0,
                24.0,
                25.0,
            ],
            start: vec![18.0, 19.0, 20.0, 21.0, -1.0, 0.0],
        };
        let b = DVector::from_vec(p.start.clone());
        let e = evaluate(&p, model(), &b);
        let factor = Factor::new(&p, model(), &e.weights).unwrap();
        let mut h = DMatrix::zeros(6, 6);
        for k in 0..6 {
            let mut bp = b.clone();
            let mut bm = b.clone();
            bp[k] += 1e-4;
            bm[k] -= 1e-4;
            h.set_column(
                k,
                &((evaluate(&p, model(), &bp).gradient - evaluate(&p, model(), &bm).gradient)
                    / 2e-4),
            );
        }
        let inv = h.try_inverse().unwrap();
        assert!((factor.solve(&e.gradient) - &inv * &e.gradient).amax() < 1e-7);
        for (j, se) in factor.se().iter().enumerate() {
            assert!((se * se - inv[(j, j)]).abs() < 1e-7);
        }
    }
}
