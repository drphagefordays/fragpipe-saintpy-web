/**
 * saintpy-engine.js
 * Exact in-browser port of FossatiLab/saintPy (https://github.com/FossatiLab/saintPy)
 * + Universal FragPipe (combined_protein.tsv, protein.tsv, reprint.tsv, experiment_annotation.tsv)
 *   and SAINTexpress / saintPy parser, validator, and Bait-Prey Designer engine.
 */

(function (root, factory) {
  if (typeof module === "object" && module.exports) {
    module.exports = factory();
  } else {
    root.SaintPyEngine = factory();
  }
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  // ============================================================================
  // 1. MATHEMATICAL & STATISTICAL PRIMITIVES (Exact match to Python saintPy)
  // ============================================================================

  const LOG_FACT_MAX = 256;
  const LOG_FACT_TABLE = new Float64Array(LOG_FACT_MAX + 1);
  for (let i = 1; i <= LOG_FACT_MAX; i++) {
    LOG_FACT_TABLE[i] = LOG_FACT_TABLE[i - 1] + Math.log(i);
  }

  function lgamma(z) {
    if (z <= 0) return Infinity;
    if (Number.isInteger(z) && z - 1 <= LOG_FACT_MAX) {
      return LOG_FACT_TABLE[z - 1];
    }
    const g = 7;
    const p = [
      0.99999999999980993,
      676.5203681218851,
      -1259.1392167224028,
      771.32342877765313,
      -176.61502916214059,
      12.507343278686905,
      -0.13857109526572012,
      9.9843695780195716e-6,
      1.5056327351493116e-7,
    ];
    let x = z - 1;
    let a = p[0];
    const t = x + g + 0.5;
    for (let i = 1; i < p.length; i++) {
      a += p[i] / (x + i);
    }
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }

  function lfact(n) {
    if (n >= 0 && n <= LOG_FACT_MAX && Number.isInteger(n)) {
      return LOG_FACT_TABLE[n];
    }
    return lgamma(n + 1);
  }

  /**
   * Generalized Poisson log PMF (GP_log_pmf1 in saintpy/inference.py)
   */
  function GP_log_pmf1(x, mean, lam2) {
    if (mean <= 0.0 || x < 0) {
      return -1e12;
    }
    const lam1 = mean * (1.0 - lam2);
    const t = lam1 + lam2 * x;
    if (lam1 <= 0.0 || t <= 0.0) {
      return -1e12;
    }
    if (x === 0) {
      return -lam1;
    }
    return Math.log(lam1) + (x - 1) * Math.log(t) - t - lfact(x);
  }

  /**
   * Sample variance (divide by n - 1), matching var1_vec in saintpy/inference.py
   */
  function var1_vec(v) {
    const n = v.length;
    if (n <= 1) return 0.0;
    let sum = 0.0;
    for (let i = 0; i < n; i++) sum += v[i];
    const m = sum / n;
    let ss = 0.0;
    for (let i = 0; i < n; i++) {
      const d = v[i] - m;
      ss += d * d;
    }
    return ss / (n - 1.0);
  }

  /**
   * MLE variance (divide by n, ignoring NaNs), matching _var_mle in saintpy/stats.py
   */
  function _var_mle(v) {
    const valid = [];
    for (let i = 0; i < v.length; i++) {
      if (!Number.isNaN(v[i])) valid.push(v[i]);
    }
    if (valid.length === 0) return 0.0;
    let sum = 0.0;
    for (let i = 0; i < valid.length; i++) sum += valid[i];
    const m = sum / valid.length;
    let ss = 0.0;
    for (let i = 0; i < valid.length; i++) {
      const d = valid[i] - m;
      ss += d * d;
    }
    return ss / valid.length;
  }

  /**
   * Golden-section 1D maximizer, matching _maximize_1d in saintpy/inference.py
   */
  function _maximize_1d(f, lo, hi, tol = 1e-4, maxiter = 100) {
    const phi = (1 + Math.sqrt(5)) / 2;
    const invphi = 1 / phi;
    let a = Number(lo);
    let b = Number(hi);
    let c = b - invphi * (b - a);
    let d = a + invphi * (b - a);
    let fc = f(c);
    let fd = f(d);
    let it = 0;
    while (b - a > tol && it < maxiter) {
      it += 1;
      if (fc >= fd) {
        b = d;
        d = c;
        fd = fc;
        c = b - invphi * (b - a);
        fc = f(c);
      } else {
        a = c;
        c = d;
        fc = fd;
        d = a + invphi * (b - a);
        fd = f(d);
      }
    }
    const xbest = fc >= fd ? c : d;
    const fbest = fc >= fd ? fc : fd;
    return [xbest, fbest];
  }

  /**
   * Gaussian log-pdf with 5-sigma clipping, matching _quant_log_pdf in saintpy/stats.py
   */
  function _quant_log_pdf(x, mu = 0.0, sigma = 1.0) {
    if (Number.isNaN(x)) {
      throw new Error("quant_log_pdf received NaN input");
    }
    if (sigma <= 0) {
      return -1e12;
    }
    if (Math.abs(x - mu) / sigma > 5) {
      const sign = x < 0 || Object.is(x, -0) ? -1.0 : 1.0;
      x = sign * 5 * sigma + mu;
    }
    return -0.5 * Math.log(2 * Math.PI * sigma * sigma) - 0.5 * Math.pow((x - mu) / sigma, 2);
  }

  /**
   * Simple linear regression (intercept, slope), matching _slr in saintpy/stats.py
   */
  function _slr(x, y) {
    const n = x.length;
    if (n === 0) return [0.0, 0.0];
    let x_sum = 0.0;
    let y_sum = 0.0;
    for (let i = 0; i < n; i++) {
      x_sum += x[i];
      y_sum += y[i];
    }
    const x_bar = x_sum / n;
    const y_bar = y_sum / n;
    let beta_num = 0.0;
    let beta_den = 0.0;
    for (let i = 0; i < n; i++) {
      beta_num += (x[i] - x_bar) * (y[i] - y_bar);
      beta_den += (x[i] - x_bar) * (x[i] - x_bar);
    }
    const beta = beta_den !== 0 ? beta_num / beta_den : 0.0;
    const alpha = y_bar - beta * x_bar;
    return [alpha, beta];
  }

  /**
   * _mean_sd_estimate in saintpy/stats.py
   */
  function _mean_sd_estimate(data, t, sd_NA, lse) {
    const detected_vals = [];
    for (let i = 0; i < data.length; i++) {
      if (!Number.isNaN(data[i])) {
        detected_vals.push(data[i]);
      }
    }
    const detected = detected_vals.length;
    if (detected === 0) {
      return [t - Math.log(1.0), sd_NA];
    }
    let sum = 0.0;
    let sq_sum = 0.0;
    for (let i = 0; i < detected; i++) {
      sum += detected_vals[i];
      sq_sum += detected_vals[i] * detected_vals[i];
    }
    const mu = sum / detected;
    if (detected >= 3) {
      const sd = Math.sqrt(Math.max(0.0, sq_sum / detected - mu * mu));
      return [mu, sd];
    }
    const sd = Math.exp(mu * lse[1] + lse[0]);
    return [mu, sd];
  }

  // ============================================================================
  // 2. SPC MODEL & INFERENCE (saintpy/stats.py & saintpy/inference.py)
  // ============================================================================

  function llik_MRF_gamma_0_spc(model, b1) {
    const beta0_ = 0.0;
    const beta1_ = b1;
    let loglik = 0.0;
    const MRFtrue = Math.exp(beta1_);
    const MRFfalse = Math.exp(beta0_);
    for (let i = 0; i < model.nprey; i++) {
      for (let j = 0; j < model.nbait; j++) {
        const y = model.test_mat1[i][j];
        let allZero = true;
        for (let r = 0; r < y.length; r++) {
          if (y[r] !== 0) {
            allZero = false;
            break;
          }
        }
        if (allZero) continue;
        const k = model.Z[i][j];
        const m = model.n_rep_vec[j];
        if (m === 0) continue;
        let prod = 1.0;
        for (let rep = 0; rep < m; rep++) {
          prod *= (k[rep] ? MRFtrue : MRFfalse) / (MRFtrue + MRFfalse);
        }
        loglik += Math.log(prod) / m;
      }
    }
    return loglik;
  }

  function llikelihood_spc(model) {
    let loglik = 0.0;
    for (let i = 0; i < model.nprey; i++) {
      for (let j = 0; j < model.nbait; j++) {
        const y = model.test_mat1[i][j];
        let allZero = true;
        for (let r = 0; r < y.length; r++) {
          if (y[r] !== 0) {
            allZero = false;
            break;
          }
        }
        if (allZero) continue;
        const k = model.Z[i][j];
        const m = model.n_rep_vec[j];
        if (m === 0) continue;
        let gsum = 0.0;
        if (model.gamma !== 0.0) {
          for (const l of model.p2p_mapping[i]) {
            let s = 0;
            for (let r = 0; r < model.Z[l][j].length; r++) s += model.Z[l][j][r] ? 1 : 0;
            gsum += s / model.Z[l][j].length;
          }
        }
        const MRFtrue = Math.exp(model.beta1 + model.gamma * gsum);
        const MRFfalse = Math.exp(model.beta0);
        let prod = 1.0;
        for (let rep = 0; rep < m; rep++) {
          const ll_true = GP_log_pmf1(
            Math.min(y[rep], model.eta[i] + model.d[i]),
            model.eta[i] + model.d[i],
            model.lambda2_true[i]
          );
          const ll_false = GP_log_pmf1(
            Math.max(y[rep], model.eta[i]),
            model.eta[i],
            model.lambda2_false[i]
          );
          const num = k[rep] ? MRFtrue : MRFfalse;
          const ll = k[rep] ? ll_true : ll_false;
          prod *= (num / (MRFtrue + MRFfalse)) * Math.exp(ll);
        }
        loglik += Math.log(prod) / m;
      }
    }
    return loglik;
  }

  function loglikelihood_Z_spc(model, i, j, rep, pre_calc_loglik) {
    const k = model.Z[i][j];
    let gsum = 0.0;
    if (model.gamma !== 0.0) {
      for (const l of model.p2p_mapping[i]) {
        let s = 0;
        for (let r = 0; r < model.Z[l][j].length; r++) s += model.Z[l][j][r] ? 1 : 0;
        gsum += s / model.Z[l][j].length;
      }
    }
    const logMRFtrue = model.beta1 + model.gamma * gsum;
    const logMRFfalse = model.beta0;
    const pcl_true = pre_calc_loglik[i][j][rep][0];
    const pcl_false = pre_calc_loglik[i][j][rep][1];
    return k[rep] ? logMRFtrue + pcl_true : logMRFfalse + pcl_false;
  }

  function icm_Z_spc(model) {
    const { nprey, nbait } = model;
    const pre_calc_loglik = new Array(nprey);
    for (let i = 0; i < nprey; i++) {
      pre_calc_loglik[i] = new Array(nbait);
      for (let j = 0; j < nbait; j++) {
        const y = model.test_mat1[i][j];
        const m = model.n_rep_vec[j];
        const arr = new Array(m);
        for (let rep = 0; rep < m; rep++) {
          arr[rep] = [
            GP_log_pmf1(
              Math.min(y[rep], model.eta[i] + model.d[i]),
              model.eta[i] + model.d[i],
              model.lambda2_true[i]
            ),
            GP_log_pmf1(
              Math.max(y[rep], model.eta[i]),
              model.eta[i],
              model.lambda2_false[i]
            ),
          ];
        }
        pre_calc_loglik[i][j] = arr;
      }
    }

    for (let i = 0; i < nprey; i++) {
      for (let j = 0; j < nbait; j++) {
        const m = model.n_rep_vec[j];
        for (let rep = 0; rep < m; rep++) {
          const first = loglikelihood_Z_spc(model, i, j, rep, pre_calc_loglik);
          model.Z[i][j][rep] = !model.Z[i][j][rep];
          const second = loglikelihood_Z_spc(model, i, j, rep, pre_calc_loglik);
          if (first > second) {
            model.Z[i][j][rep] = !model.Z[i][j][rep];
          }
        }
      }
    }
  }

  function wrt_MRF_gamma_0_spc(model) {
    const oldbeta1 = model.beta1;
    const oldf = llik_MRF_gamma_0_spc(model, oldbeta1);
    const [xopt, maxf] = _maximize_1d((b1) => llik_MRF_gamma_0_spc(model, b1), -15, 15, 1e-4, 100);
    model.beta1 = maxf >= oldf ? xopt : oldbeta1;
  }

  function statModel_spc(p2p_mapping, ubait, test_mat_DATA, ctrl_mat_DATA, ip_idx_to_bait_no, nprey, nbait, opts) {
    const bait_no_to_ip_idxes = Array.from({ length: ubait.length }, () => []);
    for (let i = 0; i < ubait.length; i++) {
      for (let j = 0; j < ip_idx_to_bait_no.length; j++) {
        if (ip_idx_to_bait_no[j] === i) {
          bait_no_to_ip_idxes[i].push(j);
        }
      }
    }
    const n_rep_vec = bait_no_to_ip_idxes.map((reps) => reps.length);
    const n_ctrl_ip = ctrl_mat_DATA[0] ? ctrl_mat_DATA[0].length : 0;
    const L_used = Math.min(n_ctrl_ip, opts.L);

    const eta = new Float64Array(nprey);
    const d = new Float64Array(nprey);

    for (let i = 0; i < nprey; i++) {
      const ctrl = Array.from(ctrl_mat_DATA[i], Number);
      ctrl.sort((a, b) => b - a);
      const topL = ctrl.slice(0, L_used);
      let tmp = 0.0;
      if (topL.length > 0) {
        let s = 0.0;
        for (let k = 0; k < topL.length; k++) s += topL[k];
        tmp = s / topL.length;
      }
      if (tmp < 0.1) tmp = 0.1;
      eta[i] = tmp;
      d[i] = 4.0 * eta[i];
    }

    const test_mat = Array.from({ length: nprey }, () => new Array(ubait.length));
    const test_mat1 = Array.from({ length: nprey }, () => new Array(ubait.length));
    const test_SS = Array.from({ length: nprey }, () => new Int32Array(ubait.length));

    for (let i = 0; i < nprey; i++) {
      for (let j = 0; j < ubait.length; j++) {
        const reps = bait_no_to_ip_idxes[j];
        const v = new Int32Array(reps.length);
        const v1 = new Int32Array(reps.length);
        let ss = 0;
        for (let k = 0; k < reps.length; k++) {
          const val = Math.trunc(test_mat_DATA[i][reps[k]]);
          v[k] = val;
          v1[k] = val;
          ss += val;
        }
        test_mat[i][j] = v;
        test_mat1[i][j] = v1;
        test_SS[i][j] = ss;
      }
    }

    const Z = Array.from({ length: nprey }, () => new Array(nbait));
    for (let i = 0; i < nprey; i++) {
      for (let j = 0; j < nbait; j++) {
        const m = n_rep_vec[j];
        const z_ij = new Array(m);
        for (let rep = 0; rep < m; rep++) {
          z_ij[rep] = test_mat1[i][j][rep] > 2.0 * eta[i];
        }
        Z[i][j] = z_ij;
      }
    }

    const ctrl_mean = new Float64Array(nprey);
    for (let i = 0; i < nprey; i++) {
      let s = 0.0;
      for (let k = 0; k < n_ctrl_ip; k++) s += ctrl_mat_DATA[i][k];
      ctrl_mean[i] = n_ctrl_ip > 0 ? s / n_ctrl_ip : 0.0;
    }

    const apply_MRF = new Array(nprey).fill(false);
    const n_test_ip = test_mat_DATA[0] ? test_mat_DATA[0].length : 0;
    for (let i = 0; i < nprey; i++) {
      let s = 0;
      for (let j = 0; j < nbait; j++) {
        for (let r = 0; r < Z[i][j].length; r++) {
          if (Z[i][j][r]) s += 1;
        }
      }
      if (s < n_test_ip * opts.f) {
        apply_MRF[i] = true;
      }
    }

    const lambda2_false = new Float64Array(nprey).fill(0.1);
    for (let i = 0; i < nprey; i++) {
      const ctrl = Array.from(ctrl_mat_DATA[i], (x) => Math.trunc(x));
      ctrl.sort((a, b) => b - a);
      const topL = ctrl.slice(0, L_used);
      const variance = var1_vec(topL);
      if (variance > eta[i]) {
        lambda2_false[i] = 1.0 - Math.sqrt(eta[i] / variance);
      }
    }
    const lambda2_true = new Float64Array(nprey).fill(0.1);

    return {
      nprey,
      nbait,
      test_mat_DATA,
      test_mat,
      test_mat1,
      test_SS,
      n_rep_vec,
      p2p_mapping,
      ctrl_mean,
      n_ctrl_ip,
      ctrl_mat_DATA,
      apply_MRF,
      Z,
      beta0: 0.0,
      beta1: 0.0,
      gamma: 0.0,
      eta,
      d,
      lambda2_true,
      lambda2_false,
      opts,
      llikelihood() {
        return llikelihood_spc(this);
      },
    };
  }

  function calculateScore_spc(model) {
    const { nprey, nbait } = model;
    const average_score = Array.from({ length: nprey }, () => new Float64Array(nbait));
    const min_log_odds_score = Array.from({ length: nprey }, () => new Float64Array(nbait));
    const maximum_score = Array.from({ length: nprey }, () => new Float64Array(nbait));
    const rep_scores_mat = Array.from({ length: nprey }, () => new Array(nbait));

    for (let i = 0; i < nprey; i++) {
      for (let j = 0; j < nbait; j++) {
        const y = model.test_mat1[i][j];
        let allZero = true;
        let maxY = -Infinity;
        for (let r = 0; r < y.length; r++) {
          if (y[r] !== 0) allZero = false;
          if (y[r] > maxY) maxY = y[r];
        }
        if (allZero) continue;

        const m = model.n_rep_vec[j];
        let gsum = 0.0;
        for (const l of model.p2p_mapping[i]) {
          let s = 0;
          for (let r = 0; r < model.Z[l][j].length; r++) s += model.Z[l][j][r] ? 1 : 0;
          gsum += s / model.Z[l][j].length;
        }

        const tmp_scores = new Array(m).fill(0.0);
        const tmp_odds_scores = new Array(m).fill(0.0);

        for (let rep = 0; rep < m; rep++) {
          const yrep = y[rep];
          let tmp_mean;
          if (yrep >= 5) {
            tmp_mean = Math.min(model.eta[i] + model.d[i], 5.0 * model.eta[i]);
          } else {
            tmp_mean = model.eta[i] + model.d[i];
          }

          const log_true = model.beta1 + model.gamma * gsum;
          const log_false =
            model.beta0 +
            GP_log_pmf1(Math.max(yrep, model.eta[i]), model.eta[i], model.lambda2_false[i]) -
            GP_log_pmf1(Math.min(yrep, tmp_mean), tmp_mean, model.lambda2_true[i]);

          tmp_odds_scores[rep] = log_true - log_false;

          if (yrep <= 1 || yrep <= model.ctrl_mean[i]) {
            tmp_scores[rep] = 0.0;
          } else {
            const delta = log_false - log_true;
            if (delta > 50) {
              tmp_scores[rep] = 0.0;
            } else if (delta < -50) {
              tmp_scores[rep] = 1.0;
            } else {
              tmp_scores[rep] = 1.0 / (1.0 + Math.exp(delta));
            }
          }
        }

        rep_scores_mat[i][j] = [...tmp_scores];
        tmp_scores.sort((a, b) => b - a);
        tmp_odds_scores.sort((a, b) => b - a);

        const max_rep = model.opts.R;
        if (m > max_rep) {
          let s = 0.0;
          for (let r = 0; r < max_rep; r++) s += tmp_scores[r];
          average_score[i][j] = s / max_rep;
        } else {
          let s = 0.0;
          for (let r = 0; r < m; r++) s += tmp_scores[r];
          average_score[i][j] = s / m;
        }

        min_log_odds_score[i][j] = tmp_odds_scores[tmp_odds_scores.length - 1];
        maximum_score[i][j] = m > 0 ? tmp_scores[0] : 0.0;

        if (maxY === 1) {
          average_score[i][j] = 0.0;
          maximum_score[i][j] = 0.0;
        }
      }
    }

    model._rep_scores_mat = rep_scores_mat;
    return [average_score, maximum_score, min_log_odds_score];
  }

  // ============================================================================
  // 3. INT MODEL & INFERENCE (saintpy/stats.py)
  // ============================================================================

  function statModel_int(p2p_mapping, ubait, test_mat_DATA, ctrl_mat_DATA, ip_idx_to_bait_no, nprey, nbait, opts) {
    const bait_no_to_ip_idxes = Array.from({ length: ubait.length }, () => []);
    for (let i = 0; i < ubait.length; i++) {
      for (let j = 0; j < ip_idx_to_bait_no.length; j++) {
        if (ip_idx_to_bait_no[j] === i) {
          bait_no_to_ip_idxes[i].push(j);
        }
      }
    }
    const n_rep_vec = bait_no_to_ip_idxes.map((reps) => reps.length);
    const n_ctrl_ip = ctrl_mat_DATA[0] ? ctrl_mat_DATA[0].length : 0;

    const ctrls_with_NA = [];
    const all_ctrl_valid = [];
    for (let i = 0; i < nprey; i++) {
      const row = ctrl_mat_DATA[i];
      let hasNA = false;
      for (let k = 0; k < row.length; k++) {
        if (Number.isNaN(row[k])) {
          hasNA = true;
        } else {
          all_ctrl_valid.push(row[k]);
        }
      }
      if (hasNA) {
        for (let k = 0; k < row.length; k++) {
          if (!Number.isNaN(row[k])) ctrls_with_NA.push(row[k]);
        }
      }
    }
    let t = 0.0;
    if (ctrls_with_NA.length > 0) {
      let s = 0.0;
      for (let k = 0; k < ctrls_with_NA.length; k++) s += ctrls_with_NA[k];
      t = s / ctrls_with_NA.length;
    } else if (all_ctrl_valid.length > 0) {
      let s = 0.0;
      for (let k = 0; k < all_ctrl_valid.length; k++) s += all_ctrl_valid[k];
      t = s / all_ctrl_valid.length;
    }

    const eta = new Float64Array(nprey);
    const d = new Float64Array(nprey);
    const sd_ctrl = new Float64Array(nprey);

    const log_sd_ctrls = [];
    const mean_ctrls = [];
    for (let idx = 0; idx < nprey; idx++) {
      const row = ctrl_mat_DATA[idx];
      let allValid = row.length > 0;
      for (let k = 0; k < row.length; k++) {
        if (Number.isNaN(row[k])) {
          allValid = false;
          break;
        }
      }
      if (allValid) {
        const variance = _var_mle(row);
        if (variance <= 0 || Number.isNaN(variance)) continue;
        log_sd_ctrls.push(Math.log(Math.sqrt(variance)));
        let s = 0.0;
        for (let k = 0; k < row.length; k++) s += row[k];
        mean_ctrls.push(s / row.length);
      }
    }

    const lse = log_sd_ctrls.length > 0 ? _slr(mean_ctrls, log_sd_ctrls) : [0.0, 0.0];
    const sorted_log_sd_ctrls = [...log_sd_ctrls].sort((a, b) => a - b);
    const median_sd =
      sorted_log_sd_ctrls.length > 0
        ? Math.exp(sorted_log_sd_ctrls[Math.floor(sorted_log_sd_ctrls.length / 2)])
        : 1.0;
    const sd_NA = median_sd;

    for (let idx = 0; idx < nprey; idx++) {
      const ctrl = ctrl_mat_DATA[idx];
      const [mu, sd_val] = _mean_sd_estimate(ctrl, t, sd_NA, lse);
      eta[idx] = mu;
      sd_ctrl[idx] = Math.max(sd_NA, sd_val);
    }

    const test_mat = Array.from({ length: nprey }, () => new Array(ubait.length));
    for (let i = 0; i < nprey; i++) {
      for (let j = 0; j < ubait.length; j++) {
        const reps = bait_no_to_ip_idxes[j];
        const v = new Float64Array(reps.length);
        for (let k = 0; k < reps.length; k++) {
          v[k] = test_mat_DATA[i][reps[k]];
        }
        test_mat[i][j] = v;
      }
    }
    const test_mat1 = test_mat;
    const sd_false = new Float64Array(sd_ctrl);

    const n_test_ip = test_mat_DATA[0] ? test_mat_DATA[0].length : 0;
    for (let i = 0; i < nprey; i++) {
      const vals = [];
      for (let ip = 0; ip < n_test_ip; ip++) {
        const val = test_mat_DATA[i][ip];
        if (!Number.isNaN(val) && eta[i] < val) {
          vals.push(val);
        }
      }
      let mu = eta[i];
      if (vals.length > 0) {
        let s = 0.0;
        for (let k = 0; k < vals.length; k++) s += vals[k];
        mu = s / vals.length;
      }
      mu = Math.max(mu, eta[i]);
      d[i] = mu - eta[i];
      if (d[i] < Math.log(4.0)) {
        d[i] = Math.log(4.0);
      }
    }

    const Z = Array.from({ length: nprey }, () => new Array(nbait));
    for (let i = 0; i < nprey; i++) {
      for (let j = 0; j < nbait; j++) {
        const m = n_rep_vec[j];
        const z_ij = new Array(m);
        for (let rep = 0; rep < m; rep++) {
          const y = test_mat1[i][j][rep];
          z_ij[rep] = Number.isNaN(y) ? false : y > eta[i] + Math.log(5.0);
        }
        Z[i][j] = z_ij;
      }
    }

    const sd_true = new Float64Array(nprey);
    for (let i = 0; i < nprey; i++) {
      const Z1_vals = [];
      for (let j = 0; j < nbait; j++) {
        for (let rep = 0; rep < n_rep_vec[j]; rep++) {
          if (Z[i][j][rep]) {
            Z1_vals.push(test_mat1[i][j][rep]);
          }
        }
      }
      if (Z1_vals.length >= 2) {
        sd_true[i] = Math.sqrt(var1_vec(Z1_vals));
      } else {
        sd_true[i] = sd_false[i];
      }
      sd_true[i] = Math.max(sd_NA, sd_true[i]);
    }

    return {
      nprey,
      nbait,
      test_mat_DATA,
      test_mat,
      test_mat1,
      n_rep_vec,
      p2p_mapping,
      t,
      n_ctrl_ip,
      Z,
      beta0: 0.0,
      beta1: 0.0,
      gamma: 0.0,
      eta,
      d,
      sd_true,
      sd_false,
      opts,
      llikelihood() {
        let loglik = 0.0;
        for (let i = 0; i < this.nprey; i++) {
          for (let j = 0; j < this.nbait; j++) {
            const y = this.test_mat1[i][j];
            let allNaN = true;
            for (let r = 0; r < y.length; r++) {
              if (!Number.isNaN(y[r])) {
                allNaN = false;
                break;
              }
            }
            if (allNaN) continue;
            const k = this.Z[i][j];
            const m = this.n_rep_vec[j];
            if (m === 0) continue;
            let gsum = 0.0;
            for (const l of this.p2p_mapping[i]) {
              let s = 0;
              for (let r = 0; r < this.Z[l][j].length; r++) s += this.Z[l][j][r] ? 1 : 0;
              gsum += s / this.Z[l][j].length;
            }
            const MRFtrue = Math.exp(this.beta1 + this.gamma * gsum);
            const MRFfalse = Math.exp(this.beta0);
            let prod = 1.0;
            for (let rep = 0; rep < m; rep++) {
              const yrep = Number.isNaN(y[rep]) ? this.t : y[rep];
              const ll_true = _quant_log_pdf(
                Math.min(yrep, this.eta[i] + this.d[i]),
                this.eta[i] + this.d[i],
                this.sd_true[i]
              );
              const ll_false = _quant_log_pdf(Math.max(yrep, this.eta[i]), this.eta[i], this.sd_false[i]);
              const num = k[rep] ? MRFtrue : MRFfalse;
              const ll = k[rep] ? ll_true : ll_false;
              prod *= (num / (MRFtrue + MRFfalse)) * Math.exp(ll);
            }
            loglik += Math.log(prod) / m;
          }
        }
        return loglik;
      },
      loglikelihood_Z(i, j, rep, pre_calc_loglik) {
        const k = this.Z[i][j];
        let gsum = 0.0;
        if (this.gamma !== 0.0) {
          for (const l of this.p2p_mapping[i]) {
            let s = 0;
            for (let r = 0; r < this.Z[l][j].length; r++) s += this.Z[l][j][r] ? 1 : 0;
            gsum += s / this.Z[l][j].length;
          }
        }
        const logMRFtrue = this.beta1 + this.gamma * gsum;
        const logMRFfalse = this.beta0;
        const pcl_true = pre_calc_loglik[i][j][rep][0];
        const pcl_false = pre_calc_loglik[i][j][rep][1];
        return k[rep] ? logMRFtrue + pcl_true : logMRFfalse + pcl_false;
      },
      icm_Z() {
        const pre_calc_loglik = new Array(this.nprey);
        for (let i = 0; i < this.nprey; i++) {
          pre_calc_loglik[i] = new Array(this.nbait);
          for (let j = 0; j < this.nbait; j++) {
            const y = this.test_mat1[i][j];
            const m = this.n_rep_vec[j];
            const arr = new Array(m);
            for (let rep = 0; rep < m; rep++) {
              const yrep = Number.isNaN(y[rep]) ? this.t : y[rep];
              arr[rep] = [
                _quant_log_pdf(
                  Math.min(yrep, this.eta[i] + this.d[i]),
                  this.eta[i] + this.d[i],
                  this.sd_true[i]
                ),
                _quant_log_pdf(Math.max(yrep, this.eta[i]), this.eta[i], this.sd_false[i]),
              ];
            }
            pre_calc_loglik[i][j] = arr;
          }
        }
        for (let i = 0; i < this.nprey; i++) {
          for (let j = 0; j < this.nbait; j++) {
            const m = this.n_rep_vec[j];
            for (let rep = 0; rep < m; rep++) {
              const first = this.loglikelihood_Z(i, j, rep, pre_calc_loglik);
              this.Z[i][j][rep] = !this.Z[i][j][rep];
              const second = this.loglikelihood_Z(i, j, rep, pre_calc_loglik);
              if (first > second) {
                this.Z[i][j][rep] = !this.Z[i][j][rep];
              }
            }
          }
        }
      },
      llik_MRF_gamma_0(b1) {
        const beta0_ = 0.0;
        const beta1_ = b1;
        let loglik = 0.0;
        const MRFtrue = Math.exp(beta1_);
        const MRFfalse = Math.exp(beta0_);
        for (let i = 0; i < this.nprey; i++) {
          for (let j = 0; j < this.nbait; j++) {
            const y = this.test_mat1[i][j];
            let allNaN = true;
            for (let r = 0; r < y.length; r++) {
              if (!Number.isNaN(y[r])) {
                allNaN = false;
                break;
              }
            }
            if (allNaN) continue;
            const k = this.Z[i][j];
            const m = this.n_rep_vec[j];
            if (m === 0) continue;
            let prod = 1.0;
            for (let rep = 0; rep < m; rep++) {
              prod *= (k[rep] ? MRFtrue : MRFfalse) / (MRFtrue + MRFfalse);
            }
            loglik += Math.log(prod) / m;
          }
        }
        return loglik;
      },
      wrt_MRF() {
        const old_beta1 = this.beta1;
        const oldf = this.llik_MRF_gamma_0(old_beta1);
        const [xopt, maxf] = _maximize_1d((b1) => this.llik_MRF_gamma_0(b1), -15.0, 15.0, 1e-8, 250);
        if (maxf >= oldf) {
          this.beta1 = xopt;
        } else {
          this.beta1 = old_beta1;
        }
      },
      calculateScore() {
        const average_score = Array.from({ length: this.nprey }, () => new Float64Array(this.nbait));
        const min_log_odds_score = Array.from({ length: this.nprey }, () => new Float64Array(this.nbait));
        const maximum_score = Array.from({ length: this.nprey }, () => new Float64Array(this.nbait));
        const rep_scores_mat = Array.from({ length: this.nprey }, () => new Array(this.nbait));

        for (let i = 0; i < this.nprey; i++) {
          for (let j = 0; j < this.nbait; j++) {
            const y = this.test_mat1[i][j];
            const m = this.n_rep_vec[j];
            let gsum = 0.0;
            for (const l of this.p2p_mapping[i]) {
              let s = 0;
              for (let r = 0; r < this.Z[l][j].length; r++) s += this.Z[l][j][r] ? 1 : 0;
              gsum += s / this.Z[l][j].length;
            }
            const MRFtrue = Math.exp(this.beta1 + this.gamma * gsum);
            const MRFfalse = Math.exp(this.beta0);
            const tmp_scores = new Array(m).fill(0.0);
            const tmp_odds_scores = new Array(m).fill(0.0);
            for (let rep = 0; rep < m; rep++) {
              const yrep = Number.isNaN(y[rep]) ? this.t : y[rep];
              const unnorm_true =
                MRFtrue *
                Math.exp(
                  _quant_log_pdf(
                    Math.min(yrep, this.eta[i] + this.d[i]),
                    this.eta[i] + this.d[i],
                    this.sd_true[i]
                  )
                );
              const unnorm_false =
                MRFfalse *
                Math.exp(_quant_log_pdf(Math.max(yrep, this.eta[i]), this.eta[i], this.sd_false[i]));
              tmp_scores[rep] = unnorm_true / (unnorm_true + unnorm_false);
              if (Number.isNaN(y[rep])) {
                tmp_scores[rep] = 0.0;
              }
              tmp_odds_scores[rep] = Math.log(unnorm_true) - Math.log(unnorm_false);
            }
            rep_scores_mat[i][j] = [...tmp_scores];
            tmp_scores.sort((a, b) => b - a);
            tmp_odds_scores.sort((a, b) => b - a);
            const max_rep = this.opts.R;
            if (m > max_rep) {
              let s = 0.0;
              for (let r = 0; r < max_rep; r++) s += tmp_scores[r];
              average_score[i][j] = s / max_rep;
            } else {
              let s = 0.0;
              for (let r = 0; r < m; r++) s += tmp_scores[r];
              average_score[i][j] = s / m;
            }
            min_log_odds_score[i][j] = tmp_odds_scores[tmp_odds_scores.length - 1];
            maximum_score[i][j] = m > 0 ? tmp_scores[0] : 0.0;
          }
        }
        this._rep_scores_mat = rep_scores_mat;
        return [average_score, maximum_score, min_log_odds_score];
      },
    };
  }

  // ============================================================================
  // 4. DELIMITED FILE PARSING & SAINT INPUT BUILDER
  // ============================================================================

  function parseDelimitedText(text, forceHeaderlessType = null) {
    const lines = String(text || "")
      .replace(/^\uFEFF/, "")
      .replace(/\r\n/g, "\n")
      .replace(/\r/g, "\n")
      .split("\n");
    const nonEmpty = [];
    for (let i = 0; i < lines.length; i++) {
      if (lines[i].trim().length > 0) nonEmpty.push(lines[i]);
    }
    if (nonEmpty.length === 0) return { headers: [], rows: [], delimiter: "," };

    const firstLine = nonEmpty[0];
    const tabCount = (firstLine.match(/\t/g) || []).length;
    const commaCount = (firstLine.match(/,/g) || []).length;
    const delimiter = tabCount > commaCount ? "\t" : ",";

    function splitRow(line) {
      if (delimiter === "\t" && !line.includes('"')) {
        return line.split("\t");
      }
      const result = [];
      let cur = "";
      let inQuotes = false;
      for (let i = 0; i < line.length; i++) {
        const ch = line[i];
        if (ch === '"') {
          if (inQuotes && line[i + 1] === '"') {
            cur += '"';
            i++;
          } else {
            inQuotes = !inQuotes;
          }
        } else if (ch === delimiter && !inQuotes) {
          result.push(cur);
          cur = "";
        } else {
          cur += ch;
        }
      }
      result.push(cur);
      return result;
    }

    const firstTokens = splitRow(nonEmpty[0]).map((h) => h.trim());

    // Detect classic headerless SAINTexpress files if needed
    let syntheticHeaders = null;
    if (forceHeaderlessType === "bait") {
      syntheticHeaders = ["ip_id", "bait_id", "is_control"];
    } else if (forceHeaderlessType === "prey") {
      syntheticHeaders = ["prey_id", "length", "gene_name"];
    } else if (forceHeaderlessType === "inter") {
      syntheticHeaders =
        firstTokens.length >= 4
          ? ["ip_id", "bait_id", "prey_id", "quant"]
          : ["ip_id", "prey_id", "quant"];
    } else {
      // Auto-detect headerless SAINT files if first line has no known header names
      const lowerTokens = firstTokens.map((t) => t.toLowerCase());
      const hasKnownHeader = lowerTokens.some((t) =>
        [
          "ip_id",
          "sample_id",
          "bait_id",
          "prey_id",
          "quant",
          "is_control",
          "tc",
          "protein",
          "protein id",
          "gene",
          "protid",
          "sample",
          "file",
        ].includes(t)
      );
      if (!hasKnownHeader && nonEmpty.length > 1) {
        if (
          firstTokens.length === 3 &&
          ["t", "c", "true", "false", "0", "1"].includes(lowerTokens[2])
        ) {
          syntheticHeaders = ["ip_id", "bait_id", "is_control"];
        } else if (firstTokens.length === 4 && !Number.isNaN(Number(firstTokens[3]))) {
          syntheticHeaders = ["ip_id", "bait_id", "prey_id", "quant"];
        } else if (
          firstTokens.length === 3 &&
          !Number.isNaN(Number(firstTokens[1])) &&
          Number.isNaN(Number(firstTokens[2]))
        ) {
          syntheticHeaders = ["prey_id", "length", "gene_name"];
        }
      }
    }

    const rawHeaders = syntheticHeaders || firstTokens;
    const startIdx = syntheticHeaders ? 0 : 1;
    const rows = [];
    for (let i = startIdx; i < nonEmpty.length; i++) {
      const vals = splitRow(nonEmpty[i]);
      const obj = {};
      for (let c = 0; c < rawHeaders.length; c++) {
        const key = rawHeaders[c] === "" && c === 0 ? "_index" : rawHeaders[c];
        obj[key] = vals[c] !== undefined ? vals[c].trim() : "";
      }
      rows.push(obj);
    }
    const headers = rawHeaders.map((h, c) => (h === "" && c === 0 ? "_index" : h));
    return { headers, rows, delimiter };
  }

  function coerceBool(x) {
    if (typeof x === "boolean") return x;
    if (typeof x === "number") return Boolean(Math.trunc(x));
    if (typeof x === "string") {
      const xl = x.trim().toLowerCase();
      if (["c", "control", "ctrl", "1", "true", "yes"].includes(xl)) return true;
      if (["t", "test", "0", "false", "no"].includes(xl)) return false;
      const num = Number(xl);
      if (!Number.isNaN(num)) return Boolean(Math.trunc(num));
    }
    throw new Error(`Cannot coerce is_control value: '${x}'`);
  }

  /**
   * Build input matrices matching saintpy/io.py build_inputs()
   * Fixes upstream KeyError when inter.csv also contains 'bait_id'
   * Supports selfInteractionMode:
   *  - "saintpy": exact saintPy (sets test[prey, sid] = 0 when gene_name == bait_id)
   *  - "smart": sets 0 for SPC and NaN for INT when prey_id or gene_name matches bait_id (or custom alias)
   *  - "keep": keeps self-interactions untouched
   */
  function build_inputs(
    samplesRows,
    interRows,
    preyRows,
    quant_type = "spc",
    value_col = "quant",
    options = {}
  ) {
    const selfInteractionMode = options.selfInteractionMode || "saintpy";
    const baitAliases = options.baitAliases || {}; // optional { bait_id: "GENE_OR_PREY_ID" }

    // 1. Validate & normalize samples
    const samples = [];
    const seenIp = new Set();
    const sid2sample = new Map();

    for (const r of samplesRows) {
      if (r.excluded === true || r.excluded === "true") continue;
      const ip_id = String(r.ip_id ?? r.sample_id ?? "").trim();
      const bait_id = String(r.bait_id ?? r.condition ?? "").trim();
      if (!ip_id) continue;
      let rawCtrl = r.is_control;
      if (rawCtrl === undefined || rawCtrl === "") {
        if (r.tc !== undefined && r.tc !== "") {
          rawCtrl = r.tc;
        } else {
          throw new Error("Bait/Samples table must contain 'is_control' or 'tc' column.");
        }
      }
      const is_control = coerceBool(rawCtrl);
      if (seenIp.has(ip_id)) {
        throw new Error(`Duplicate ip_id in samples: '${ip_id}'`);
      }
      seenIp.add(ip_id);
      const item = { ip_id, bait_id, is_control };
      samples.push(item);
      sid2sample.set(ip_id, item);
    }

    const numTestSamples = samples.filter((s) => !s.is_control).length;
    const numCtrlSamples = samples.filter((s) => s.is_control).length;
    if (numTestSamples === 0) {
      throw new Error("At least 1 Test sample (is_control = FALSE) is required.");
    }
    if (numCtrlSamples === 0) {
      throw new Error(
        "At least 1 Negative Control sample (is_control = TRUE) is required by SAINT to model background."
      );
    }

    // 2. Validate & normalize interactions
    const df = [];
    const unknownIps = new Set();
    for (const r of interRows) {
      const ip_id = String(r.ip_id ?? r.sample_id ?? "").trim();
      const prey_id = String(r.prey_id ?? "").trim();
      if (!ip_id || !prey_id) continue;
      if (!sid2sample.has(ip_id)) {
        // If sample was explicitly excluded by user in UI, skip its interactions gracefully;
        // otherwise if strict foreign key check is active, record unknown IP
        if (options.skipExcludedSamples) continue;
        unknownIps.add(ip_id);
        continue;
      }
      const rawVal = r[value_col] !== undefined ? r[value_col] : r.quant ?? r.spec;
      if (rawVal === undefined) {
        throw new Error(`Interactions table must have column '${value_col}'`);
      }
      let q = rawVal === "" || rawVal === null ? NaN : Number(rawVal);
      if (q < 0) {
        throw new Error(`Negative quantification value (${q}) found for ${ip_id} / ${prey_id}`);
      }
      if (quant_type === "int") {
        if (Number.isNaN(q) || q <= 0) {
          q = NaN;
        } else {
          q = Math.log(q);
        }
      }
      const samp = sid2sample.get(ip_id);
      df.push({
        ip_id,
        prey_id,
        quant: q,
        bait_id: samp.bait_id,
        is_control: samp.is_control,
      });
    }

    if (unknownIps.size > 0 && !options.skipExcludedSamples) {
      const ex = Array.from(unknownIps).sort().slice(0, 5).join(", ");
      throw new Error(`Unknown ip_id in interactions not found in bait.csv: [${ex}]`);
    }

    // 3. Prey metadata (preserving first-seen prey_id order in df)
    const preyMap = new Map();
    const preyLengthMap = new Map();
    if (Array.isArray(preyRows)) {
      for (const r of preyRows) {
        const pid = String(r.prey_id ?? "").trim();
        if (!pid) continue;
        const gname = String(r.gene_name ?? r.prey_gene ?? pid).trim() || pid;
        const len = r.length ?? r.prey_length ?? "";
        preyMap.set(pid, gname);
        if (len !== "" && len !== undefined && !Number.isNaN(Number(len))) {
          preyLengthMap.set(pid, Number(len));
        }
      }
    }

    const all_preys = [];
    const prey_index = new Map();
    const PDATA = [];
    const prey_to_gene = new Map();

    for (const row of df) {
      if (!prey_index.has(row.prey_id)) {
        const rowId = all_preys.length;
        all_preys.push(row.prey_id);
        prey_index.set(row.prey_id, rowId);
        const gene_name = preyMap.get(row.prey_id) || row.prey_id;
        prey_to_gene.set(row.prey_id, gene_name);
        PDATA.push({
          prey_id: row.prey_id,
          gene_name,
          length: preyLengthMap.get(row.prey_id) ?? null,
          get_preyId() {
            return this.prey_id;
          },
          get_preyGeneId() {
            return this.gene_name;
          },
        });
      }
    }

    // 4. Split test_df and control_df, sort column ip_ids lexicographically like pandas pivot_table
    const testIpSet = new Set();
    const ctrlIpSet = new Set();
    const ubait = [];
    const seenTestIpForUbait = new Set();
    const seenBaitForUbait = new Set();

    for (const row of df) {
      if (!row.is_control) {
        testIpSet.add(row.ip_id);
        if (!seenTestIpForUbait.has(row.ip_id)) {
          seenTestIpForUbait.add(row.ip_id);
          if (!seenBaitForUbait.has(row.bait_id)) {
            seenBaitForUbait.add(row.bait_id);
            ubait.push(row.bait_id);
          }
        }
      } else {
        ctrlIpSet.add(row.ip_id);
      }
    }

    const test_ip_ids = Array.from(testIpSet).sort();
    const ctrl_ip_ids = Array.from(ctrlIpSet).sort();
    const testIpIndex = new Map(test_ip_ids.map((id, idx) => [id, idx]));
    const ctrlIpIndex = new Map(ctrl_ip_ids.map((id, idx) => [id, idx]));

    const nprey = all_preys.length;
    const nexpr = test_ip_ids.length;
    const nctrl = ctrl_ip_ids.length;
    const nbait = ubait.length;

    const zero_fill = quant_type === "spc";
    const fill_val = zero_fill ? 0 : NaN;

    const test_mat_DATA = Array.from({ length: nprey }, () => new Float64Array(nexpr).fill(fill_val));
    const test_seen = Array.from({ length: nprey }, () => new Uint8Array(nexpr));
    const ctrl_mat_DATA = Array.from({ length: nprey }, () => new Float64Array(nctrl).fill(fill_val));
    const ctrl_seen = Array.from({ length: nprey }, () => new Uint8Array(nctrl));

    for (const row of df) {
      const pIdx = prey_index.get(row.prey_id);
      if (!row.is_control) {
        const cIdx = testIpIndex.get(row.ip_id);
        if (!Number.isNaN(row.quant)) {
          if (test_seen[pIdx][cIdx] === 0) {
            test_mat_DATA[pIdx][cIdx] = row.quant;
            test_seen[pIdx][cIdx] = 1;
          } else {
            test_mat_DATA[pIdx][cIdx] += row.quant;
          }
        }
      } else {
        const cIdx = ctrlIpIndex.get(row.ip_id);
        if (!Number.isNaN(row.quant)) {
          if (ctrl_seen[pIdx][cIdx] === 0) {
            ctrl_mat_DATA[pIdx][cIdx] = row.quant;
            ctrl_seen[pIdx][cIdx] = 1;
          } else {
            ctrl_mat_DATA[pIdx][cIdx] += row.quant;
          }
        }
      }
    }

    const bait_pos = new Map(ubait.map((b, i) => [b, i]));
    const ip_idx_to_bait_no = test_ip_ids.map((sid) => bait_pos.get(sid2sample.get(sid).bait_id));

    // 5. Self-interaction handling
    const maskedSelfInteractions = [];
    if (selfInteractionMode !== "keep") {
      for (let cIdx = 0; cIdx < test_ip_ids.length; cIdx++) {
        const sid = test_ip_ids[cIdx];
        const bait = sid2sample.get(sid).bait_id;
        const aliasTarget = baitAliases[bait] ? String(baitAliases[bait]).toLowerCase() : null;
        for (let pIdx = 0; pIdx < nprey; pIdx++) {
          const prey = all_preys[pIdx];
          const gene = prey_to_gene.get(prey) || prey;
          let isSelf = false;
          if (selfInteractionMode === "saintpy") {
            isSelf = gene === bait;
          } else {
            const bLow = bait.toLowerCase();
            isSelf =
              gene.toLowerCase() === bLow ||
              prey.toLowerCase() === bLow ||
              (aliasTarget && (gene.toLowerCase() === aliasTarget || prey.toLowerCase() === aliasTarget));
          }
          if (isSelf) {
            const oldVal = test_mat_DATA[pIdx][cIdx];
            if (selfInteractionMode === "saintpy") {
              test_mat_DATA[pIdx][cIdx] = 0;
            } else {
              test_mat_DATA[pIdx][cIdx] = quant_type === "int" ? NaN : 0;
            }
            if (!Number.isNaN(oldVal) && oldVal !== 0) {
              maskedSelfInteractions.push({ ip_id: sid, bait_id: bait, prey_id: prey, gene_name: gene, oldVal });
            }
          }
        }
      }
    }

    if (quant_type === "spc") {
      for (let pIdx = 0; pIdx < nprey; pIdx++) {
        for (let cIdx = 0; cIdx < nexpr; cIdx++) {
          test_mat_DATA[pIdx][cIdx] = Math.trunc(test_mat_DATA[pIdx][cIdx]);
        }
        for (let cIdx = 0; cIdx < nctrl; cIdx++) {
          ctrl_mat_DATA[pIdx][cIdx] = Math.trunc(ctrl_mat_DATA[pIdx][cIdx]);
        }
      }
    }

    return {
      PDATA,
      ubait,
      ip_idx_to_bait_no,
      test_mat_DATA,
      ctrl_mat_DATA,
      nprey,
      nbait,
      nexpr,
      nctrl,
      test_index: all_preys,
      test_columns: test_ip_ids,
      ctrl_columns: ctrl_ip_ids,
      maskedSelfInteractions,
    };
  }

  function _build_test_mask(nprey, nbait, test_mat_DATA, ip_idx_to_bait_no) {
    const bait_no_to_ip_idxes = Array.from({ length: nbait }, () => []);
    for (let ip_idx = 0; ip_idx < ip_idx_to_bait_no.length; ip_idx++) {
      bait_no_to_ip_idxes[ip_idx_to_bait_no[ip_idx]].push(ip_idx);
    }
    const mask = Array.from({ length: nprey }, () => new Uint8Array(nbait));
    for (let j = 0; j < nbait; j++) {
      const ip_indices = bait_no_to_ip_idxes[j];
      if (ip_indices.length === 0) continue;
      for (let i = 0; i < nprey; i++) {
        for (let r = 0; r < ip_indices.length; r++) {
          const val = test_mat_DATA[i][ip_indices[r]];
          if (!Number.isNaN(val) && val > 0) {
            mask[i][j] = 1;
            break;
          }
        }
      }
    }
    return mask;
  }

  /**
   * Compute Bayesian False Discovery Rate (BFDR) for each output row:
   * BFDR(i) = sum_{k: AvgP_k >= AvgP_i} (1 - AvgP_k) / count_{k: AvgP_k >= AvgP_i}
   */
  function computeBFDR(rows) {
    const n = rows.length;
    if (n === 0) return;
    const indices = Array.from({ length: n }, (_, i) => i);
    indices.sort((a, b) => rows[b].AvgP - rows[a].AvgP);

    let cumComp = 0.0;
    const bfdrBySortedPos = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      const p = rows[indices[k]].AvgP;
      cumComp += 1.0 - p;
      bfdrBySortedPos[k] = cumComp / (k + 1);
    }
    let k = 0;
    while (k < n) {
      let end = k;
      while (end + 1 < n && rows[indices[end + 1]].AvgP === rows[indices[k]].AvgP) {
        end++;
      }
      const tieBfdr = bfdrBySortedPos[end];
      for (let m = k; m <= end; m++) {
        rows[indices[m]].computedBFDR = tieBfdr;
      }
      k = end + 1;
    }
  }

  function build_output_rows(
    model,
    avgp,
    maxp,
    min_logodds,
    ip_idx_to_bait_no,
    test_mat_DATA,
    ctrl_mat_DATA,
    test_mat_mask,
    quant_type = "spc",
    test_columns = [],
    ctrl_columns = []
  ) {
    const ubait = model.ubait;
    const PDATA = model.PDATA;
    const nprey = test_mat_DATA.length;
    const nbait = ubait.length;

    const bait_no_to_ip_idxes = Array.from({ length: nbait }, () => []);
    for (let ip_idx = 0; ip_idx < ip_idx_to_bait_no.length; ip_idx++) {
      bait_no_to_ip_idxes[ip_idx_to_bait_no[ip_idx]].push(ip_idx);
    }

    const rows = [];
    for (let i = 0; i < nprey; i++) {
      const gene_name = PDATA[i].get_preyGeneId();
      const prey_id = PDATA[i].get_preyId();
      const prey_length = PDATA[i].length;
      for (let bait_idx = 0; bait_idx < nbait; bait_idx++) {
        if (test_mat_mask && test_mat_mask[i][bait_idx] === 0) {
          continue;
        }
        const avg_p = avgp[i][bait_idx];
        const max_p = maxp[i][bait_idx];
        const log_odds = min_logodds[i][bait_idx];
        const rep_indices = bait_no_to_ip_idxes[bait_idx];
        const spec = rep_indices.map((idx) => test_mat_DATA[i][idx]);
        const rep_ip_names = rep_indices.map((idx) => test_columns[idx] || `Rep_${idx + 1}`);

        let spec_str;
        if (quant_type === "spc") {
          spec_str = spec.map((x) => (Number.isNaN(x) ? "" : String(Math.trunc(x)))).join(",");
        } else {
          spec_str = spec.map((x) => (Number.isNaN(x) ? "" : Number(x).toFixed(4))).join(",");
        }

        const valid_spec = spec.filter((x) => !Number.isNaN(x));
        let spec_sum = 0.0;
        for (let r = 0; r < valid_spec.length; r++) spec_sum += valid_spec[r];
        const avg_spec = valid_spec.length > 0 ? spec_sum / valid_spec.length : 0.0;
        const num_reps = rep_indices.length;

        const ctrl_counts = Array.from(ctrl_mat_DATA[i]);
        let ctrl_str;
        if (quant_type === "spc") {
          ctrl_str = ctrl_counts.map((x) => (Number.isNaN(x) ? "" : String(Math.trunc(x)))).join(",");
        } else {
          ctrl_str = ctrl_counts.map((x) => (Number.isNaN(x) ? "" : Number(x).toFixed(4))).join(",");
        }

        const valid_ctrl = ctrl_counts.filter((x) => !Number.isNaN(x));
        let ctrl_sum = 0.0;
        for (let r = 0; r < valid_ctrl.length; r++) ctrl_sum += valid_ctrl[r];
        const ctrl_mean = valid_ctrl.length > 0 ? ctrl_sum / valid_ctrl.length : 0.0;

        // Exact saintPy FoldChange (spec_sum / ctrl_mean if ctrl_mean > 0 else NaN)
        const fold_change = ctrl_mean > 0 ? spec_sum / ctrl_mean : NaN;

        // Regularized FoldChange & log2FC for seamless Volcano visualization (even when ctrl_mean == 0 or INT log-space)
        let empiricalFC;
        let log2FC;
        if (quant_type === "spc") {
          const denom = ctrl_mean > 0 ? ctrl_mean : Math.max(model.eta[i], 0.1);
          empiricalFC = avg_spec / denom;
          log2FC = Math.log2(Math.max(avg_spec, 0.05) / denom);
        } else {
          // In INT mode, spec and ctrl_counts are natural-log intensities
          const baseLog = valid_ctrl.length > 0 ? ctrl_mean : model.eta[i];
          const diffLn = avg_spec - baseLog;
          log2FC = diffLn / Math.LN2;
          empiricalFC = Math.pow(2, log2FC);
        }

        rows.push({
          Bait: ubait[bait_idx],
          Prey: prey_id,
          PreyGene: gene_name,
          Spec: spec_str,
          SpecSum: spec_sum,
          AvgSpec: avg_spec,
          NumReplicates: num_reps,
          ctrlCounts: ctrl_str,
          AvgP: avg_p,
          MaxP: max_p,
          TopoAvgP: avg_p,
          TopoMaxP: max_p,
          SaintScore: avg_p,
          logOddsScore: log_odds,
          FoldChange: fold_change,
          BFDR: NaN,
          boosted_by: "",
          // Computed fields for UI & plotting
          computedBFDR: NaN,
          empiricalFC,
          log2FC,
          _zeroCtrl: ctrl_mean <= 0,
          _preyLength: prey_length,
          _ctrlMean: ctrl_mean,
          _eta: model.eta[i],
          _d: model.d[i],
          _replicatesRaw: spec,
          _replicateNames: rep_ip_names,
          _controlsRaw: ctrl_counts,
          _controlNames: ctrl_columns,
          _zBits: Array.from(model.Z[i][bait_idx]),
          _repScores:
            model._rep_scores_mat && model._rep_scores_mat[i] && model._rep_scores_mat[i][bait_idx]
              ? model._rep_scores_mat[i][bait_idx]
              : [],
        });
      }
    }

    computeBFDR(rows);
    return rows;
  }

  /**
   * Main pipeline entry point matching saintpy/main.py run_pipeline()
   */
  function run_pipeline({
    samplesRows,
    interRows,
    preyRows = null,
    opts = { f: 0.0, R: 100, L: 100 },
    quant_type = "spc",
    selfInteractionMode = "saintpy",
    baitAliases = {},
    skipExcludedSamples = true,
  }) {
    const startTime = Date.now();
    const T = build_inputs(samplesRows, interRows, preyRows, quant_type, "quant", {
      selfInteractionMode,
      baitAliases,
      skipExcludedSamples,
    });

    const warnings = [];
    if (opts.L > T.nctrl) {
      warnings.push(
        `Parameter L (${opts.L}) is larger than the number of control IPs (${T.nctrl}); using L = ${T.nctrl}.`
      );
    }

    const p2p_mapping = Array.from({ length: T.nprey }, () => []);
    let model;
    if (quant_type === "int") {
      model = statModel_int(
        p2p_mapping,
        T.ubait,
        T.test_mat_DATA,
        T.ctrl_mat_DATA,
        T.ip_idx_to_bait_no,
        T.nprey,
        T.nbait,
        opts
      );
    } else {
      model = statModel_spc(
        p2p_mapping,
        T.ubait,
        T.test_mat_DATA,
        T.ctrl_mat_DATA,
        T.ip_idx_to_bait_no,
        T.nprey,
        T.nbait,
        opts
      );
    }
    model.PDATA = T.PDATA;
    model.ubait = T.ubait;

    const icmTrace = [];
    for (let iter = 0; iter < 15; iter++) {
      const old_ll = model.llikelihood();
      if (quant_type === "int") {
        model.icm_Z();
        model.wrt_MRF();
      } else {
        icm_Z_spc(model);
        wrt_MRF_gamma_0_spc(model);
      }
      const new_ll = model.llikelihood();
      const diff = new_ll - old_ll;
      icmTrace.push({
        iteration: iter + 1,
        old_ll,
        new_ll,
        diff,
        beta1: model.beta1,
      });
      if (diff >= 0 && diff < 1e-3) {
        break;
      }
    }

    let avgp, maxp, min_logodds;
    if (quant_type === "int") {
      [avgp, maxp, min_logodds] = model.calculateScore();
    } else {
      [avgp, maxp, min_logodds] = calculateScore_spc(model);
    }

    const test_mat_mask = _build_test_mask(T.nprey, T.nbait, T.test_mat_DATA, T.ip_idx_to_bait_no);
    const rows = build_output_rows(
      model,
      avgp,
      maxp,
      min_logodds,
      T.ip_idx_to_bait_no,
      T.test_mat_DATA,
      T.ctrl_mat_DATA,
      test_mat_mask,
      quant_type,
      T.test_columns,
      T.ctrl_columns
    );

    const elapsedMs = Date.now() - startTime;
    return {
      rows,
      warnings,
      summary: {
        nbait: T.nbait,
        nprey: T.nprey,
        nexpr: T.nexpr,
        nctrl: T.nctrl,
        totalPairs: rows.length,
        quant_type,
        beta1: model.beta1,
        gamma: model.gamma,
        icmIterations: icmTrace.length,
        icmTrace,
        elapsedMs,
        test_columns: T.test_columns,
        ctrl_columns: T.ctrl_columns,
        ubait: T.ubait,
        maskedSelfInteractions: T.maskedSelfInteractions,
      },
    };
  }

  // ============================================================================
  // 5. FRAGPIPE PARSER & SMART BAIT-PREY DESIGNER ENGINE
  // ============================================================================

  const CONTROL_KEYWORDS = [
    "ctrl",
    "control",
    "igg",
    "gfp",
    "empty",
    "mock",
    "vec",
    "vector",
    "bead",
    "beads",
    "bira",
    "nls",
    "flag_only",
    "flag_empty",
    "wt_ctrl",
    "neg",
    "negative",
    "untransfected",
    "parental",
  ];

  function inferBaitAndControlFromSampleName(sampleName, annotationMap = null) {
    const clean = String(sampleName || "").trim();
    if (annotationMap && annotationMap.has(clean)) {
      const ann = annotationMap.get(clean);
      const cond = ann.condition || ann.bait_id || clean;
      const isCtrl =
        ann.is_control !== undefined
          ? ann.is_control
          : CONTROL_KEYWORDS.some((kw) => cond.toLowerCase().includes(kw) || clean.toLowerCase().includes(kw));
      return {
        ip_id: clean,
        bait_id: cond,
        is_control: Boolean(isCtrl),
        replicate: ann.bioreplicate || "",
      };
    }

    const lower = clean.toLowerCase();
    const is_control = CONTROL_KEYWORDS.some((kw) => {
      const re = new RegExp(`(^|[_\\-\\.\\s])${kw}([0-9_\\-\\.\\s]|$)`, "i");
      return re.test(lower) || lower.startsWith(kw);
    });

    // Strip trailing replicate patterns like _1, _2, _rep1, -R1, _IP1, .1
    let bait_id = clean
      .replace(/[_\-\.](?:rep|r|ip|bio|run)?[0-9]+$/i, "")
      .replace(/\s+\([0-9]+\)$/, "")
      .trim();

    if (!bait_id) bait_id = clean;
    if (is_control && /^ctrl[0-9]*$/i.test(clean)) {
      bait_id = "CTRL";
    }

    return {
      ip_id: clean,
      bait_id,
      is_control,
      replicate: "",
    };
  }

  /**
   * Parse FragPipe experiment_annotation.tsv or fp-manifest
   */
  function parseAnnotationTable(text) {
    const { headers, rows } = parseDelimitedText(text);
    const map = new Map();
    for (const r of rows) {
      const sample = (r.sample || r.sample_name || r.experiment || r.ip_id || "").trim();
      if (!sample) continue;
      const condition = (r.condition || r.bait || r.bait_id || sample).trim();
      const bioreplicate = (r.bioreplicate || r.replicate || "").trim();
      let is_control = undefined;
      if (r.is_control !== undefined && r.is_control !== "") {
        try {
          is_control = coerceBool(r.is_control);
        } catch (_) {}
      }
      map.set(sample, { sample, condition, bioreplicate, is_control });
    }
    return map;
  }

  /**
   * Parse FragPipe combined_protein.tsv, protein.tsv, or reprint.spc.tsv / reprint.int.tsv
   */
  function parseFragPipeTable(text, annotationText = null) {
    const { headers, rows } = parseDelimitedText(text);
    if (headers.length === 0 || rows.length === 0) {
      throw new Error("The uploaded FragPipe file is empty or could not be parsed.");
    }

    const annotationMap = annotationText ? parseAnnotationTable(annotationText) : null;

    // Check if this is a REPRINT matrix format (PROTID, GENEID, LENGTH, <sample1>, <sample2>, ...)
    const upperHeaders = headers.map((h) => h.toUpperCase());
    const isReprint =
      upperHeaders.includes("PROTID") ||
      (upperHeaders.includes("GENEID") && upperHeaders.includes("LENGTH"));

    // Ordered longest-suffix-first so " Unique Spectral Count" and " Total Spectral Count"
    // always match before " Spectral Count" (and similarly for Intensity variants).
    const METRIC_SUFFIXES = [
      { suffix: " Unique Spectral Count", label: "Unique Spectral Count", quant_type: "spc" },
      { suffix: " Total Spectral Count", label: "Total Spectral Count", quant_type: "spc" },
      { suffix: " Razor Spectral Count", label: "Razor Spectral Count", quant_type: "spc" },
      { suffix: " Spectral Count", label: "Spectral Count", quant_type: "spc" },
      { suffix: " MaxLFQ Intensity", label: "MaxLFQ Intensity", quant_type: "int" },
      { suffix: " Unique Intensity", label: "Unique Intensity", quant_type: "int" },
      { suffix: " Total Intensity", label: "Total Intensity", quant_type: "int" },
      { suffix: " Razor Intensity", label: "Razor Intensity", quant_type: "int" },
      { suffix: " Intensity", label: "Intensity", quant_type: "int" },
    ];

    const EXCLUDED_COMBINED_PREFIXES = ["combined ", "total ", "razor "];

    const metricsFound = {}; // label -> { quant_type, sampleToCol: Map<sample, colName> }

    if (!isReprint) {
      for (const h of headers) {
        const hLow = h.toLowerCase();
        if (EXCLUDED_COMBINED_PREFIXES.some((pref) => hLow.startsWith(pref))) {
          continue;
        }
        // Check longest suffixes first so " Unique Spectral Count" matches before " Spectral Count"
        for (const m of METRIC_SUFFIXES) {
          if (h.endsWith(m.suffix)) {
            const sampleName = h.slice(0, h.length - m.suffix.length).trim();
            if (!sampleName) continue;
            if (!metricsFound[m.label]) {
              metricsFound[m.label] = {
                label: m.label,
                quant_type: m.quant_type,
                sampleToCol: new Map(),
              };
            }
            metricsFound[m.label].sampleToCol.set(sampleName, h);
            break;
          }
        }
      }
    }

    // If no suffixed columns found, check if it's a matrix format (like reprint.spc.tsv or custom wide matrix)
    const availableMetrics = Object.keys(metricsFound);
    if (availableMetrics.length === 0) {
      const metaCols = new Set([
        "_index",
        "protein",
        "protein id",
        "entry name",
        "gene",
        "gene_name",
        "geneid",
        "protid",
        "prey_id",
        "length",
        "prey_length",
        "organism",
        "description",
        "protein existence",
        "coverage",
        "protein probability",
        "top peptide probability",
        "combined total peptides",
        "combined spectral count",
        "combined unique spectral count",
        "combined total spectral count",
        "indistinguishable proteins",
        "num_reps",
      ]);
      const sampleCols = headers.filter((h) => !metaCols.has(h.toLowerCase()));
      if (sampleCols.length === 0) {
        throw new Error(
          "Could not detect sample quantification columns in the file. Expected '<Sample> Spectral Count', '<Sample> Intensity', or sample matrix columns."
        );
      }
      // Infer whether matrix values are integers (SPC) or floats (INT)
      let hasNonInteger = false;
      for (let i = 0; i < Math.min(rows.length, 40); i++) {
        for (const col of sampleCols) {
          const v = Number(rows[i][col]);
          if (!Number.isNaN(v) && v > 0 && !Number.isInteger(v)) {
            hasNonInteger = true;
            break;
          }
        }
      }
      const defaultLabel = hasNonInteger ? "Matrix Intensity" : "Matrix Spectral Count";
      metricsFound[defaultLabel] = {
        label: defaultLabel,
        quant_type: hasNonInteger ? "int" : "spc",
        sampleToCol: new Map(sampleCols.map((c) => [c, c])),
      };
    }

    const metricList = Object.values(metricsFound);
    const primaryMetric =
      metricList.find((m) => m.label === "Spectral Count") ||
      metricList.find((m) => m.label === "Unique Spectral Count") ||
      metricList[0];

    // Collect union of all sample names in deterministic order
    const allSampleNames = Array.from(primaryMetric.sampleToCol.keys());
    const samplesMeta = allSampleNames.map((s) => inferBaitAndControlFromSampleName(s, annotationMap));

    // Inspect available ID columns
    const idColCandidates = [];
    for (const candidate of ["Gene", "Protein ID", "Protein", "Entry Name", "PROTID", "GENEID", "prey_id"]) {
      const actual = headers.find((h) => h.toLowerCase() === candidate.toLowerCase());
      if (actual) idColCandidates.push(actual);
    }
    const defaultPreyIdCol =
      idColCandidates.find((c) => c.toLowerCase() === "gene") ||
      idColCandidates.find((c) => c.toLowerCase() === "protein id") ||
      idColCandidates[0] ||
      headers[0];
    const defaultGeneCol =
      idColCandidates.find((c) => c.toLowerCase().includes("gene")) || defaultPreyIdCol;
    const lengthCol =
      headers.find((h) => ["length", "prey_length"].includes(h.toLowerCase())) || null;
    const probCol =
      headers.find((h) => h.toLowerCase() === "protein probability") || null;
    const peptidesCol =
      headers.find((h) => h.toLowerCase() === "combined total peptides" || h.toLowerCase() === "total peptides") ||
      null;
    const descCol =
      headers.find((h) => h.toLowerCase() === "description") || null;

    return {
      headers,
      rows,
      metricsFound,
      metricLabels: metricList.map((m) => m.label),
      defaultMetricLabel: primaryMetric.label,
      samplesMeta,
      idColCandidates,
      defaultPreyIdCol,
      defaultGeneCol,
      lengthCol,
      probCol,
      peptidesCol,
      descCol,
    };
  }

  /**
   * Given parsed FragPipe table and user's experimental design & prey filter settings,
   * generate { samplesRows, preyRows, interRows, filterStats, baitSelfPulldown }
   */
  function extractSaintInputsFromFragPipe(parsedFP, config) {
    const {
      metricLabel = parsedFP.defaultMetricLabel,
      preyIdCol = parsedFP.defaultPreyIdCol,
      geneCol = parsedFP.defaultGeneCol,
      samplesMeta = parsedFP.samplesMeta,
      excludeContaminants = true,
      excludeDecoys = true,
      minProteinProb = 0.0,
      minPeptides = 0,
      minMaxQuant = 0,
      excludedPreyIds = new Set(),
      baitAliases = {},
    } = config;

    const metricInfo = parsedFP.metricsFound[metricLabel];
    if (!metricInfo) {
      throw new Error(`Metric '${metricLabel}' not found in FragPipe table.`);
    }

    const activeSamples = samplesMeta.filter((s) => !s.excluded);
    const samplesRows = activeSamples.map((s) => ({
      ip_id: s.ip_id,
      bait_id: s.bait_id,
      is_control: Boolean(s.is_control),
    }));

    let totalRawRows = parsedFP.rows.length;
    let filteredContam = 0;
    let filteredDecoy = 0;
    let filteredProb = 0;
    let filteredPeptides = 0;
    let filteredZeroOrLow = 0;
    let filteredManual = 0;

    const preyRows = [];
    const interRows = [];
    const preyPreviewList = [];

    for (const row of parsedFP.rows) {
      const proteinFull = String(row["Protein"] || row[preyIdCol] || "").trim();
      const rawPreyId = String(row[preyIdCol] || row["Protein ID"] || proteinFull).trim();
      const rawGene = String(row[geneCol] || rawPreyId).trim() || rawPreyId;

      if (!rawPreyId) continue;

      // Contaminant check (FragPipe marks contaminants with contam_ or cRAP or CON__)
      const isContam =
        /^contam_|^con__|crap/i.test(proteinFull) || /^contam_|^con__/i.test(rawPreyId);
      if (excludeContaminants && isContam) {
        filteredContam++;
        continue;
      }

      // Decoy check (rev_ or REV__)
      const isDecoy = /^rev_|^rev__/i.test(proteinFull) || /^rev_|^rev__/i.test(rawPreyId);
      if (excludeDecoys && isDecoy) {
        filteredDecoy++;
        continue;
      }

      // Protein Probability check
      if (parsedFP.probCol && minProteinProb > 0) {
        const prob = Number(row[parsedFP.probCol]);
        if (!Number.isNaN(prob) && prob < minProteinProb) {
          filteredProb++;
          continue;
        }
      }

      // Minimum Peptides check
      if (parsedFP.peptidesCol && minPeptides > 0) {
        const pep = Number(row[parsedFP.peptidesCol]);
        if (!Number.isNaN(pep) && pep < minPeptides) {
          filteredPeptides++;
          continue;
        }
      }

      if (excludedPreyIds && excludedPreyIds.has(rawPreyId)) {
        filteredManual++;
        continue;
      }

      // Extract sample quantifications
      let maxQuantAcrossActive = 0;
      let testSum = 0;
      let ctrlSum = 0;
      const sampleQuants = [];

      for (const s of activeSamples) {
        const colName = metricInfo.sampleToCol.get(s.ip_id);
        if (!colName) continue;
        const rawVal = row[colName];
        const val = rawVal === "" || rawVal === undefined ? 0 : Number(rawVal);
        const cleanVal = Number.isNaN(val) ? 0 : val;
        if (cleanVal > maxQuantAcrossActive) maxQuantAcrossActive = cleanVal;
        if (s.is_control) ctrlSum += cleanVal;
        else testSum += cleanVal;
        sampleQuants.push({ ip_id: s.ip_id, quant: cleanVal, is_control: s.is_control });
      }

      if (maxQuantAcrossActive <= 0 || maxQuantAcrossActive < minMaxQuant) {
        filteredZeroOrLow++;
        continue;
      }

      const lenVal =
        parsedFP.lengthCol && row[parsedFP.lengthCol] !== ""
          ? Number(row[parsedFP.lengthCol])
          : "";

      preyRows.push({
        prey_id: rawPreyId,
        gene_name: rawGene,
        length: lenVal,
      });

      preyPreviewList.push({
        prey_id: rawPreyId,
        gene_name: rawGene,
        protein: proteinFull,
        length: lenVal,
        description: parsedFP.descCol ? row[parsedFP.descCol] : "",
        maxQuant: maxQuantAcrossActive,
        testSum,
        ctrlSum,
      });

      for (const sq of sampleQuants) {
        if (metricInfo.quant_type === "spc") {
          // In SPC mode, only emit non-zero interactions (or all, since saintPy pivots with fill_value=0)
          if (sq.quant > 0) {
            interRows.push({
              ip_id: sq.ip_id,
              prey_id: rawPreyId,
              quant: sq.quant,
            });
          }
        } else {
          // In INT mode, emit positive intensities (<=0 becomes NaN in saintPy)
          if (sq.quant > 0) {
            interRows.push({
              ip_id: sq.ip_id,
              prey_id: rawPreyId,
              quant: sq.quant,
            });
          }
        }
      }
    }

    // Detect Bait Self-Pulldown QC for each unique Test Bait
    const baitSelfPulldown = detectBaitSelfPreys(samplesRows, interRows, preyRows, baitAliases);

    return {
      quant_type: metricInfo.quant_type,
      samplesRows,
      preyRows,
      interRows,
      preyPreviewList,
      baitSelfPulldown,
      filterStats: {
        totalRawRows,
        retainedPreys: preyRows.length,
        retainedInteractions: interRows.length,
        filteredContam,
        filteredDecoy,
        filteredProb,
        filteredPeptides,
        filteredZeroOrLow,
        filteredManual,
      },
    };
  }

  /**
   * Inspect whether each test bait_id matches a prey_id or gene_name in the dataset,
   * and compute its total abundance in its own replicates vs controls.
   */
  function detectBaitSelfPreys(samplesRows, interRows, preyRows, baitAliases = {}) {
    const activeBaits = [];
    const baitToIps = new Map();
    const ctrlIps = new Set();

    for (const s of samplesRows) {
      if (s.excluded) continue;
      const isCtrl = coerceBool(s.is_control ?? s.tc ?? false);
      if (isCtrl) {
        ctrlIps.add(s.ip_id);
      } else {
        if (!baitToIps.has(s.bait_id)) {
          baitToIps.set(s.bait_id, new Set());
          activeBaits.push(s.bait_id);
        }
        baitToIps.get(s.bait_id).add(s.ip_id);
      }
    }

    const preyLookup = new Map(); // lowercase prey_id or gene_name -> { prey_id, gene_name }
    if (Array.isArray(preyRows)) {
      for (const p of preyRows) {
        const pid = String(p.prey_id || "").trim();
        const gene = String(p.gene_name || p.prey_gene || pid).trim();
        if (pid) preyLookup.set(pid.toLowerCase(), { prey_id: pid, gene_name: gene });
        if (gene) preyLookup.set(gene.toLowerCase(), { prey_id: pid, gene_name: gene });
      }
    }

    // Also scan interRows in case preyRows was omitted
    for (const r of interRows) {
      const pid = String(r.prey_id || "").trim();
      if (pid && !preyLookup.has(pid.toLowerCase())) {
        preyLookup.set(pid.toLowerCase(), { prey_id: pid, gene_name: pid });
      }
    }

    // Index interactions by prey_id
    const preyQuants = new Map();
    for (const r of interRows) {
      const pid = String(r.prey_id || "").trim();
      if (!preyQuants.has(pid)) preyQuants.set(pid, []);
      preyQuants.get(pid).push({ ip_id: r.ip_id, quant: Number(r.quant) || 0 });
    }

    const report = [];
    for (const bait of activeBaits) {
      const alias = baitAliases[bait] || bait;
      const match =
        preyLookup.get(String(alias).toLowerCase()) ||
        preyLookup.get(String(bait).toLowerCase()) ||
        null;
      const ownIps = baitToIps.get(bait);
      let ownSum = 0;
      let ctrlSum = 0;
      if (match) {
        const entries = preyQuants.get(match.prey_id) || [];
        for (const e of entries) {
          if (ownIps.has(e.ip_id)) ownSum += e.quant;
          else if (ctrlIps.has(e.ip_id)) ctrlSum += e.quant;
        }
      }
      report.push({
        bait_id: bait,
        numReplicates: ownIps.size,
        replicates: Array.from(ownIps),
        matchedPreyId: match ? match.prey_id : null,
        matchedGene: match ? match.gene_name : null,
        exactGeneMatch: match ? match.gene_name === bait : false,
        ownAbundanceSum: ownSum,
        ctrlAbundanceSum: ctrlSum,
      });
    }
    return report;
  }

  /**
   * Format output rows to standard saintPy CSV string
   */
  function formatOutputCSV(rows, includeComputedBFDR = true) {
    const cols = [
      "Bait",
      "Prey",
      "PreyGene",
      "Spec",
      "SpecSum",
      "AvgSpec",
      "NumReplicates",
      "ctrlCounts",
      "AvgP",
      "MaxP",
      "TopoAvgP",
      "TopoMaxP",
      "SaintScore",
      "logOddsScore",
      "FoldChange",
      "BFDR",
      "boosted_by",
    ];
    const lines = [cols.join(",")];
    for (const r of rows) {
      const vals = cols.map((c) => {
        let v = r[c];
        if (c === "BFDR" && includeComputedBFDR) {
          v = r.computedBFDR;
        }
        if (v === undefined || v === null || (typeof v === "number" && Number.isNaN(v))) {
          return "";
        }
        const s = String(v);
        if (s.includes(",") || s.includes('"') || s.includes("\n")) {
          return `"${s.replace(/"/g, '""')}"`;
        }
        return s;
      });
      lines.push(vals.join(","));
    }
    return lines.join("\n");
  }

  return {
    run_pipeline,
    build_inputs,
    parseDelimitedText,
    parseFragPipeTable,
    parseAnnotationTable,
    inferBaitAndControlFromSampleName,
    extractSaintInputsFromFragPipe,
    detectBaitSelfPreys,
    formatOutputCSV,
    GP_log_pmf1,
    _quant_log_pdf,
  };
});
