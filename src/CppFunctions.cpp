#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// @description Computes the distance matrix.
// @param coords \eqn{n \times 2} coordinates matrix.
// [[Rcpp::export]]
arma::mat dist1(arma::mat coords) {
  int n = coords.n_rows;
  arma::mat D(n, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      D(i, j) = arma::norm(coords.row(i) - coords.row(j));
      D(j, i) = D(i, j);
    }
  }
  return D;
}

// @description Computes the squared distance matrix.
// @param coords \eqn{n \times m} coordinates matrix.
// [[Rcpp::export]]
arma::mat dist2(arma::mat coords) {
  int n = coords.n_rows;
  int m = coords.n_cols;
  arma::mat D(n, n, arma::fill::zeros);
  arma::rowvec v(m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      v = coords.row(i) - coords.row(j);
      D(i, j) = arma::dot(v, v);
      D(j, i) = D(i, j);
    }
  }
  return D;
}

// Online Update of Mean and Variance
//
// Update mu and Sigma for tunning the Metropolis step
//
// @param x Vector New value
// @param mu Vector Old mean
// @param Sigma Matrix Old variance
// @param n Number of data with x
// @return List with updated mu and Sigma
// [[Rcpp::export]]
Rcpp::List MuSigmaUpdate(arma::vec x, arma::vec mu, arma::mat Sigma, int n) {
  arma::vec u = (x - mu) / n;
  mu += u;
  Sigma = (n - 2) * Sigma / (n - 1) + u * u.t() * n;
  return Rcpp::List::create(Rcpp::Named("mu")=mu, Rcpp::Named("Sigma")=Sigma);
}

// [[Rcpp::export]]
Rcpp::List MuSigma(arma::mat x) {
  arma::rowvec mu = arma::mean(x, 0);
  arma::mat Sigma = arma::cov(x);
  return Rcpp::List::create(Rcpp::Named("mu")=mu.t(), Rcpp::Named("Sigma")=Sigma);
}

// [[Rcpp::export]]
Rcpp::List MuSigmaUpdate2(
    arma::vec x_new, arma::vec x_old,
    arma::vec mu_old, arma::mat Sigma_old,
    int fn_new, int fn_old,
    int n_new, int d) {

  arma::vec mu_new(d);
  arma::mat Sigma_new(d, d);
  int diff_n_fn = n_new - fn_new;

  if (fn_new == fn_old) {
    mu_new = (diff_n_fn * mu_old + x_new) / (diff_n_fn + 1);
    Sigma_new =
      mu_old * mu_old.t() + (
        (diff_n_fn - 1) * Sigma_old +
        x_new * x_new.t() -
        (diff_n_fn + 1) * mu_new * mu_new.t()
      ) / diff_n_fn;
  } else if (fn_new == fn_old + 1) {
    mu_new = mu_old + (x_new - x_old) / (diff_n_fn + 1);
    Sigma_new = Sigma_old + (
      x_new * x_new.t() -
      x_old * x_old.t() +
      (diff_n_fn + 1) * (mu_old * mu_old.t() -
                         mu_new * mu_new.t())
    ) / diff_n_fn;
  }

  return Rcpp::List::create(Rcpp::Named("mu")=mu_new, Rcpp::Named("Sigma")=Sigma_new);
}

// [[Rcpp::export]]
Rcpp::List MuSigmaUpdate3(
    arma::vec x_new, arma::vec x_old,
    arma::vec mu_old, arma::mat Sigma_old,
    int fn_new, int fn_old,
    int n_new, int d) {

  int nu0 = 2 * d;

  arma::vec mu_new(d);
  arma::mat Sigma_new(d, d);
  int diff_n_fn = n_new - fn_new;

  if (n_new == 1) {
    mu_new = (mu_old + x_new) / 2;
    Sigma_new = (
      mu_old * mu_old.t() +
      x_new * x_new.t() -
      2 * mu_new * mu_new.t() +
      (nu0 + d + 1) * Sigma_old
    ) / (nu0 + d + 3);
  } else if (fn_new == fn_old) {
    mu_new = (diff_n_fn * mu_old + x_new) / (diff_n_fn + 1);
    Sigma_new = (
      (diff_n_fn + nu0 + d + 1) * Sigma_old +
      x_new * x_new.t() +
      diff_n_fn * mu_old * mu_old.t() -
      (diff_n_fn + 1) * mu_new * mu_new.t()
    ) / (diff_n_fn + nu0 + d + 2);
  } else if (fn_new == fn_old + 1) {
    mu_new = mu_old + (x_new - x_old) / (diff_n_fn + 1);
    Sigma_new = Sigma_old + (
      x_new * x_new.t() -
      x_old * x_old.t() +
      (diff_n_fn + 1) * (mu_old * mu_old.t() -
                         mu_new * mu_new.t())
    ) / (diff_n_fn + nu0 + d + 2);
  }

  return Rcpp::List::create(Rcpp::Named("mu")=mu_new, Rcpp::Named("Sigma")=Sigma_new);
}

// @description Random generation for the truncated normal distribution with
//   mean equal to mu, standard deviation equal to 1 and lower and upper
//   limits (0, Inf).
// @param N number of observations.
// @param mu vector of means.
// [[Rcpp::export]]
arma::vec rtnorm(const int N, arma::vec mu) {

  double a, alpha, z, g;
  arma::vec x(N);

  for (int i = 0; i < N; ++i) {
    a = - mu(i);
    if (a > 0.67448975) { // Robert 1995
      alpha = (a + sqrt(a * a + 4)) / 2;
      do {
        z = R::rexp(1 / alpha) + a;
        g = exp(-pow(z - alpha, 2) / 2);
      } while (R::runif(0, 1) > g);
    } else { // inverse method
      z = -R::qnorm(R::pnorm(-a, 0, 1, 1, 0) * R::runif(0, 1), 0, 1, 1, 0);
    }

    x(i) = z + mu(i);
  }

  return x;
}

// As rtnorm but for numbers and with sigma (sd). Limits (0, Inf).
// [[Rcpp::export]]
double rtnorm_numeric(double mu, double sigma) {

  double a, alpha, z, g;

  a = - mu / sigma;
  if (a > 0.67448975) { // Robert 1995
    alpha = (a + sqrt(a * a + 4)) / 2;
    do {
      z = R::rexp(1 / alpha) + a;
      g = exp(-pow(z - alpha, 2) / 2);
    } while (R::runif(0, 1) > g);
  } else { // inverse method
    z = -R::qnorm(R::pnorm(-a, 0, 1, 1, 0) * R::runif(0, 1), 0, 1, 1, 0);
  }

  return z * sigma + mu;
}

// @description Not vectorized truncated normal in (a, b).
// [[Rcpp::export]]
double rtnorm1(double mu, double sigma, double a, double b) {

  double pAlpha = R::pnorm(a, mu, sigma, 1, 0);
  double pBeta  = R::pnorm(b, mu, sigma, 1, 0);
  double x = R::qnorm(pAlpha + R::runif(0, 1) * (pBeta - pAlpha), mu, sigma, 1, 0);

  return x;
}

// [[Rcpp::export]]
double dtnorm1(double x, double mu, double sigma, double a, double b) { // log
  return R::dnorm(x, mu, sigma, 1) - std::log(R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0));
}

// [[Rcpp::export]]
double dnorm1(double x, double mu, double prec) { // log
  x = x - mu;
  return - x * x * prec / 2;
}

// [[Rcpp::export]]
double dA(arma::vec la, arma::vec a,
          arma::vec mu, arma::mat Sigma_inv,
          arma::vec mu_likelihood, arma::mat Sigma_inv_likelihood) { // log
  double value;
  arma::vec process = la - mu;
  value = arma::as_scalar(process.t() * Sigma_inv * process);
  process = a - mu_likelihood;
  value += arma::as_scalar(process.t() * Sigma_inv_likelihood * process);
  return - value / 2;
}



// Univariate model (non-spatial coregionalization)
// [[Rcpp::export]]
arma::mat uglm(
    const int N, const int k, const arma::vec Nweak,
    arma::mat X, const arma::mat X0, const arma::mat X1,
    arma::vec beta, arma::vec Wtls,
    double a11,
    arma::vec decay, const int M,
    const arma::cube dist, const int decayPrior,
    const int n, const int T, const int L,
    const arma::uvec site, const arma::uvec year, const arma::uvec day,
    const arma::uvec weak, const arma::uvec Weak,
    const arma::vec noNA, const arma::uvec indLag,
    const arma::vec prob, arma::vec abLim,
    const arma::uvec weakl0, const arma::vec probl0,
    const arma::mat V,
    const double na, const double nb,
    const double diag_b,
    const arma::vec da, const arma::vec db,
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {

  const int TL = T * L;
  arma::mat Xb = X * beta;

  arma::vec Z(N, arma::fill::zeros); // the first sampled (no need to specify)
  arma::vec auxZ = Z - Xb;
  auxZ -= a11 * Wtls;

  arma::mat In(n, n, arma::fill::eye);
  arma::mat IM(M, M, arma::fill::eye);
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaW(n, n);

  // GPs
  double delta, chi;
  arma::vec vn(n);

  arma::uvec indWUvec(n);
  arma::umat indWUmat(n, TL);
  for (int l = 0; l < L; ++l) {
    for (int t = 0; t < T; ++t) {
      indWUmat.col(l * T + t) = arma::find((year == t) && (day == l));
    }
  }

  // decay (starter)
  arma::vec accept(M, arma::fill::zeros);
  int total = 0;
  arma::vec r(M);
  arma::vec sd(M, arma::fill::ones);
  arma::vec lsd(M, arma::fill::zeros); //log(sd);

  arma::mat Rinv(n, n);
  double Rlogdet;

  arma::mat dd_aux(n, n);
  arma::mat dd(n, n, arma::fill::zeros);
  for (int m = 0; m < M; ++m) {
    dd -= decay(m) * dist.slice(m);
  }
  Rinv = arma::inv_sympd(exp(dd));
  Rlogdet = arma::log_det_sympd(Rinv);

  arma::vec decay_aux(M);
  arma::vec ldecay_aux(M, arma::fill::zeros);
  arma::vec ldecay = log(decay);
  arma::mat Rinv_aux(n, n);
  double Rlogdet_aux;

  arma::mat WtRW(1, 1), WtRW_aux(1, 1);
  double A = 0;

  // weak
  arma::uvec ind;

  // LOOP
  for (int b = 1 - nBurnin; b <= nSims; ++b) {

    // check
    if (b % nReport == 0) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
      if ((b < 1) && (total != 0)) {
        Rcpp::Rcout << "Acceptance rate : " << accept.t() / total << "\n";
      } else if (b < 1) {
        Rcpp::Rcout << "Acceptance rate : " << 0 << "\n";
      } else {
        Rcpp::Rcout << "Acceptance rate : " <<  accept.t() / b << "\n";
      }
    }

    // weak records
    auxZ += Xb;
    for (int i = 0; i < Nweak(0); i++) { // Y and lag
      if (R::runif(0, 1) < prob(weak(i))) {
        abLim(weak(i)) = 1;
        if (noNA(i)) {
          ind = Weak(i);
          X.submat(ind, indLag) = X1.submat(ind, indLag);
        }
      } else {
        abLim(weak(i)) = -1;
        if (noNA(i)) {
          ind = Weak(i);
          X.submat(ind, indLag)  = X0.submat(ind, indLag);
        }
      }
    }
    if (Nweak(1) > 0) {
      for (int i = 0; i < Nweak(1); i++) { // lag l = 0
        ind = weakl0(i);
        if (R::runif(0, 1) < probl0(i)) {
          X.submat(ind, indLag) = X1.submat(ind, indLag);
        } else {
          X.submat(ind, indLag) = X0.submat(ind, indLag);
        }
      }
    }
    Xb = X * beta;
    auxZ -= Xb;

    // Z
    auxZ -= Z;
    Z = abLim % rtnorm(N, -abLim % auxZ);
    auxZ += Z;

    // beta
    auxZ += Xb;
    OmegaBeta = arma::inv_sympd(X.t() * X + V);
    beta = OmegaBeta * (X.t() * auxZ) + // prior 0k
      arma::chol(OmegaBeta, "lower") * arma::randn(k);
    Xb = X * beta;
    auxZ -= Xb;

    // Wtl(s)
    auxZ += a11 * Wtls;
    OmegaW = arma::inv_sympd(a11 * a11 * In + Rinv);
    for (int l = 0; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        Wtls(indWUvec) = a11 * auxZ(indWUvec);
        Wtls(indWUvec) = OmegaW * Wtls(indWUvec) + arma::chol(OmegaW, "lower") * arma::randn(n);
      }
    }
    auxZ -= a11 * Wtls;

    // a11
    auxZ += a11 * Wtls;
    delta = 1 / (arma::dot(Wtls, Wtls) + diag_b);
    chi   = arma::dot(Wtls, auxZ); // + 0 * diag_b;
    a11   = rtnorm_numeric(delta * chi, sqrt(delta));
    auxZ -= a11 * Wtls;

    // decay
    if (decayPrior != 1) {
      for (int m = 0; m < M; ++m) {
        if (decayPrior == 2) {
          decay_aux(m) = rtnorm1(decay(m), sd(m), da(m), db(m));
          dd_aux       = dd + (decay(m) - decay_aux(m)) * dist.slice(m);
          Rinv_aux     = arma::inv_sympd(exp(dd_aux));
          Rlogdet_aux  = arma::log_det_sympd(Rinv_aux);
          WtRW_aux.zeros();
          WtRW.zeros();
          for (int l = 0; l < L; ++l) {
            for (int t = 0; t < T; ++t) {
              indWUvec = indWUmat.col(l * T + t);
              vn = Wtls(indWUvec);
              WtRW_aux += vn.t() * Rinv_aux * vn;
              WtRW += vn.t() * Rinv * vn;
            }
          }
          A =
            ((TL * Rlogdet_aux - arma::as_scalar(WtRW_aux)) -
             (TL * Rlogdet     - arma::as_scalar(WtRW)    )) / 2 +
            dtnorm1(decay(m)    , decay_aux(m), sd(m), da(m), db(m)) -
            dtnorm1(decay_aux(m), decay(m)    , sd(m), da(m), db(m));
        } else if (decayPrior == 3) {
          ldecay_aux(m) = R::rnorm(ldecay(m), sd(m));
          decay_aux(m)  = exp(ldecay_aux(m));
          dd_aux       = dd + (decay(m) - decay_aux(m)) * dist.slice(m);
          Rinv_aux     = arma::inv_sympd(exp(dd_aux));
          Rlogdet_aux  = arma::log_det_sympd(Rinv_aux);
          WtRW_aux.zeros();
          WtRW.zeros();
          for (int l = 0; l < L; ++l) {
            for (int t = 0; t < T; ++t) {
              indWUvec = indWUmat.col(l * T + t);
              vn = Wtls(indWUvec);
              WtRW_aux += vn.t() * Rinv_aux * vn;
              WtRW += vn.t() * Rinv * vn;
            }
          }
          A =
            (TL * Rlogdet_aux - arma::as_scalar(WtRW_aux)) / 2 +
            da(m) * ldecay_aux(m) - db(m) * decay_aux(m) -
            ((TL * Rlogdet - arma::as_scalar(WtRW)) / 2 +
            da(m) * ldecay(m) - db(m) * decay(m));
        }

        if (log(R::runif(0, 1)) <= A) {
          ++accept(m);
          decay(m) = decay_aux(m);
          ldecay(m) = ldecay_aux(m);
          Rinv = Rinv_aux;
          Rlogdet = Rlogdet_aux;
          dd = dd_aux;
        }
      }

      // tune sd of the proposal for decay
      if (b == 0) {
        accept.zeros();
        total = 0;
      } else if ((b < 1) && (++total % 25 == 0)) {
        r = accept / total;
        for (int m = 0; m < M; ++m) {
          if (r(m) > 0.33) {
            lsd(m) += 1 / sqrt(total / 25);
          } else {
            lsd(m) -= 1 / sqrt(total / 25);
          }
          sd(m) = exp(lsd(m));
        }
      }
    }

    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, arma::span(k, k + n * TL - 1)) = Wtls.t();
      keep(b / nThin - 1, k + n * TL) = a11;
      keep(b / nThin - 1, arma::span(k + n * TL + 1, k + n * TL + M)) = decay.t();
    }
  }

  return keep;
}



// @description Iterations of the Metropolis-within-Gibbs algorithm for the
//   logistic regression model.
// @param N number of observations.
// @param k number of covariates.
// @param Nl number of observations in l = 1 or l = 2.
// @param kl number of covariates for l = 1 or l = 2 models.
// @param Nweak number of weak records.
// @param X,X0,X1 matrix of covariates (X0 and X1 when weak records are all 0
//   or all 1, respectively).
// @param beta,betal1,betal2 vector initial value of beta, betal1 and betal2.
// @param weak,weak1,weak2 vector of positions of weak records (weak1 and weak2
//   position in lag1 and lag2).
// @param indLag1,indLag2,indLag12 vector of positions of columns with
//   lag1, lag2, or lag1:lag2.
// @param prob vector of probabilities that the observed indicator is 1.
// @param abLim vector indicating a record (1) or non-record (-1).
// @param V matrix prior var-cov of beta.
// @param keep matrix of beta's to keep.
// @param ... other arguments.
// [[Rcpp::export]]
arma::mat bglm(
    const int N, const int k, const arma::vec Nweak,
    arma::mat X, const arma::mat X0, const arma::mat X1,
    arma::mat beta, arma::mat Wtls,
    double a11, double a22, double a21,
    arma::vec decay, const int M,
    const arma::cube dist, const int decayPrior,
    const int n, const int T, const int L,
    const arma::uvec site, const arma::uvec year, const arma::uvec day,
    const arma::uvec weak1, const arma::uvec weak2,
    const arma::uvec Weak1, const arma::uvec Weak2,
    const arma::vec noNA1, const arma::vec noNA2,
    const arma::uvec indLag1, const arma::uvec indLag2,
    const arma::mat prob, arma::mat abLim,
    const arma::uvec weak1l0, const arma::uvec weak2l0,
    const arma::vec prob1l0, const arma::vec prob2l0,
    const arma::mat V,
    const double na, const double nb,
    const double diag_b,
    const arma::vec da, const arma::vec db,
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {

  const int TL = T * L;
  arma::mat Xb = X * beta;

  arma::mat Z(N, 2, arma::fill::zeros); // the first sampled (no need to specify)
  arma::mat auxZ = Z - Xb;
  auxZ.col(0) -= a11 * Wtls.col(0);
  auxZ.col(1) -= a21 * Wtls.col(0) + a22 * Wtls.col(1);

  arma::mat In(n, n, arma::fill::eye);
  arma::mat IM(M, M, arma::fill::eye);
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaW(n, n);

  // GPs
  double delta, chi;
  arma::vec vn(n);

  arma::uword ind0w = 0;
  arma::uword ind1w = 1;
  arma::uvec ind0 = { 0 };
  arma::uvec ind1 = { 1 };
  arma::uvec indWUvec(n);
  arma::umat indWUmat(n, TL);
  for (int l = 0; l < L; ++l) {
    for (int t = 0; t < T; ++t) {
      indWUmat.col(l * T + t) = arma::find((year == t) && (day == l));
    }
  }

  // decay (starter)
  arma::vec accept(M, arma::fill::zeros);
  int total = 0;
  arma::vec r(M);
  arma::vec sd(M, arma::fill::ones);
  arma::vec lsd(M, arma::fill::zeros); //log(sd);

  arma::mat Rinv(n, n);
  double Rlogdet;

  arma::mat dd_aux(n, n);
  arma::mat dd(n, n, arma::fill::zeros);
  for (int m = 0; m < M; ++m) {
    dd -= decay(m) * dist.slice(m);
  }
  Rinv = arma::inv_sympd(exp(dd));
  Rlogdet = arma::log_det_sympd(Rinv);

  arma::vec decay_aux(M);
  arma::vec ldecay_aux(M, arma::fill::zeros);
  arma::vec ldecay = log(decay);
  arma::mat Rinv_aux(n, n);
  double Rlogdet_aux;

  arma::mat WtRW(1, 1), WtRW_aux(1, 1);
  double A = 0;

  // weak
  arma::uword indw;
  arma::uvec ind;

  // LOOP
  for (int b = 1 - nBurnin; b <= nSims; ++b) {

    // check
    if (b % nReport == 0) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
      if ((b < 1) && (total != 0)) {
        Rcpp::Rcout << "Acceptance rate : " << accept.t() / total << "\n";
      } else if (b < 1) {
        Rcpp::Rcout << "Acceptance rate : " << 0 << "\n";
      } else {
        Rcpp::Rcout << "Acceptance rate : " <<  accept.t() / b << "\n";
      }
    }

    // weak records
    auxZ += Xb;
    for (int i = 0; i < Nweak(0); i++) { // tx and lag.tx
      indw = weak1(i);
      if (R::runif(0, 1) < prob(indw, ind0w)) {
        abLim(indw, ind0w) = 1;
        if (noNA1(i)) {
          ind = Weak1(i);
          X.submat(ind, indLag1) = X1.submat(ind, indLag1);
        }
      } else {
        abLim(indw, ind0w) = -1;
        if (noNA1(i)) {
          ind = Weak1(i);
          X.submat(ind, indLag1)  = X0.submat(ind, indLag1);
        }
      }
    }
    if (Nweak(2) > 0) {
    for (int i = 0; i < Nweak(2); i++) { // lag.tx l = 0
      ind = weak1l0(i);
      if (R::runif(0, 1) < prob1l0(i)) {
        X.submat(ind, indLag1) = X1.submat(ind, indLag1);
      } else {
        X.submat(ind, indLag1) = X0.submat(ind, indLag1);
      }
    }
    }
    for (int i = 0; i < Nweak(1); i++) { // tn and lag.tn
      indw = weak2(i);
      if (R::runif(0, 1) < prob(indw, ind1w)) {
        abLim(indw, ind1w) = 1;
        if (noNA2(i)) {
          ind = Weak2(i);
          X.submat(ind, indLag2) = X1.submat(ind, indLag2);
        }
      } else {
        abLim(indw, ind1w) = -1;
        if (noNA2(i)) {
          ind = Weak2(i);
          X.submat(ind, indLag2)  = X0.submat(ind, indLag2);
        }
      }
    }
    if (Nweak(3) > 0) {
    for (int i = 0; i < Nweak(3); i++) { // lag.tn l = 0
      ind = weak2l0(i);
      if (R::runif(0, 1) < prob2l0(i)) {
        X.submat(ind, indLag2) = X1.submat(ind, indLag2);
      } else {
        X.submat(ind, indLag2) = X0.submat(ind, indLag2);
      }
    }
    }
    Xb = X * beta;
    auxZ -= Xb;

    // Z x
    auxZ.col(0) -= Z.col(0);
    Z.col(0) = abLim.col(0) % rtnorm(N, -abLim.col(0) % auxZ.col(0));
    auxZ.col(0) += Z.col(0);

    // Z n
    auxZ.col(1) -= Z.col(1);
    Z.col(1) = abLim.col(1) % rtnorm(N, -abLim.col(1) % auxZ.col(1));
    auxZ.col(1) += Z.col(1);

    // beta x
    auxZ.col(0) += Xb.col(0);
    OmegaBeta = arma::inv_sympd(X.t() * X + V);
    beta.col(0) = OmegaBeta * (X.t() * auxZ.col(0)) + // prior 0k
      arma::chol(OmegaBeta, "lower") * arma::randn(k);
    Xb.col(0) = X * beta.col(0);
    auxZ.col(0) -= Xb.col(0);

    // beta n
    auxZ.col(1) += Xb.col(1);
    //// OmegaBeta = arma::inv_sympd(X.t() * X + V); // the same
    beta.col(1) = OmegaBeta * (X.t() * auxZ.col(1)) + // prior 0k
      arma::chol(OmegaBeta, "lower") * arma::randn(k);
    Xb.col(1) = X * beta.col(1);
    auxZ.col(1) -= Xb.col(1);

    // W1tl(s)
    auxZ.col(0) += a11 * Wtls.col(0);
    auxZ.col(1) += a21 * Wtls.col(0);
    OmegaW = arma::inv_sympd((a11 * a11 + a21 * a21) * In + Rinv);
    for (int l = 0; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        Wtls(indWUvec, ind0) = a11 * auxZ(indWUvec, ind0) + a21 * auxZ(indWUvec, ind1);
        Wtls(indWUvec, ind0) = OmegaW * Wtls(indWUvec, ind0) + arma::chol(OmegaW, "lower") * arma::randn(n);
      }
    }
    auxZ.col(0) -= a11 * Wtls.col(0);
    auxZ.col(1) -= a21 * Wtls.col(0);

    // W2tl(s)
    auxZ.col(1) += a22 * Wtls.col(1);
    OmegaW = arma::inv_sympd(a22 * a22 * In + Rinv);
    for (int l = 0; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        Wtls(indWUvec, ind1) = a22 * auxZ(indWUvec, ind1);
        Wtls(indWUvec, ind1) = OmegaW * Wtls(indWUvec, ind1) + arma::chol(OmegaW, "lower") * arma::randn(n);
      }
    }
    auxZ.col(1) -= a22 * Wtls.col(1);

    // a11
    auxZ.col(0) += a11 * Wtls.col(0);
    delta = 1 / (arma::dot(Wtls.col(0), Wtls.col(0)) + diag_b);
    chi   = arma::dot(Wtls.col(0), auxZ.col(0)); // + 0 * diag_b;
    a11   = rtnorm_numeric(delta * chi, sqrt(delta));
    auxZ.col(0) -= a11 * Wtls.col(0);

    // a22
    auxZ.col(1) += a22 * Wtls.col(1);
    delta = 1 / (arma::dot(Wtls.col(1), Wtls.col(1)) + diag_b);
    chi   = arma::dot(Wtls.col(1), auxZ.col(1)); // + 0 * diag_b;
    a22   = rtnorm_numeric(delta * chi, sqrt(delta));
    auxZ.col(1) -= a22 * Wtls.col(1);

    // a21
    auxZ.col(1) += a21 * Wtls.col(0);
    delta = 1 / (arma::dot(Wtls.col(0), Wtls.col(0)) + nb);
    chi   = arma::dot(Wtls.col(0), auxZ.col(1)) + na * nb;
    a21   = R::rnorm(delta * chi, sqrt(delta));
    auxZ.col(1) -= a21 * Wtls.col(0);

    // decay
    if (decayPrior != 1) {
    for (int m = 0; m < M; ++m) {
      if (decayPrior == 2) {
        decay_aux(m) = rtnorm1(decay(m), sd(m), da(m), db(m));
        dd_aux       = dd + (decay(m) - decay_aux(m)) * dist.slice(m);
        Rinv_aux     = arma::inv_sympd(exp(dd_aux));
        Rlogdet_aux  = arma::log_det_sympd(Rinv_aux);
        WtRW_aux.zeros();
        WtRW.zeros();
        for (int l = 0; l < L; ++l) {
          for (int t = 0; t < T; ++t) {
            indWUvec = indWUmat.col(l * T + t);
            vn = Wtls(indWUvec, ind0);
            WtRW_aux += vn.t() * Rinv_aux * vn;
            WtRW += vn.t() * Rinv * vn;
            vn = Wtls(indWUvec, ind1);
            WtRW_aux += vn.t() * Rinv_aux * vn;
            WtRW += vn.t() * Rinv * vn;
          }
        }
        A =
          (TL * Rlogdet_aux - arma::as_scalar(WtRW_aux) / 2) -
          (TL * Rlogdet     - arma::as_scalar(WtRW)     / 2) +
          dtnorm1(decay(m)    , decay_aux(m), sd(m), da(m), db(m)) -
          dtnorm1(decay_aux(m), decay(m)    , sd(m), da(m), db(m));
      } else if (decayPrior == 3) {
        ldecay_aux(m) = R::rnorm(ldecay(m), sd(m));
        decay_aux(m)  = exp(ldecay_aux(m));
        dd_aux       = dd + (decay(m) - decay_aux(m)) * dist.slice(m);
        Rinv_aux     = arma::inv_sympd(exp(dd_aux));
        Rlogdet_aux  = arma::log_det_sympd(Rinv_aux);
        WtRW_aux.zeros();
        WtRW.zeros();
        for (int l = 0; l < L; ++l) {
          for (int t = 0; t < T; ++t) {
            indWUvec = indWUmat.col(l * T + t);
            vn = Wtls(indWUvec, ind0);
            WtRW_aux += vn.t() * Rinv_aux * vn;
            WtRW += vn.t() * Rinv * vn;
            vn = Wtls(indWUvec, ind1);
            WtRW_aux += vn.t() * Rinv_aux * vn;
            WtRW += vn.t() * Rinv * vn;
          }
        }
        A =
          (TL * Rlogdet_aux - arma::as_scalar(WtRW_aux) / 2) +
          da(m) * ldecay_aux(m) - db(m) * decay_aux(m) -
          ((TL * Rlogdet - arma::as_scalar(WtRW) / 2) +
          da(m) * ldecay(m) - db(m) * decay(m));
      }

      if (log(R::runif(0, 1)) <= A) {
        ++accept(m);
        decay(m) = decay_aux(m);
        ldecay(m) = ldecay_aux(m);
        Rinv = Rinv_aux;
        Rlogdet = Rlogdet_aux;
        dd = dd_aux;
      }
    }

      // tune sd of the proposal for decay
      if (b == 0) {
        accept.zeros();
        total = 0;
      } else if ((b < 1) && (++total % 25 == 0)) {
        r = accept / total;
        for (int m = 0; m < M; ++m) {
          if (r(m) > 0.33) {
            lsd(m) += 1 / sqrt(total / 25);
          } else {
            lsd(m) -= 1 / sqrt(total / 25);
          }
          sd(m) = exp(lsd(m));
        }
      }
    }

    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.col(0).t();
      keep(b / nThin - 1, arma::span(k, 2 * k - 1)) = beta.col(1).t();
      keep(b / nThin - 1, arma::span(2 * k, 2 * k + n * TL - 1)) = Wtls.col(0).t();
      keep(b / nThin - 1, arma::span(2 * k + n * TL, 2 * k + 2 * n * TL - 1)) = Wtls.col(1).t();
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL, 2 * k + 2 * n * TL + 2)) = { a11, a22, a21 };
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL + 3, 2 * k + 2 * n * TL + M + 2)) = decay.t();
    }
  }

  return keep;
}



// [[Rcpp::export]]
arma::mat bglm2(
    const int N, const int k, const arma::vec Nweak,
    arma::mat X, const arma::mat X0, const arma::mat X1,
    arma::mat beta, arma::mat Wtls,
    arma::vec a11, arma::vec a22, arma::vec a21,
    arma::vec hpA, // mu , sigma , decay for 11 22 21
    const arma::mat onen, const int q, // onen is the matrix of covariates in A(s)
    arma::vec decay, const int M,
    const arma::cube dist, const int decayPrior,
    const int n, const int T, const int L,
    const arma::uvec site, const arma::uvec year, const arma::uvec day,
    const arma::uvec weak1, const arma::uvec weak2,
    const arma::uvec Weak1, const arma::uvec Weak2,
    const arma::vec noNA1, const arma::vec noNA2,
    const arma::uvec indLag1, const arma::uvec indLag2,
    const arma::mat prob, arma::mat abLim,
    const arma::uvec weak1l0, const arma::uvec weak2l0,
    const arma::vec prob1l0, const arma::vec prob2l0,
    const arma::mat Vk, const arma::mat Vq,
    const double na, const double nb,
    const double ga, const double gb,
    const arma::vec da, const arma::vec db,
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {

  const int TL = T * L;
  arma::mat Xb = X * beta;

  arma::mat Z(N, 2, arma::fill::zeros); // the first sampled (no need to specify)
  arma::mat auxZ = Z - Xb;
  auxZ.col(0) -= a11(site) % Wtls.col(0);
  auxZ.col(1) -= a21(site) % Wtls.col(0) + a22(site) % Wtls.col(1);

  arma::mat In(n, n, arma::fill::eye);
  arma::mat IM(M, M, arma::fill::eye);
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaW(n, n);
  arma::mat cholOmegaW(n, n);

  // A GPs
  int Nn = N / n;
  int count0, count1;
  arma::vec la11 = log(a11);
  arma::vec la22 = log(a22);
  arma::vec process(n);
  arma::vec rinv_A(n);
  //arma::vec onen(n, arma::fill::ones); // now intercept + covariates
  arma::mat Rinv_A = arma::inv_sympd(exp(- hpA(2) * dist.slice(0)));
  arma::mat oneR_Aone = onen.t() * Rinv_A * onen;

  double la_aux;
  double a_aux;

  arma::vec accept_A(2 * n, arma::fill::zeros);
  int total_A = 0;
  arma::vec r_A(2 * n);
  arma::vec sd_A(2 * n, arma::fill::ones);
  arma::vec lsd_A(2 * n, arma::fill::zeros); //log(sd_A);

  // GPs
  arma::mat Delta(q, q);
  double delta, chi;
  arma::vec vn(n);
  arma::vec mua11 = onen * hpA(arma::span(0, q - 1));
  arma::vec mua22 = onen * hpA(arma::span(q + 2, 2 * q + 1));
  arma::vec mua21 = onen * hpA(arma::span(2 * q + 4, 3 * q + 3));

  arma::uword ind0w = 0;
  arma::uword ind1w = 1;
  arma::uvec ind0 = { 0 };
  arma::uvec ind1 = { 1 };
  arma::uvec indWUvec(n);
  arma::umat indWUmat(n, TL);
  for (int l = 0; l < L; ++l) {
    for (int t = 0; t < T; ++t) {
      indWUmat.col(l * T + t) = arma::find((year == t) && (day == l));
    }
  }

  // decay (starter)
  arma::vec accept(M, arma::fill::zeros);
  int total = 0;
  arma::vec r(M);
  arma::vec sd(M, arma::fill::ones);
  arma::vec lsd(M, arma::fill::zeros); //log(sd);

  arma::mat Rinv(n, n);
  double Rlogdet;

  arma::mat dd_aux(n, n);
  arma::mat dd(n, n, arma::fill::zeros);
  for (int m = 0; m < M; ++m) {
    dd -= decay(m) * dist.slice(m);
  }
  Rinv = arma::inv_sympd(exp(dd));
  Rlogdet = arma::log_det_sympd(Rinv);

  arma::vec decay_aux(M);
  arma::vec ldecay_aux(M, arma::fill::zeros);
  arma::vec ldecay = log(decay);
  arma::mat Rinv_aux(n, n);
  double Rlogdet_aux;

  arma::mat WtRW(1, 1), WtRW_aux(1, 1);
  double A = 0;

  // weak
  arma::uword indw;
  arma::uvec ind;

  // LOOP
  for (int b = 1 - nBurnin; b <= nSims; ++b) {

    // check
    if (b % nReport == 0) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
      if ((b < 1) && (total != 0)) {
        Rcpp::Rcout << "Acceptance rate for decay : " << accept.t() / total << ", for A : " << accept_A.t() / total_A << "\n";
      } else if (b == 1) {
        Rcpp::Rcout << "Acceptance rate : -- " << "\n";
      } else {
        Rcpp::Rcout << "Acceptance rate for decay : " << accept.t() / (b - 1) << ", for A : " << accept_A.t() / (b - 1) << "\n";
      }
    }

    // weak records
    auxZ += Xb;
    for (int i = 0; i < Nweak(0); i++) { // tx and lag.tx
      indw = weak1(i);
      if (R::runif(0, 1) < prob(indw, ind0w)) {
        abLim(indw, ind0w) = 1;
        if (noNA1(i)) {
          ind = Weak1(i);
          X.submat(ind, indLag1) = X1.submat(ind, indLag1);
        }
      } else {
        abLim(indw, ind0w) = -1;
        if (noNA1(i)) {
          ind = Weak1(i);
          X.submat(ind, indLag1)  = X0.submat(ind, indLag1);
        }
      }
    }
    if (Nweak(2) > 0) {
    for (int i = 0; i < Nweak(2); i++) { // lag.tx l = 0
      ind = weak1l0(i);
      if (R::runif(0, 1) < prob1l0(i)) {
        X.submat(ind, indLag1) = X1.submat(ind, indLag1);
      } else {
        X.submat(ind, indLag1) = X0.submat(ind, indLag1);
      }
    }
    }
    for (int i = 0; i < Nweak(1); i++) { // tn and lag.tn
      indw = weak2(i);
      if (R::runif(0, 1) < prob(indw, ind1w)) {
        abLim(indw, ind1w) = 1;
        if (noNA2(i)) {
          ind = Weak2(i);
          X.submat(ind, indLag2) = X1.submat(ind, indLag2);
        }
      } else {
        abLim(indw, ind1w) = -1;
        if (noNA2(i)) {
          ind = Weak2(i);
          X.submat(ind, indLag2)  = X0.submat(ind, indLag2);
        }
      }
    }
    if (Nweak(3) > 0) {
    for (int i = 0; i < Nweak(3); i++) { // lag.tn l = 0
      ind = weak2l0(i);
      if (R::runif(0, 1) < prob2l0(i)) {
        X.submat(ind, indLag2) = X1.submat(ind, indLag2);
      } else {
        X.submat(ind, indLag2) = X0.submat(ind, indLag2);
      }
    }
    }
    Xb = X * beta;
    auxZ -= Xb;

    // Z x
    auxZ.col(0) -= Z.col(0);
    Z.col(0) = abLim.col(0) % rtnorm(N, -abLim.col(0) % auxZ.col(0));
    auxZ.col(0) += Z.col(0);

    // Z n
    auxZ.col(1) -= Z.col(1);
    Z.col(1) = abLim.col(1) % rtnorm(N, -abLim.col(1) % auxZ.col(1));
    auxZ.col(1) += Z.col(1);

    // beta x
    auxZ.col(0) += Xb.col(0);
    OmegaBeta = arma::inv_sympd(X.t() * X + Vk);
    beta.col(0) = OmegaBeta * (X.t() * auxZ.col(0)) + // prior 0k
      arma::chol(OmegaBeta, "lower") * arma::randn(k);
    Xb.col(0) = X * beta.col(0);
    auxZ.col(0) -= Xb.col(0);

    // beta n
    auxZ.col(1) += Xb.col(1);
    //// OmegaBeta = arma::inv_sympd(X.t() * X + Vk); // the same
    beta.col(1) = OmegaBeta * (X.t() * auxZ.col(1)) + // prior 0k
      arma::chol(OmegaBeta, "lower") * arma::randn(k);
    Xb.col(1) = X * beta.col(1);
    auxZ.col(1) -= Xb.col(1);

    // W1tl(s)
    auxZ.col(0) += a11(site) % Wtls.col(0);
    auxZ.col(1) += a21(site) % Wtls.col(0);
    OmegaW = arma::inv_sympd(arma::diagmat(a11 % a11 + a21 % a21) + Rinv);
    cholOmegaW = arma::chol(OmegaW, "lower");
    for (int l = 0; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        Wtls(indWUvec, ind0) = a11 % auxZ(indWUvec, ind0) + a21 % auxZ(indWUvec, ind1);
        Wtls(indWUvec, ind0) = OmegaW * Wtls(indWUvec, ind0) + cholOmegaW * arma::randn(n);
      }
    }
    auxZ.col(0) -= a11(site) % Wtls.col(0);
    auxZ.col(1) -= a21(site) % Wtls.col(0);

    // W2tl(s)
    auxZ.col(1) += a22(site) % Wtls.col(1);
    OmegaW = arma::inv_sympd(arma::diagmat(a22 % a22) + Rinv);
    cholOmegaW = arma::chol(OmegaW, "lower");
    for (int l = 0; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        Wtls(indWUvec, ind1) = a22 % auxZ(indWUvec, ind1);
        Wtls(indWUvec, ind1) = OmegaW * Wtls(indWUvec, ind1) + cholOmegaW * arma::randn(n);
      }
    }
    auxZ.col(1) -= a22(site) % Wtls.col(1);

    // a11
    auxZ.col(0) += a11(site) % Wtls.col(0);
    count1 = 0;
    // IMPORTANT: data ordered by site (first all s1, then s2, etc.)
    for (int i = 0; i < n; ++i) {
      la_aux = R::rnorm(la11(i), sd_A(i));
      a_aux  = exp(la_aux);
      count0 = count1;
      count1 += Nn;

      delta = arma::accu(Wtls(arma::span(count0, count1 - 1), 0) %
                         Wtls(arma::span(count0, count1 - 1), 0));
      chi   = arma::accu(Wtls(arma::span(count0, count1 - 1), 0) %
                         auxZ(arma::span(count0, count1 - 1), 0)) / delta;

      A = dnorm1(a_aux, chi, delta) - dnorm1(a11(i), chi, delta);

      process = mua11 - la11;
      process.shed_row(i);
      rinv_A = Rinv_A.col(i);
      rinv_A.shed_row(i);

      delta = hpA(q) * Rinv_A(i, i);
      chi   = mua11(i) + arma::dot(rinv_A, process) / Rinv_A(i, i);

      A += dnorm1(la_aux, chi, delta) - dnorm1(la11(i), chi, delta);

      if (log(R::runif(0, 1)) <= A) {
        ++accept_A(i);
        a11(i) = a_aux;
        la11(i) = la_aux;
      }
    }
    auxZ.col(0) -= a11(site) % Wtls.col(0);

    // hp a11
    //// mu
    Delta  = arma::inv_sympd(oneR_Aone * hpA(q) + Vq);
    hpA(arma::span(0, q - 1)) =
      Delta * (onen.t() * Rinv_A * la11 * hpA(q)) + // prior 0k
      arma::chol(Delta, "lower") * arma::randn(q);
    mua11 = onen * hpA(arma::span(0, q - 1));

    //// prec
    process = la11 - mua11;
    hpA(q) = R::rgamma(n / 2 + ga,
        1 / (arma::as_scalar(process.t() * Rinv_A * process) / 2 + gb));

    // a22
    auxZ.col(1) += a22(site) % Wtls.col(1);
    count1 = 0;
    // IMPORTANT: data ordered by site (first all s1, then s2, etc.)
    for (int i = 0; i < n; ++i) {
      la_aux = R::rnorm(la22(i), sd_A(i + n));
      a_aux  = exp(la_aux);
      count0 = count1;
      count1 += Nn;

      delta = arma::accu(Wtls(arma::span(count0, count1 - 1), 1) %
                         Wtls(arma::span(count0, count1 - 1), 1));
      chi   = arma::accu(Wtls(arma::span(count0, count1 - 1), 1) %
                         auxZ(arma::span(count0, count1 - 1), 1)) / delta;

      A = dnorm1(a_aux, chi, delta) - dnorm1(a22(i), chi, delta);

      process = mua22 - la22;
      process.shed_row(i);
      rinv_A = Rinv_A.col(i);
      rinv_A.shed_row(i);

      delta = hpA(2 * q + 2) * Rinv_A(i, i);
      chi   = mua22(i) + arma::dot(rinv_A, process) / Rinv_A(i, i);

      A += dnorm1(la_aux, chi, delta) - dnorm1(la22(i), chi, delta);

      if (log(R::runif(0, 1)) <= A) {
        ++accept_A(i + n);
        a22(i) = a_aux;
        la22(i) = la_aux;
      }
    }
    auxZ.col(1) -= a22(site) % Wtls.col(1);

    // tuning A
    if (b == 0) {
      accept_A.zeros();
      total_A = 0;
    } else if ((b < 1) && (++total_A % 25 == 0)) {
      r_A = accept_A / total_A;
      for (int i = 0; i < 2 * n; ++i) {
        if (r_A(i) > 0.33) {
          lsd_A(i) += 1 / sqrt(total_A / 25);
        } else {
          lsd_A(i) -= 1 / sqrt(total_A / 25);
        }
        sd_A(i) = exp(lsd_A(i));
      }
    }

    // hp a22
    //// mu
    Delta  = arma::inv_sympd(oneR_Aone * hpA(2 * q + 2) + Vq);
    hpA(arma::span(q + 2, 2 * q + 1)) =
      Delta * (onen.t() * Rinv_A * la22 * hpA(2 * q + 2)) + // prior 0k
      arma::chol(Delta, "lower") * arma::randn(q);
    mua22 = onen * hpA(arma::span(q + 2, 2 * q + 1));

    //// prec
    process = la22 - mua22;
    hpA(2 * q + 2) = R::rgamma(n / 2 + ga,
        1 / (arma::as_scalar(process.t() * Rinv_A * process) / 2 + gb));

    // a21
    auxZ.col(1) += a21(site) % Wtls.col(0);
    OmegaW = hpA(3 * q + 4) * Rinv_A;
    a21 = OmegaW * mua21;
    count1 = 0;
    // IMPORTANT: data ordered by site (first all s1, then s2, etc.)
    for (int i = 0; i < n; ++i) {
      count0 = count1;
      count1 += Nn;
      OmegaW(i, i) += arma::accu(Wtls(arma::span(count0, count1 - 1), 0) % Wtls(arma::span(count0, count1 - 1), 0));
      a21(i) += arma::accu(Wtls(arma::span(count0, count1 - 1), 0) % auxZ(arma::span(count0, count1 - 1), 1));
    }
    OmegaW = arma::inv_sympd(OmegaW);
    a21 = OmegaW * a21 + arma::chol(OmegaW, "lower") * arma::randn(n);
    auxZ.col(1) -= a21(site) % Wtls.col(0);

    // hp a21
    //// mu
    Delta  = arma::inv_sympd(oneR_Aone * hpA(3 * q + 4) + Vq);
    hpA(arma::span(2 * q + 4, 3 * q + 3)) =
      Delta * (onen.t() * Rinv_A * a21 * hpA(3 * q + 4)) + // prior 0k
      arma::chol(Delta, "lower") * arma::randn(q);
    mua21 = onen * hpA(arma::span(2 * q + 4, 3 * q + 3));

    //// prec
    process = a21 - mua21;
    hpA(3 * q + 4) = R::rgamma(n / 2 + ga,
       1 / (arma::as_scalar(process.t() * Rinv_A * process) / 2 + gb));

    // decay
    if (decayPrior != 1) {
      for (int m = 0; m < M; ++m) {
        if (decayPrior == 2) {
          decay_aux(m) = rtnorm1(decay(m), sd(m), da(m), db(m));
          dd_aux       = dd + (decay(m) - decay_aux(m)) * dist.slice(m);
          Rinv_aux     = arma::inv_sympd(exp(dd_aux));
          Rlogdet_aux  = arma::log_det_sympd(Rinv_aux);
          WtRW_aux.zeros();
          WtRW.zeros();
          for (int l = 0; l < L; ++l) {
            for (int t = 0; t < T; ++t) {
              indWUvec = indWUmat.col(l * T + t);
              vn = Wtls(indWUvec, ind0);
              WtRW_aux += vn.t() * Rinv_aux * vn;
              WtRW += vn.t() * Rinv * vn;
              vn = Wtls(indWUvec, ind1);
              WtRW_aux += vn.t() * Rinv_aux * vn;
              WtRW += vn.t() * Rinv * vn;
            }
          }
          A =
            (TL * Rlogdet_aux - arma::as_scalar(WtRW_aux) / 2) -
            (TL * Rlogdet     - arma::as_scalar(WtRW)     / 2) +
            dtnorm1(decay(m)    , decay_aux(m), sd(m), da(m), db(m)) -
            dtnorm1(decay_aux(m), decay(m)    , sd(m), da(m), db(m));
        } else if (decayPrior == 3) {
          ldecay_aux(m) = R::rnorm(ldecay(m), sd(m));
          decay_aux(m)  = exp(ldecay_aux(m));
          dd_aux       = dd + (decay(m) - decay_aux(m)) * dist.slice(m);
          Rinv_aux     = arma::inv_sympd(exp(dd_aux));
          Rlogdet_aux  = arma::log_det_sympd(Rinv_aux);
          WtRW_aux.zeros();
          WtRW.zeros();
          for (int l = 0; l < L; ++l) {
            for (int t = 0; t < T; ++t) {
              indWUvec = indWUmat.col(l * T + t);
              vn = Wtls(indWUvec, ind0);
              WtRW_aux += vn.t() * Rinv_aux * vn;
              WtRW += vn.t() * Rinv * vn;
              vn = Wtls(indWUvec, ind1);
              WtRW_aux += vn.t() * Rinv_aux * vn;
              WtRW += vn.t() * Rinv * vn;
            }
          }
          A =
            (TL * Rlogdet_aux - arma::as_scalar(WtRW_aux) / 2) +
            da(m) * ldecay_aux(m) - db(m) * decay_aux(m) -
            ((TL * Rlogdet - arma::as_scalar(WtRW) / 2) +
            da(m) * ldecay(m) - db(m) * decay(m));
        }

        if (log(R::runif(0, 1)) <= A) {
          ++accept(m);
          decay(m) = decay_aux(m);
          ldecay(m) = ldecay_aux(m);
          Rinv = Rinv_aux;
          Rlogdet = Rlogdet_aux;
          dd = dd_aux;
        }
      }

      // tune sd of the proposal for decay
      if (b == 0) {
        accept.zeros();
        total = 0;
      } else if ((b < 1) && (++total % 25 == 0)) {
        r = accept / total;
        for (int m = 0; m < M; ++m) {
          if (r(m) > 0.33) {
            lsd(m) += 1 / sqrt(total / 25);
          } else {
            lsd(m) -= 1 / sqrt(total / 25);
          }
          sd(m) = exp(lsd(m));
        }
      }
    }

    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.col(0).t();
      keep(b / nThin - 1, arma::span(k, 2 * k - 1)) = beta.col(1).t();
      keep(b / nThin - 1, arma::span(2 * k, 2 * k + n * TL - 1)) = Wtls.col(0).t();
      keep(b / nThin - 1, arma::span(2 * k + n * TL, 2 * k + 2 * n * TL - 1)) = Wtls.col(1).t();
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL, 2 * k + 2 * n * TL + n - 1)) = a11.t();
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL + n, 2 * k + 2 * n * TL + 2 * n - 1)) = a22.t();
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL + 2 * n, 2 * k + 2 * n * TL + 3 * n - 1)) = a21.t();
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL + 3 * n, 2 * k + 2 * n * TL + 3 * n + 3 * (q + 2) - 1)) = hpA.t();
      keep(b / nThin - 1, arma::span(2 * k + 2 * n * TL + 3 * n + 3 * (q + 2), 2 * k + 2 * n * TL + 3 * n + 3 * (q + 2) + M - 1)) = decay.t();
    }
  }

  return keep;
}



// [[Rcpp::export]]
arma::cube predGlmBerKFCV(
    const int B, const arma::vec Nweak, const int newN,
    arma::mat X, const arma::mat X0, const arma::mat X1,
    const arma::mat beta1, const arma::mat beta2,
    const arma::uvec weak1, const arma::uvec weak2,
    const arma::uvec indLag1, const arma::uvec indLag2,
    const arma::vec prob1, const arma::vec prob2
) {

  arma::cube Xb(B, newN, 2);
  arma::uvec ind;

  for (int b = 0; b < B; ++b) {

    // weak records
    if (Nweak(0) > 0) {
    for (int i = 0; i < Nweak(0); i++) { // lag.tx
      ind = weak1(i);
      if (R::runif(0, 1) < prob1(i)) {
        X.submat(ind, indLag1) = X1.submat(ind, indLag1);
      } else {
        X.submat(ind, indLag1) = X0.submat(ind, indLag1);
      }
    }
    }
    if (Nweak(1) > 0) {
    for (int i = 0; i < Nweak(1); i++) { // lag.tn
      ind = weak2(i);
      if (R::runif(0, 1) < prob2(i)) {
        X.submat(ind, indLag2) = X1.submat(ind, indLag2);
      } else {
        X.submat(ind, indLag2) = X0.submat(ind, indLag2);
      }
    }
    }

    // Xb
    Xb.slice(0).row(b) = beta1.row(b) * X.t();
    Xb.slice(1).row(b) = beta2.row(b) * X.t();
  }

  return Xb;
}



// [[Rcpp::export]]
arma::cube predictBglm(
    const int B, const arma::vec Nweak, const int newN,
    const int n, const int newn, const int T, const int L,
    arma::mat X, const arma::mat X0, const arma::mat X1,
    arma::vec year, arma::vec day,
    arma::uvec newsite, arma::vec newyear, arma::vec newday,
    arma::cube dist,
    const arma::mat beta1, const arma::mat beta2,
    const arma::mat W1tls, const arma::mat W2tls,
    const arma::mat a, bool spA, arma::mat hpA,
    const arma::mat onen, const arma::mat onenew,
    const int q,
    const arma::mat decay, const int M,
    const arma::uvec weak1, const arma::uvec weak2,
    const arma::uvec indLag1, const arma::uvec indLag2,
    const arma::vec prob1, const arma::vec prob2
) {

  arma::cube Xb(B, newN, 2);

  arma::rowvec W1tls0(newN);
  arma::rowvec W2tls0(newN);

  arma::mat dd(newn + n, newn + n);
  arma::mat Sigma11(newn, newn);
  arma::mat Sigma12(newn, n);
  arma::mat Sigma22inv(n, n);
  arma::mat SigmaAux(n, newn); // transposed
  arma::mat Sigma(newn, newn);

  arma::uvec indWUvec(n);
  arma::umat indWUmat(n, T * L);
  arma::uvec indW0Uvec(newn);
  arma::umat indW0Umat(newn, T * L);
  for (int l = 0; l < L; ++l) {
    for (int t = 0; t < T; ++t) {
      indWUmat.col(l * T + t) = arma::find((year == t) && (day == l));
      indW0Umat.col(l * T + t) = arma::find((newyear == t) && (newday == l));
    }
  }

  arma::vec a11s(newn);
  arma::vec a22s(newn);
  arma::vec a21s(newn);
  arma::vec aAux(n);
  arma::mat Sigma11_A(newn, newn);
  arma::mat Sigma12_A(newn, n);
  arma::mat Sigma22inv_A(n, n);
  arma::mat SigmaAux_A(newn, n);
  arma::mat Sigma_A(newn, newn);

  if (spA) {
    dd = - hpA(0, q + 1) * dist.slice(0);
    Sigma11_A    =
      exp(dd.submat(0, 0, newn - 1, newn - 1));
    Sigma12_A    =
      exp(dd.submat(0, newn, newn - 1, newn + n - 1));
    Sigma22inv_A =
      arma::inv_sympd(exp(dd.submat(newn, newn, newn + n - 1, newn + n - 1)));
    SigmaAux_A   =
      Sigma12_A * Sigma22inv_A;
    Sigma_A      =
      arma::chol(Sigma11_A - SigmaAux_A * Sigma12_A.t(), "lower");
  }

  arma::uvec ind;
  arma::uvec indAux;

  for (int b = 0; b < B; ++b) {
    indAux = b;

    // weak records
    if (Nweak(0) > 0) {
    for (int i = 0; i < Nweak(0); i++) { // lag.tx
      ind = weak1(i);
      if (R::runif(0, 1) < prob1(i)) {
        X.submat(ind, indLag1) = X1.submat(ind, indLag1);
      } else {
        X.submat(ind, indLag1) = X0.submat(ind, indLag1);
      }
    }
    }
    if (Nweak(1) > 0) {
    for (int i = 0; i < Nweak(1); i++) { // lag.tn
      ind = weak2(i);
      if (R::runif(0, 1) < prob2(i)) {
        X.submat(ind, indLag2) = X1.submat(ind, indLag2);
      } else {
        X.submat(ind, indLag2)  = X0.submat(ind, indLag2);
      }
    }
    }

    // Wtls
    dd.zeros();
    for (int m = 0; m < M; ++m) {
      dd -= decay(b, m) * dist.slice(m);
    }
    Sigma11    = exp(dd.submat(0, 0, newn - 1, newn - 1));
    Sigma12    = exp(dd.submat(0, newn, newn - 1, newn + n - 1));
    Sigma22inv = arma::inv_sympd(exp(dd.submat(newn, newn, newn + n - 1, newn + n - 1)));
    SigmaAux   = (Sigma12 * Sigma22inv).t();
    Sigma      = Sigma11 - SigmaAux.t() * Sigma12.t();
    Sigma      = arma::chol(Sigma, "upper");

    for (int l = 0; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        indW0Uvec = indW0Umat.col(l * T + t);
        W1tls0(indW0Uvec) =
          W1tls.submat(indAux, indWUvec) * SigmaAux +
          arma::randn(newn).t() * Sigma;
        W2tls0(indW0Uvec) =
          W2tls.submat(indAux, indWUvec) * SigmaAux +
          arma::randn(newn).t() * Sigma;
      }
    }

    if (spA) {
      // a11s
      aAux = a(b, arma::span(0, n - 1)).t();
      a11s = onenew * hpA(b, arma::span(0, q - 1)).t() +
        SigmaAux_A * (log(aAux) - onen * hpA(b, arma::span(0, q - 1)).t()) +
        Sigma_A * arma::randn(newn) / sqrt(hpA(b, q));
      a11s = exp(a11s);

      // a22s
      aAux = a(b, arma::span(n, 2 * n - 1)).t();
      a22s = onenew * hpA(b, arma::span(q + 2, 2 * q + 1)).t() +
        SigmaAux_A * (log(aAux) - onen * hpA(b, arma::span(q + 2, 2 * q + 1)).t()) +
        Sigma_A * arma::randn(newn) / sqrt(hpA(b, 2 * q + 2));
      a22s = exp(a22s);

      // a21s
      aAux = a(b, arma::span(2 * n, 3 * n - 1)).t();
      a21s = onenew * hpA(b, arma::span(2 * q + 4, 3 * q + 3)).t() +
        SigmaAux_A * (aAux - onen * hpA(b, arma::span(2 * q + 4, 3 * q + 3)).t()) +
        Sigma_A * arma::randn(newn) / sqrt(hpA(b, 3 * q + 4));

      // Xb
      Xb.slice(0).row(b) = beta1.row(b) * X.t() + a11s(newsite).t() % W1tls0;
      Xb.slice(1).row(b) = beta2.row(b) * X.t() + a21s(newsite).t() % W1tls0 + a22s(newsite).t() % W2tls0;
    } else {
      // Xb
      Xb.slice(0).row(b) = beta1.row(b) * X.t() + a(b, 0) * W1tls0;
      Xb.slice(1).row(b) = beta2.row(b) * X.t() + a(b, 2) * W1tls0 + a(b, 1) * W2tls0;
    }
  }

  return Xb;
}
