functions {
  vector sump(vector c, vector r, int order) {
    int N = num_elements(r);
    vector[N] s;
    for (n in 1:N) {
      s[n] = c[1];
      for (i in 2:order) {
        s[n] = s[n] + c[i]*r[n]^(i-1);
      }
    }
    return s;
  }

  vector gp_pred_rng(row_vector[] xy_pred,
                     vector v_obs,
                     row_vector[] xy_obs,
                     real alpha,
                     real rho,
                     real sigma) {

    int N_pred;
    int N;
    N_pred = dims(xy_pred)[1];
    N = rows(v_obs);
    {
      vector[N_pred] v_pred;
          
      vector[N] K_div_v;
      matrix[N, N_pred] k_xy_xy_pred;
      matrix[N, N_pred] qv_pred;
      vector[N_pred] v_pred_mu;
      matrix[N_pred, N_pred] cov_v_pred;
      matrix[N_pred, N_pred] nug_pred;
      matrix[N, N] cov;
      matrix[N, N] L_cov;
      
      cov = cov_exp_quad(xy_obs, alpha, rho);
      for (n in 1:N) {
        cov[n, n] = cov[n,n] + square(sigma);
      }
      L_cov = cholesky_decompose(cov);
      K_div_v = mdivide_left_tri_low(L_cov, v_obs);
      K_div_v = mdivide_right_tri_low(K_div_v',L_cov)';
      k_xy_xy_pred = cov_exp_quad(xy_obs, xy_pred, alpha, rho);
      v_pred_mu = k_xy_xy_pred' * K_div_v; 
      qv_pred = mdivide_left_tri_low(L_cov, k_xy_xy_pred);
      cov_v_pred = cov_exp_quad(xy_pred, alpha, rho) - qv_pred' * qv_pred;
      nug_pred = diag_matrix(rep_vector(1e-10,N_pred));
      
      v_pred = multi_normal_rng(v_pred_mu, cov_v_pred + nug_pred);
      
      return v_pred;
    }
  }
  vector inv_gamma_tail(vector pars, vector theta, real[] x_r, int[] x_i) {
    vector[2] deltas;
    deltas[1] = inv_gamma_cdf(theta[1], exp(pars[1]), exp(pars[2])) - 0.01;
    deltas[2] = 1 - inv_gamma_cdf(theta[2], exp(pars[1]), exp(pars[2])) - 0.01;
    return deltas;
  }

}

data {
  int<lower=1> order; //order for polynomial representation
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
  vector[N] v;  //measured los velocity
  vector[N] dv; //error estimate on v
  real phi0;
  real<lower=0.> sd_phi0;
  real ci0;
  real<lower=0.> sd_ci0;
  real<lower=0.> sd_kc0;
  real<lower=0.> r_norm;
  real<lower=0.> v_norm;
  int<lower=1> N_xy; //number of points to fill in
  int<lower=1> N_r;
  vector[N_xy] x_pred;
  vector[N_xy] y_pred;
}
transformed data {
  vector[N_xy] r_pred = sqrt(x_pred .* x_pred + y_pred .* y_pred);
  vector[N_xy] theta_pred;
  row_vector[2] xy_pred[N_xy];
  vector[2] pars;
  vector[2] p_guess;
  vector[2] tails;
  vector[2] p_sol;
  real x_r[0];
  int x_i[0];
  real a;
  real b;
    
  xy_pred[:, 1] = to_array_1d(x_pred);
  xy_pred[:, 2] = to_array_1d(y_pred);
  for (n in 1:N_xy) {
    theta_pred[n] = atan2(-x_pred[n], y_pred[n]);
  }
  
  tails[1] = 2./r_norm;
  tails[2] = 1.;
  p_guess = [log(r_norm/2.), log(2)]';
  p_sol = algebra_solver(inv_gamma_tail, p_guess, tails, x_r, x_i);
  a = exp(p_sol[1]);
  b = exp(p_sol[2]);
  print("a = ",a);
  print("b = ",b);
  
}
parameters {
  real phi;
  real<lower=0., upper=1.> cos_i;  // cosine disk inclination
  real x_c;  //kinematic centers
  real y_c;
  real v_sys;     //system velocity offset (should be small)
  vector[order] c_r;
  vector[N] v_los;  //latent "real" los velocity
  real<lower=0.> sigma_los;
  real<lower=0.> alpha;
  real<lower=0.> rho;
}

transformed parameters {
  real sin_i = sqrt(1.-cos_i^2);
  vector[N] xc = x - x_c;
  vector[N] yc = y - y_c;
  real sin_phi = sin(phi);
  real cos_phi = cos(phi);
  matrix[N, 2] X;
  matrix[N, 2] Xhat;
  matrix[2, 2] Rot;
  matrix[2, 2] Stretch;
  vector[N] r;
  vector[N] xhat;
  vector[N] yhat;
  row_vector[2] xyhat[N];
  
  X = append_col(xc, yc);
  Rot = [ [-cos_phi, -sin_phi],
          [-sin_phi,  cos_phi]];
  Stretch = diag_matrix([1./cos_i, 1.]');
  Xhat = X * Rot * Stretch;
  xhat = Xhat[ : , 1];
  yhat = Xhat[ : , 2];
  r = sqrt(xhat .* xhat + yhat .* yhat);
  
  xyhat[:, 1] = to_array_1d(xhat);
  xyhat[:, 2] = to_array_1d(yhat);
}

model {
  matrix[N, N] cov;
  matrix[N, N] L_cov;
  
  // priors
  
  phi ~ normal(phi0, sd_phi0);
  cos_i ~ normal(ci0, sd_ci0);
  x_c ~ normal(0, sd_kc0/r_norm);
  y_c ~ normal(0, sd_kc0/r_norm);
  v_sys ~ normal(0, 150./v_norm);
  c_r ~ normal(0, 500./v_norm);
  sigma_los ~ normal(0, 50./v_norm);
  alpha ~ normal(0, 150./v_norm);
  rho ~ inv_gamma(a, b);
  
  // likelihood
  
  v ~ normal(v_los, dv);
  
  cov = cov_exp_quad(xyhat, alpha, rho);
  for (n in 1:N) {
    cov[n, n] = cov[n, n] + square(sigma_los);
  }
  L_cov = cholesky_decompose(cov);
  v_los ~ multi_normal_cholesky(v_sys + sin_i * sump(c_r, r, order) .* yhat, L_cov);
}
generated quantities {
  vector[N] v_rot;
  vector[N_xy] v_res;
  vector[N_xy] v_pred;
  vector[N_r+1] vrot_pred;
  vector[N_r+1] vexp_pred;
  
  v_rot = fabs(r .* sump(c_r, r, order));
  v_res = gp_pred_rng(xy_pred, v_los-v_sys-sin_i*sump(c_r, r, order) .* yhat, xyhat, alpha, rho, sigma_los);
  v_pred = (v_res + sin_i * sump(c_r, r_pred, order) .* y_pred) *v_norm/sin_i;
  v_res = v_res * v_norm/sin_i;
  vrot_pred[1] = 0.;
  vexp_pred[1] = 0.;
  for (n in 1:N_r) {
    vrot_pred[n+1] = dot_product(v_pred[((n-1)*N_r+1):(n*N_r)], cos(theta_pred[((n-1)*N_r+1):(n*N_r)]));
    vexp_pred[n+1] = dot_product(v_pred[((n-1)*N_r+1):(n*N_r)], sin(theta_pred[((n-1)*N_r+1):(n*N_r)]));
  }
  vrot_pred = fabs(vrot_pred * 2./N_r);
  vexp_pred = vexp_pred * 2./N_r;
}

