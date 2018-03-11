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
}

data {
  int<lower=1> order;  //order for polynomial representation of velocities
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
  vector[N] v;  //measured los velocity
  vector[N] dv; //error estimate on v
  real phi0;
  real<lower=0.> sd_phi0;
  real si0;
  real<lower=0.> sd_si0;
  real<lower=0.> sd_kc0;
  real<lower=0.> r_norm;
  real<lower=0.> v_norm;
  int<lower=1> N_r; //number of points to fill in
  vector[N_r] r_post;
}

parameters {
  real phi;
  real<lower=0., upper=1.> sin_i;  // sine disk inclination
  real x_c;  //kinematic centers
  real y_c;
  real v_sys;     //system velocity offset (should be small)
  vector[N] v_los;  //latent "real" los velocity
  vector[order] c_rot;
  vector[order] c_exp;
  real<lower=0.> sigma_los;
}

transformed parameters {
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
  
  X = append_col(xc, yc);
  Rot = [ [-cos_phi, -sin_phi],
          [-sin_phi, cos_phi]];
  Stretch = diag_matrix([1./sqrt(1.-sin_i^2), 1.]');
  Xhat = X * Rot * Stretch;
  xhat = Xhat[ : , 1];
  yhat = Xhat[ : , 2];
  r = sqrt(yhat .* yhat + xhat .* xhat);
}

model {
  phi ~ normal(phi0, sd_phi0);
  sin_i ~ normal(si0, sd_si0);
  x_c ~ normal(0, sd_kc0/r_norm);
  y_c ~ normal(0, sd_kc0/r_norm);
  v_sys ~ normal(0, 150./v_norm);
  sigma_los ~ normal(0, 50./v_norm);
  c_rot ~ normal(0, 1000./v_norm);                                                                        
  c_exp ~ normal(0, 1000./v_norm);
  
  v ~ normal(v_los, dv);
  v_los ~ normal(v_sys + sin_i * (
                  sump(c_rot, r, order) .* yhat + sump(c_exp, r, order) .* xhat),
                  sigma_los);
}

generated quantities {
  vector[N] v_rot;
  vector[N] v_exp;
  vector[N] v_model;
  vector[N] v_res;
  vector[N] log_lik;
  vector[N_r] vrot_post;
  
  v_rot = fabs(r .* sump(c_rot, r, order));
  v_exp = r .* sump(c_exp, r, order);
  v_model = v_sys + sin_i * (
                      sump(c_rot, r, order) .* yhat + 
                      sump(c_exp, r, order) .* xhat);
  v_res = v - v_model;
  vrot_post = fabs(r_post .* sump(c_rot, r_post, order));
  for (i in 1:N) {
  	  log_lik[i] = normal_lpdf(v[i] | v_los[i], dv[i]) +
  	               normal_lpdf(v_los[i] | v_model[i], sigma_los);
  }
}

