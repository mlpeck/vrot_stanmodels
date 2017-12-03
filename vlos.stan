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
  real<lower=0., upper=1.> si;  // sine disk inclination
  real x_c;  //kinematic centers
  real y_c;
  real v_sys;     //system velocity offset (should be small)
  vector[N] v_los;  //latent "real" los velocity
  vector[order] c_rot;
  vector[order] c_exp;
  real<lower=0.> sigma_los;
}

transformed parameters {
  vector[N] dx = x - x_c;
  vector[N] dy = y - y_c;
  vector[N] yhat = -dx * sin(phi) + dy * cos(phi);
  vector[N] xhat = -(dx * cos(phi) + dy * sin(phi))/sqrt(1-si^2);
  vector[N] r;

  r = sqrt(yhat .* yhat + xhat .* xhat);
}

model {
  phi ~ normal(phi0, sd_phi0);
  si ~ normal(si0, sd_si0);
  x_c ~ normal(0, sd_kc0/r_norm);
  y_c ~ normal(0, sd_kc0/r_norm);
  v_sys ~ normal(0, 150./v_norm);
  sigma_los ~ normal(0, 50./v_norm);
  c_rot ~ normal(0, 5);                                                                        
  c_exp ~ normal(0, 5);
  
  v ~ normal(v_los, dv);
  v_los ~ normal(v_sys + si * (
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
  v_model = v_sys + si * (
                      sump(c_rot, r, order) .* yhat + 
                      sump(c_exp, r, order) .* xhat);
  v_res = v - v_model;
  vrot_post = fabs(r_post .* sump(c_rot, r_post, order));
  for (i in 1:N) {
  	  log_lik[i] = normal_lpdf(v[i] | v_los[i], dv[i]) +
  	               normal_lpdf(v_los[i] | v_model[i], sigma_los);
  }
}

