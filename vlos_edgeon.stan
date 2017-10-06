functions {
  real sump(vector c, real r, int order) {
    real s = c[1];
    for (i in 2:order) {
      s = s + c[i]*r^(i-1);
    }
    return s;
  }
}

data {
  int<lower=1> order;  //order for polynomial representation of velocities
  int<lower=1> N;
  real x[N];
  real y[N];
  real v[N];  //measured los velocity
  real dv[N]; //error estimate on v
  real phi0;
  real<lower=0.> sd_phi0;
  int<lower=1> N_r; //number of points to fill in
  real r_post[N_r];
}

parameters {
  real phi;
  real x_c;  //kinematic centers
  real y_c;
  real v_los[N];  //latent "real" los velocity
  vector[order] c_rot;
  real<lower=0.> sigma_los;
}

transformed parameters {
  real xp[N];
  
  for (i in 1:N) {
    xp[i] = -(x[i] - x_c) * sin(phi) + (y[i] - y_c) * cos(phi);
  }
}

model {
  phi ~ normal(phi0, sd_phi0);
  x_c ~ normal(0, 2);
  y_c ~ normal(0, 2);
  sigma_los ~ normal(0, 20);
  for (i in 1:order) {
    c_rot[i] ~ normal(0, 100);
  }
  
  for (i in 1:N) {
    v[i] ~ normal(v_los[i], dv[i]);
    v_los[i] ~ normal(sump(c_rot, fabs(xp[i]), order) * xp[i], sigma_los);
  }                                           
}

generated quantities {
  real v_rot[N];
  real v_model[N];
  real v_res[N];
  real vrot_post[N_r];
  
  for (i in 1:N) {
    v_model[i] = xp[i]*sump(c_rot, fabs(xp[i]), order);
    v_rot[i] = fabs(v_model[i]);
    v_res[i] = v[i] - v_model[i];
  }
  for (i in 1:N_r) {
    vrot_post[i] = fabs(r_post[i] * sump(c_rot, r_post[i], order));
  }
}

