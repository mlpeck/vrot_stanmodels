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
  int<lower=1> N_r; //number of points to fill in
  vector[N_r] r_post;
}

parameters {
  real phi;
  real<lower=0., upper=1.> si;  //sine disk inclination
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
  vector[N] r;

  for (i in 1:N) {
    r[i] = sqrt((-dx[i] * sin(phi) + dy[i] * cos(phi))^2 +
                ((dx[i] * cos(phi) + dy[i] * sin(phi))/sqrt(1-si^2))^2);
  }
  
}

model {
  phi ~ normal(phi0, sd_phi0);
  si ~ normal(si0, sd_si0);
  x_c ~ normal(0, 2);
  y_c ~ normal(0, 2);
  v_sys ~ normal(0, 50);
  sigma_los ~ normal(0, 20);
  c_rot ~ normal(0, 100);
  c_exp ~ normal(0, 100);
  
  v ~ normal(v_los, dv);
  v_los ~ normal(v_sys + si * (
                  sump(c_rot, r, order) .* (-dx * sin(phi) + dy * cos(phi))  + 
                  sump(c_exp, r, order) .* (-dx * cos(phi) - dy * sin(phi)) / sqrt(1-si^2)),
                  sigma_los);
}

generated quantities {
  vector[N] v_rot;
  vector[N] v_exp;
  vector[N] v_model;
  vector[N] v_res;
  vector[N_r] vrot_post;
  
  v_rot = fabs(r .* sump(c_rot, r, order));
  v_exp = r .* sump(c_exp, r, order);
  v_model = v_sys + si * (
                      sump(c_rot, r, order) .* (-dx * sin(phi) + dy * cos(phi))  + 
                      sump(c_exp, r, order) .* (-dx * cos(phi) - dy * sin(phi)) / sqrt(1-si^2));
  v_res = v - v_model;
  vrot_post = fabs(r_post .* sump(c_rot, r_post, order));
}

