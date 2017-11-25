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
  int<lower=1> N_r; //number of points to fill in
  vector[N_r] r_post;
}

parameters {
  real phi;
  real x_c;  //kinematic centers
  real y_c;
  real v_sys;
  real v_los[N];  //latent "real" los velocity
  vector[order] c_rot;
  real<lower=0.> sigma_los;
}

transformed parameters {
  vector[N] r;
  
  r = -(x - x_c) * sin(phi) + (y - y_c) * cos(phi);
}

model {
  phi ~ normal(phi0, sd_phi0);
  x_c ~ normal(0, 2);
  y_c ~ normal(0, 2);
  v_sys ~ normal(0, 50);
  sigma_los ~ normal(0, 20);
  c_rot ~ normal(0, 100);
  
  v ~ normal(v_los, dv);
  v_los ~ normal(v_sys + sump(c_rot, fabs(r), order) .* r, sigma_los);
}

generated quantities {
  vector[N] v_rot;
  vector[N] v_exp;   //dummy variable so I can use a single function
  vector[N] v_model;
  vector[N] v_res;
  vector[N_r] vrot_post;
 
  v_model = v_sys + r .* sump(c_rot, fabs(r), order);
  v_rot = fabs(r .* sump(c_rot, fabs(r), order));
  v_exp = v_rot;
  v_res = v - v_model;

  vrot_post = fabs(r_post .* sump(c_rot, r_post, order));
}

