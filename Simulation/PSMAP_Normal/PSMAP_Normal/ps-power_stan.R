powerps<-"
//
//  STRATIFIED POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter As follows 
//
data {
  int<lower = 2> S;
  
  //existing data
  int<lower = 1>  N0[S];
  real            YBAR0[S];
  real<lower = 0> SD0[S];
  
  //current data
  int<lower = 1> N1[S];
  int<lower = 1> TN1;
  real           Y1[TN1];
  int<lower = 1> INX1[TN1];
  
  //prior of vs
  vector<lower=0>[S] RS;
  //fix vs
  int<lower = 0, upper = 1> FIXVS;
  
  //target borrowing
  real<lower = 0> A;
}

transformed data {
  row_vector<lower = 0, upper = 1>[S] WS1;
  for (i in 1:S) {
    WS1[i] = N1[i];
    WS1[i] = WS1[i]/TN1;
  }
}

parameters {
  simplex[S]    vs;
  vector[S]     thetas;
  real<lower=0> taus[S];
}

transformed parameters {
  real<lower = 0, upper = 1> as[S];
  real<lower = 0> sds[S];
  
  for (i in 1:S) {
    if (0 == FIXVS) {
      as[i]  = 1 < A*vs[i]/N0[i] ? 1:A*vs[i]/N0[i];
    } else {
      as[i]  = 1 < A*RS[i]/N0[i] ? 1:A*RS[i]/N0[i];
    }
    
    sds[i] = 0 == as[i] ? 0 : SD0[i] / sqrt(as[i] * N0[i]);
  }
}

model {
  //prior
  if (A > 0) {
    target += normal_lpdf(YBAR0 | thetas, sds);
  } else {
    thetas ~ normal(0, 1000);
  }
  vs      ~ dirichlet(RS);
  taus    ~ cauchy(0, 2.5);
  
  //likelihood
  Y1 ~ normal(thetas[INX1], taus[INX1]);
}

generated quantities {
  real theta;
  theta = WS1 * thetas;
} 
"







powerp<-"//
//  TYPICAL POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter A is fixed
//

data {
  //existing data
  int<lower = 1>  N0;
  real            YBAR0;
  real<lower = 0> SD0;

  //current data
  int<lower = 1>  TN1;
  real            Y1[TN1];

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  real<lower = 0> a0;
  real<lower = 0> sn0;

  a0  = 1 < A/N0 ? 1 : A/N0;
  sn0 = SD0/sqrt(N0*1.0);
}

parameters {
  real          theta;
  real<lower=0> tau1;
}

model {
  //prior
  theta ~ normal(0, 1000);
  tau1  ~ cauchy(0, 2.5);

  //likelihood
  if (N0 > 0) {
    target +=  normal_lpdf(YBAR0 | theta, sn0) * a0;
  }

  Y1 ~ normal(theta, tau1);
}
"




powerpsbinary<-"//
//  STRATIFIED POWER PRIOR FOR BINARY DATA
//  Power parameter As follows dirichlet
//
data {
  int<lower = 2> S;

  //existing data
  int<lower = 1> N0[S];
  real<lower = 0, upper = 1> YBAR0[S];

  //current data
  int<lower = 1> N1[S];
  int<lower = 0> YSUM1[S];

  //prior of vs
  vector<lower=0>[S] RS;

  //fix vs
  int<lower = 0, upper = 1> FIXVS;

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  row_vector<lower = 0, upper = 1>[S] WS1;
  int<lower = 0> sn1;

  sn1 = sum(N1);
  for (i in 1:S) {
    WS1[i] = N1[i];
    WS1[i] = WS1[i]/sn1;
  }
}

parameters {
  simplex[S] vs;
  vector<lower=0, upper=1>[S]  thetas;
}

transformed parameters {
  real<lower = 0, upper = 1> as[S];
  real<lower = 0> alphas[S];
  real<lower = 0> betas[S];

  for (i in 1:S) {
    if (0 == FIXVS) {
      as[i]  = 1 < A*vs[i]/N0[i] ? 1:A*vs[i]/N0[i];
    } else {
      as[i]  = 1 < A*RS[i]/N0[i] ? 1:A*RS[i]/N0[i];
    }
    alphas[i] = as[i] * N0[i] * YBAR0[i]  + 1;
    betas[i]  = as[i] * N0[i] * (1-YBAR0[i]) + 1;
  }
}

model {
  //prior
  if (A > 0) {
    target += beta_lpdf(thetas | alphas, betas);
  } else {
    thetas ~ uniform(0,1);
  }
  vs ~ dirichlet(RS);

  //likelihood
  YSUM1 ~ binomial(N1, thetas);
}

generated quantities {
  real theta;
  theta = WS1 * thetas;
}
"


prior<-"//
//  POWER PRIOR ONLY FOR CONTINUOUS DATA
//  Power parameter A is fixed
//

data {
  //existing data
  int<lower = 1>  N0;
  real            YBAR0;
  real<lower = 0> SD0;

  //target borrowing
  real<lower = 0> ALPHA;
}

transformed data {
  real<lower = 0> sn0;
  sn0 = SD0/sqrt(N0*1.0);
}

parameters {
  real theta;
}

model {
  //prior
  theta ~ normal(0, 1000);
  target +=  normal_lpdf(YBAR0 | theta, sn0) * ALPHA;
}
"



