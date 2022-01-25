##This is the code for the models for the analysis of foraging children's participation in play and work. 
##WARNING: These models take weeks to fit. 
##NOTE: This code will not run without requesting the Tsimane dataset from Jonathan Stieglitz (jonathan.stieglitz@iast.fr)
##Make sure to set working directory to 'data'
##IF you don't have the necessary packages, load them following the instructions here https://xcelab.net/rm/software/

library(rethinking) ##load rethinking package

##read csv file with complete data into workspace
d<-read.csv("dataset.csv")

##Dataset without Dukha
ds<-subset(d,Society!="Dukha")

##data
N<-nrow(d)
d$y<-coerce_index(d$variable)
K<-max(d$y)
d$id<-coerce_index(d$newID)
N_id<-max(d$id)
d$society<-coerce_index(d$Society)
N_society<-max(d$society)
d$NPP_z<-(d$NPP-mean(d$NPP))/sd(d$NPP)
d$nonforaged_z<-(d$nonforaged-mean(d$nonforaged))/sd(d$nonforaged)
d$temp_z<-(d$meanAnnualTemp-mean(d$meanAnnualTemp))/sd(d$meanAnnualTemp)
d$prec_z<-(d$totalAnnualPrec-mean(d$totalAnnualPrec))/sd(d$totalAnnualPrec)
d$CV_z<-(d$precSeasonality-mean(d$precSeasonality))/sd(d$precSeasonality)
d$snake_z<-(d$snake_count-mean(d$snake_count))/sd(d$snake_count)
d$dens<-ifelse(d$sum_density>1,1,0)
d$Middle<-ifelse(d$newAge<13&d$newAge>6,1,0)
d$Ado<-ifelse(d$newAge>12,1,0)

##Dataset without Dukha
N_ds<-nrow(ds)
ds$y<-coerce_index(ds$variable)
K_ds<-max(ds$y)
ds$id<-coerce_index(ds$newID)
N_id_ds<-max(ds$id)
ds$society<-coerce_index(ds$Society)
N_society_ds<-max(ds$society)
ds$NPP_z<-(ds$NPP-mean(ds$NPP))/sd(ds$NPP)
ds$nonforaged_z<-(ds$nonforaged-mean(ds$nonforaged))/sd(ds$nonforaged)
ds$temp_z<-(ds$meanAnnualTemp-mean(ds$meanAnnualTemp))/sd(ds$meanAnnualTemp)
ds$prec_z<-(ds$totalAnnualPrec-mean(ds$totalAnnualPrec))/sd(ds$totalAnnualPrec)
ds$Middle<-ifelse(ds$newAge<13&ds$newAge>6,1,0)
ds$Ado<-ifelse(ds$newAge>12,1,0)

##dataset without adolescents
dc<-subset(d,Age_cat!=3)
N_dc<-nrow(dc)
dc$y<-coerce_index(dc$variable)
K_dc<-max(dc$y)
dc$id<-coerce_index(dc$newID)
N_id_dc<-max(dc$id)
dc$society<-coerce_index(dc$Society)
N_society_dc<-max(dc$society)
dc$NPP_z<-(dc$NPP-mean(dc$NPP))/sd(dc$NPP)
dc$nonforaged_z<-(dc$nonforaged-mean(dc$nonforaged))/sd(dc$nonforaged)
dc$temp_z<-(dc$meanAnnualTemp-mean(dc$meanAnnualTemp))/sd(dc$meanAnnualTemp)
dc$prec_z<-(dc$totalAnnualPrec-mean(dc$totalAnnualPrec))/sd(dc$totalAnnualPrec)
dc$Middle<-ifelse(dc$newAge<13&dc$newAge>6,1,0)
dc$dens<-ifelse(dc$sum_density>1,1,0)

##specify data for model 1
data_list1<-list(
  K=K,
  N=N,
  N_id=N_id,
  N_society=N_society,
  y=d$y,
  id=d$id,
  society=d$society
)

model1 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] +v_society[society[i],k]+
v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] +v_society[society[i],k]+
v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model1 <- list(
  a = rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model1 <- 3
init_model1 <- list()
for ( i in 1:n_chains_model1 ) init_model1[[i]] <- start_model1

mfit_model1 <- cstan( model_code=model1 , data=data_list1 , chains=n_chains_model1 , cores= 3 , warmup=1000, iter=2000, init=init_model1 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA1<-precis(mfit_model1,dept=3,prob=0.95)
write.csv(TA1,"TA1.csv")
post1<-extract.samples(mfit_model1)
save(post1, file = "post1.rda")

##specify data for model 2
data_list2<-list(
  K=K,
  N=N,
  N_id=N_id,
  N_society=N_society,
  y=d$y,
  id=d$id,
  society=d$society,
  sex=d$Sex,
  middle=d$Middle,
  ado=d$Ado,
  nonforaged=d$nonforaged_z)

model2 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model2 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model2 <- 3
init_model2 <- list()
for ( i in 1:n_chains_model2 ) init_model2[[i]] <- start_model2

mfit_model2 <- cstan( model_code=model2 , data=data_list2 , chains=n_chains_model2 , cores= 3 , warmup=1000, iter=2000, init=init_model2 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA2<-precis(mfit_model2,dept=3,prob=0.95)
write.csv(TA2,"TA2.csv")
post2<-extract.samples(mfit_model2)
save(post2, file = "post2.rda")

##specify data for model 5, with division of labor 
data_list5<-list(
  K=K,
  N=N,
  N_id=N_id,
  N_society=N_society,
  y=d$y,
  id=d$id,
  society=d$society,
  sex=d$Sex,
  div=d$div,
  middle=d$Middle,
  ado=d$Ado,
  nonforaged=d$nonforaged_z)

model5 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real div[N];
real nonforaged[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_div[K-1];
real b_sex_div[K-1];
real b_nonforaged[K-1];
matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_div~ normal(0,1);
b_sex_div~ normal(0,1);
b_nonforaged~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_div[k] * div[i]
+ b_sex_div[k] * sex[i] * div[i]
+ b_nonforaged[k] * nonforaged[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_div[k] * div[i]
+ b_sex_div[k] * sex[i] * div[i]
+ b_nonforaged[k] * nonforaged[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model5 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_div= rep(0,K-1),
  b_sex_div= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model5 <- 3
init_model5 <- list()
for ( i in 1:n_chains_model5 ) init_model5[[i]] <- start_model5

mfit_model5 <- cstan( model_code=model5 , data=data_list5 , chains=n_chains_model5 , cores= 3 , warmup=1000, iter=2000, init=init_model5 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA5<-precis(mfit_model5,dept=3,prob=0.95)
write.csv(TA5,"TA5.csv")
post5<-extract.samples(mfit_model5)
save(post5, file = "post5.rda")

##specify data for model 3
data_list3<-list(
  K=K,
  N=N,
  N_id=N_id,
  N_society=N_society,
  y=d$y,
  id=d$id,
  society=d$society,
  sex=d$Sex,
  middle=d$Middle,
  ado=d$Ado,
  nonforaged=d$nonforaged_z,
  NPP=d$NPP_z,
  temp=d$temp_z,
  prec=d$prec_z,
  CV=d$CV_z)

model3 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real NPP[N];
real temp[N];
real prec[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_NPP[K-1];
real b_temp[K-1];
real b_prec[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_NPP~ normal(0,1);
b_temp~ normal(0,1);
b_prec~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ b_temp[k] * temp[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ b_temp[k] * temp[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_NPP= rep(0,K-1),
  b_temp= rep(0,K-1),
  b_prec= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model3 <- 3
init_model3 <- list()
for ( i in 1:n_chains_model3 ) init_model3[[i]] <- start_model3

mfit_model3 <- cstan( model_code=model3 , data=data_list3 , chains=n_chains_model3 , cores= 3 , warmup=1000, iter=2000, init=init_model3 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3<-precis(mfit_model3,dept=3,prob=0.95)
write.csv(TA3,"TA3.csv")
post3<-extract.samples(mfit_model3)
save(post3, file = "post3.rda")

##specify data for model 4
data_list4<-list(
  K=K,
  N=N,
  N_id=N_id,
  N_society=N_society,
  y=d$y,
  id=d$id,
  society=d$society,
  sex=d$Sex,
  middle=d$Middle,
  ado=d$Ado,
  nonforaged=d$nonforaged_z,
  water=d$water_rating,
  dens=d$dens,
  snake=d$snake_z)

model4 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
int water[N];
int dens[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_water[K-1];
real b_dens[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_water~ normal(0,1);
b_dens~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_water[k] * water[i]
+ b_dens[k] * dens[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_water[k] * water[i]
+ b_dens[k] * dens[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model4 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_water= rep(1,K-1),
  b_dens= rep(1,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model4 <- 3
init_model4 <- list()
for ( i in 1:n_chains_model4 ) init_model4[[i]] <- start_model4

mfit_model4 <- cstan( model_code=model4 , data=data_list4 , chains=n_chains_model4 , cores= 3 , warmup=1000, iter=2000, init=init_model4, control = list(adapt_delta = 0.99, max_treedepth = 15))

TA4<-precis(mfit_model4,dept=3,prob=0.95)
write.csv(TA4,"TA4.csv")
post4<-extract.samples(mfit_model4)
save(post4, file = "post4.rda")

##calculate WAIC
waicTA<-compare(mfit_model1,mfit_model2,mfit_model3, mfit_model4, mfit_model5)
plot(waicTA)

##Supplementary models
model3.1 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real NPP[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_NPP[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_NPP~ normal(0,1);


// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.1 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_NPP= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model3.1 <- 3
init_model3.1 <- list()
for ( i in 1:n_chains_model3.1 ) init_model3.1[[i]] <- start_model3.1

mfit_model3.1 <- cstan( model_code=model3.1 , data=data_list3 , chains=n_chains_model3.1, cores= 3 , warmup=1000, iter=2000, init=init_model3.1 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.1<-precis(mfit_model3.1,dept=3,prob=0.95)
write.csv(TA3.1,"TA3.1.csv")
post3.1<-extract.samples(mfit_model3.1)
save(post3.1, file = "post3.1.rda")

model3.2 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real temp[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_temp[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_temp~ normal(0,1);


// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_temp[k] * temp[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_temp[k] * temp[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.2 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_temp= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model3.2 <- 3
init_model3.2 <- list()
for ( i in 1:n_chains_model3.2 ) init_model3.2[[i]] <- start_model3.2

mfit_model3.2 <- cstan( model_code=model3.2 , data=data_list3 , chains=n_chains_model3.2, cores= 3 , warmup=1000, iter=2000, init=init_model3.2 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.2<-precis(mfit_model3.2,dept=3,prob=0.95)
write.csv(TA3.2,"TA3.2.csv")
post3.2<-extract.samples(mfit_model3.2)
save(post3.2, file = "post3.2.rda")

model3.3 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real prec[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_prec[K-1];


matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_prec~ normal(0,1);


// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.3 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_prec= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model3.3 <- 3
init_model3.3 <- list()
for ( i in 1:n_chains_model3.3 ) init_model3.3[[i]] <- start_model3.3

mfit_model3.3 <- cstan( model_code=model3.3 , data=data_list3 , chains=n_chains_model3.3, cores= 3 , warmup=1000, iter=2000, init=init_model3.3 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.3<-precis(mfit_model3.3,dept=3,prob=0.95)
write.csv(TA3.3,"TA3.3.csv")
post3.3<-extract.samples(mfit_model3.3)
save(post3.3, file = "post3.3.rda")

model3.4 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real CV[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_CV[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_CV~ normal(0,1);


// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_CV[k] * CV[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_CV[k] * CV[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.4 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_CV= rep(0,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model3.4 <- 3
init_model3.4 <- list()
for ( i in 1:n_chains_model3.4 ) init_model3.4[[i]] <- start_model3.4

mfit_model3.4 <- cstan( model_code=model3.4 , data=data_list3 , chains=n_chains_model3.4, cores= 3 , warmup=1000, iter=2000, init=init_model3.4 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.4<-precis(mfit_model3.4,dept=3,prob=0.95)
write.csv(TA3.4,"TA3.4.csv")
post3.4<-extract.samples(mfit_model3.4)
save(post3.4, file = "post3.4.rda")

data_list3.5<-list(
  K=K_ds,
  N=N_ds,
  N_id=N_id_ds,
  N_society=N_society_ds,
  y=ds$y,
  id=ds$id,
  society=ds$society,
  sex=ds$Sex,
  middle=ds$Middle,
  ado=ds$Ado,
  nonforaged=ds$nonforaged_z,
  NPP=ds$NPP_z,
  temp=ds$temp_z,
  prec=ds$prec_z)

model3.5 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real NPP[N];
real temp[N];
real prec[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_NPP[K-1];
real b_temp[K-1];
real b_prec[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_NPP~ normal(0,1);
b_temp~ normal(0,1);
b_prec~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ b_temp[k] * temp[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ b_temp[k] * temp[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.5 <- list(
  a = rep(0,K_ds-1),
  b_sex = rep(0,K_ds-1),
  b_middle= rep(0,K_ds-1),
  b_ado= rep(0,K_ds-1),
  b_sex_middle= rep(0,K_ds-1),
  b_sex_ado= rep(0,K_ds-1),
  b_nonforaged= rep(0,K_ds-1),
  b_NPP= rep(0,K_ds-1),
  b_temp= rep(0,K_ds-1),
  b_prec= rep(0,K_ds-1),
  sigma_id = rep(1,K_ds-1),
  L_Rho_id = diag(K_ds-1),
  z_id = matrix(0,nrow=K_ds-1,ncol=N_id_ds),
  sigma_society = rep(1,K_ds-1),
  L_Rho_society = diag(K_ds-1),
  z_society = matrix(0,nrow=K_ds-1,ncol=N_society_ds)
)

n_chains_model3.5 <- 3
init_model3.5 <- list()
for ( i in 1:n_chains_model3.5 ) init_model3.5[[i]] <- start_model3.5

mfit_model3.5 <- cstan( model_code=model3.5 , data=data_list3.5 , chains=n_chains_model3.5 , cores= 3 , warmup=1000, iter=2000, init=init_model3.5 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.5<-precis(mfit_model3.5,dept=3,prob=0.95)
write.csv(TA3.5,"TA3.5.csv")
post3.5<-extract.samples(mfit_model3.5)
save(post3.5, file = "post3.5.rda")

##specify data for model 3.6
data_list3.6<-list(
  K=K_dc,
  N=N_dc,
  N_id=N_id_dc,
  N_society=N_society_dc,
  y=dc$y,
  id=dc$id,
  society=dc$society,
  sex=dc$Sex,
  middle=dc$Middle,
  nonforaged=dc$nonforaged_z,
  NPP=dc$NPP_z,
  temp=dc$temp_z,
  prec=dc$prec_z)

model3.6 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
real nonforaged[N];
real NPP[N];
real temp[N];
real prec[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_sex_middle[K-1];
real b_nonforaged[K-1];
real b_NPP[K-1];
real b_temp[K-1];
real b_prec[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_sex_middle~ normal(0,1);
b_nonforaged~ normal(0,1);
b_NPP~ normal(0,1);
b_temp~ normal(0,1);
b_prec~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ b_temp[k] * temp[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_NPP[k] * NPP[i]
+ b_temp[k] * temp[i]
+ b_prec[k] * prec[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.6 <- list(
  a = rep(0,K_dc-1),
  b_sex = rep(0,K_dc-1),
  b_middle= rep(0,K_dc-1),
  b_sex_middle= rep(0,K_dc-1),
  b_nonforaged= rep(0,K_dc-1),
  b_NPP= rep(0,K_dc-1),
  b_temp= rep(0,K_dc-1),
  b_prec= rep(0,K_dc-1),
  sigma_id = rep(1,K_dc-1),
  L_Rho_id = diag(K_dc-1),
  z_id = matrix(0,nrow=K_dc-1,ncol=N_id_dc),
  sigma_society = rep(1,K_dc-1),
  L_Rho_society = diag(K_dc-1),
  z_society = matrix(0,nrow=K_dc-1,ncol=N_society_dc)
)

n_chains_model3.6 <- 3
init_model3.6 <- list()
for ( i in 1:n_chains_model3.6 ) init_model3.6[[i]] <- start_model3.6

mfit_model3.6 <- cstan( model_code=model3.6 , data=data_list3.6 , chains=n_chains_model3.6 , cores= 3 , warmup=1000, iter=2000, init=init_model3.6 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.6<-precis(mfit_model3.6,dept=3,prob=0.95)
write.csv(TA3.6,"TA3.6.csv")
post3.6<-extract.samples(mfit_model3.6)
save(post3.6, file = "post3.6.rda")

model3.7 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real temp[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_temp[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_temp~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_temp[k] * temp[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_temp[k] * temp[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model3.7 <- list(
  a = rep(0,K_ds-1),
  b_sex = rep(0,K_ds-1),
  b_middle= rep(0,K_ds-1),
  b_ado= rep(0,K_ds-1),
  b_sex_middle= rep(0,K_ds-1),
  b_sex_ado= rep(0,K_ds-1),
  b_nonforaged= rep(0,K_ds-1),
  b_temp= rep(0,K_ds-1),
  sigma_id = rep(1,K_ds-1),
  L_Rho_id = diag(K_ds-1),
  z_id = matrix(0,nrow=K_ds-1,ncol=N_id_ds),
  sigma_society = rep(1,K_ds-1),
  L_Rho_society = diag(K_ds-1),
  z_society = matrix(0,nrow=K_ds-1,ncol=N_society_ds)
)

n_chains_model3.7 <- 3
init_model3.7 <- list()
for ( i in 1:n_chains_model3.7 ) init_model3.7[[i]] <- start_model3.7

mfit_model3.7 <- cstan( model_code=model3.7 , data=data_list3.5 , chains=n_chains_model3.7 , cores= 3 , warmup=1000, iter=2000, init=init_model3.7 , control = list(adapt_delta = 0.99, max_treedepth = 15))

TA3.7<-precis(mfit_model3.7,dept=3,prob=0.95)
write.csv(TA3.7,"TA3.7.csv")
post3.7<-extract.samples(mfit_model3.7)
save(post3.7, file = "post3.7.rda")

model4.1 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
int water[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_water[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_water~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_water[k] * water[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_water[k] * water[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model4.1 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_water= rep(1,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model4.1 <- 3
init_model4.1 <- list()
for ( i in 1:n_chains_model4.1 ) init_model4.1[[i]] <- start_model4.1

mfit_model4.1 <- cstan( model_code=model4.1 , data=data_list4 , chains=n_chains_model4.1 , cores= 3 , warmup=1000, iter=2000, init=init_model4.1, control = list(adapt_delta = 0.99, max_treedepth = 15))

TA4.1<-precis(mfit_model4.1,dept=3,prob=0.95)
write.csv(TA4.1,"TA4.1.csv")
post4.1<-extract.samples(mfit_model4.1)
save(post4.1, file = "post4.1.rda")

model4.2 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
int dens[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_dens[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_dens~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_dens[k] * dens[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_dens[k] * dens[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model4.2 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_dens= rep(1,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model4.2 <- 3
init_model4.2 <- list()
for ( i in 1:n_chains_model4.2 ) init_model4.2[[i]] <- start_model4.2

mfit_model4.2 <- cstan( model_code=model4.2 , data=data_list4 , chains=n_chains_model4.2 , cores= 3 , warmup=1000, iter=2000, init=init_model4.2, control = list(adapt_delta = 0.99, max_treedepth = 15))

TA4.2<-precis(mfit_model4.2,dept=3,prob=0.95)
write.csv(TA4.2,"TA4.2.csv")
post4.2<-extract.samples(mfit_model4.2)
save(post4.2, file = "post4.2.rda")

model4.3 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
int ado[N];
real nonforaged[N];
real snake[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_ado[K-1];
real b_sex_middle[K-1];
real b_sex_ado[K-1];
real b_nonforaged[K-1];
real b_snake[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_ado~ normal(0,1);
b_sex_middle~ normal(0,1);
b_sex_ado~ normal(0,1);
b_nonforaged~ normal(0,1);
b_snake~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_snake[k] * snake[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_ado[k] * ado[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_sex_ado[k] * sex[i] * ado[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_snake[k] * snake[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model4.3 <- list(
  a = rep(0,K-1),
  b_sex = rep(0,K-1),
  b_middle= rep(0,K-1),
  b_ado= rep(0,K-1),
  b_sex_middle= rep(0,K-1),
  b_sex_ado= rep(0,K-1),
  b_nonforaged= rep(0,K-1),
  b_snake= rep(1,K-1),
  sigma_id = rep(1,K-1),
  L_Rho_id = diag(K-1),
  z_id = matrix(0,nrow=K-1,ncol=N_id),
  sigma_society = rep(1,K-1),
  L_Rho_society = diag(K-1),
  z_society = matrix(0,nrow=K-1,ncol=N_society)
)

n_chains_model4.3 <- 3
init_model4.3 <- list()
for ( i in 1:n_chains_model4.3 ) init_model4.3[[i]] <- start_model4.3

mfit_model4.3 <- cstan( model_code=model4.3 , data=data_list4 , chains=n_chains_model4.3 , cores= 3 , warmup=1000, iter=2000, init=init_model4.3, control = list(adapt_delta = 0.99, max_treedepth = 15))

TA4.3<-precis(mfit_model4.3,dept=3,prob=0.95)
write.csv(TA4.3,"TA4.3.csv")
post4.3<-extract.samples(mfit_model4.3)
save(post4.3, file = "post4.3.rda")

data_list4.4<-list(
  K=K_dc,
  N=N_dc,
  N_id=N_id_dc,
  N_society=N_society_dc,
  y=dc$y,
  id=dc$id,
  society=dc$society,
  sex=dc$Sex,
  middle=dc$Middle,
  nonforaged=dc$nonforaged_z,
  water=dc$water_rating,
  dens=dc$dens)

model4.4 <- "
data{
int N;
int N_id;
int N_society;
int y[N];
int id[N];
int society[N];
int sex[N];
int middle[N];
real nonforaged[N];
int water[N];
int dens[N];
int K;
}
parameters{
real a[K-1];                // intercepts for each behavior
real b_sex[K-1];
real b_middle[K-1];
real b_sex_middle[K-1];
real b_nonforaged[K-1];
real b_water[K-1];
real b_dens[K-1];

matrix[K-1,N_id] z_id;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_id;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_id;         
matrix[K-1,N_society] z_society;      // matrix of standardized random effects
vector<lower=0>[K-1] sigma_society;   // stddev of random effects
cholesky_factor_corr[K-1] L_Rho_society;        
}
transformed parameters{
matrix[N_id,K-1] v_id;      // matrix of scaled random effects
matrix[N_society,K-1] v_society;      // matrix of scaled random effects
v_id = (diag_pre_multiply(sigma_id,L_Rho_id) * z_id)';
v_society = (diag_pre_multiply(sigma_society,L_Rho_society) * z_society)';
}
model{

// priors
a ~ normal(0,1);
b_sex~ normal(0,1);
b_middle~ normal(0,1);
b_sex_middle~ normal(0,1);
b_nonforaged~ normal(0,1);
b_water~ normal(0,1);
b_dens~ normal(0,1);

// hyper-prior
to_vector(z_id) ~ normal(0,1);
sigma_id ~ exponential(1);
L_Rho_id ~ lkj_corr_cholesky(2);

to_vector(z_society) ~ normal(0,1);
sigma_society ~ exponential(1);
L_Rho_society ~ lkj_corr_cholesky(2);

// likelihood
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_water[k] * water[i]
+ b_dens[k] * dens[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
y[i] ~ categorical_logit( p );
}
}
generated quantities{
matrix[K-1,K-1] Rho_id;
matrix[K-1,K-1] Rho_society;
vector[N] log_lik;
Rho_id = L_Rho_id * L_Rho_id';
Rho_society = L_Rho_society * L_Rho_society';
for ( i in 1:N ) {
vector[K] p;
for ( k in 1:(K-1) ) 
p[k] = a[k] 
+ b_sex[k] * sex[i]
+ b_middle[k] * middle[i]
+ b_sex_middle[k] * sex[i] * middle[i]
+ b_nonforaged[k] * nonforaged[i]
+ b_water[k] * water[i]
+ b_dens[k] * dens[i]
+ v_society[society[i],k]
+ v_id[id[i],k];
p[K] = 0;
log_lik[i] = categorical_logit_lpmf( y[i] | p );
}
}
"

start_model4.4 <- list(
  a = rep(0,K_dc-1),
  b_sex = rep(0,K_dc-1),
  b_middle= rep(0,K_dc-1),
  b_sex_middle= rep(0,K_dc-1),
  b_nonforaged= rep(0,K_dc-1),
  b_water= rep(1,K_dc-1),
  b_dens= rep(1,K_dc-1),
  sigma_id = rep(1,K_dc-1),
  L_Rho_id = diag(K_dc-1),
  z_id = matrix(0,nrow=K_dc-1,ncol=N_id_dc),
  sigma_society = rep(1,K_dc-1),
  L_Rho_society = diag(K_dc-1),
  z_society = matrix(0,nrow=K_dc-1,ncol=N_society_dc)
)

n_chains_model4.4 <- 3
init_model4.4 <- list()
for ( i in 1:n_chains_model4.4 ) init_model4.4[[i]] <- start_model4.4

mfit_model4.4 <- cstan( model_code=model4.4 , data=data_list4.4 , chains=n_chains_model4.4 , cores= 3 , warmup=1000, iter=2000, init=init_model4.4, control = list(adapt_delta = 0.99, max_treedepth = 15))

TA4.4<-precis(mfit_model4.4,dept=3,prob=0.95)
write.csv(TA4.4,"TA4.4.csv")
post4.4<-extract.samples(mfit_model4.4)
save(post4.4, file = "post4.4.rda")