library(R2jags)
#library(runjags)
#library(MCMCpack)
library(rjags)
library(coda)
library(bayesplot)
library(distr)
library(mcmcplots)
library(dplyr)
library(tidyr)
library(ggplot2)


beta_cons<- beta(2.5,2.5)
c_alpha<- sd(qnorm(qbeta(1:10000/10001,2.5,2.5)))


###########################################
#####HuProt Generating=Model 1############
##########################################

#generate (5 arrays): 

model_1_Y_HP<- matrix(NA, nrow=34610,ncol=5)

S_1<- matrix(NA,15385,5)
mu_ctrl<- matrix(NA,17,5)
tau_active<- rep(NA,5)
mu_active<- rep(NA,5)
sigma_1<- rep(NA,5)
error<- matrix(NA, 34610,5)
signal<- matrix(NA,34610,5)

for(i in 1:5){
  sigma_1[i]<- runif(1,0,10)
  error[,i]<- sigma_1[i]*(1/c_alpha)*qnorm(rbeta(34610,2.5,2.5))
  mu_active[i]<- runif(1,-10^2,10^2)
  tau_active[i]<- runif(1,0,10^3)
  mu_ctrl[,i]<- runif(17,-10^2,10^2)
  S_1[,i]<- rnorm(15385,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],48),rep(mu_ctrl[2,i],48),rep(mu_ctrl[3,i],48),rep(mu_ctrl[4,i],48),
                 rep(mu_ctrl[5,i],48),rep(mu_ctrl[6,i],48),rep(mu_ctrl[7,i],48),rep(mu_ctrl[8,i],48),
                 rep(mu_ctrl[9,i],48),rep(mu_ctrl[10,i],48),rep(mu_ctrl[11,i],48),rep(mu_ctrl[12,i],48),
                 rep(mu_ctrl[13,i],48),rep(mu_ctrl[14,i],48),rep(mu_ctrl[15,i],48),rep(mu_ctrl[16,i],48),
                 rep(mu_ctrl[17,i],3072),rep(S_1[,i],each=2))
  model_1_Y_HP[,i]<- signal[,i]+error[,i]
}


model_2_Y_HP<- matrix(NA, nrow=34610,ncol=5)

S_2<- matrix(NA,15385,5)
mu_ctrl<- matrix(NA,17,5)
tau_active<- rep(NA,5)
mu_active<- rep(NA,5)
sigma_2<- rep(NA,5)
error<- matrix(NA, 34610,5)
signal<- matrix(NA,34610,5)

for(i in 1:5){
  sigma_2[i]<- runif(1,0,10)
  error[,i]<- rnorm(34610,0,sigma_2[i])
  mu_active[i]<- runif(1,-10^2,10^2)
  tau_active[i]<- runif(1,0,10^3)
  mu_ctrl[,i]<- runif(17,-10^2,10^2)
  S_2[,i]<- rnorm(15385,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],48),rep(mu_ctrl[2,i],48),rep(mu_ctrl[3,i],48),rep(mu_ctrl[4,i],48),
                 rep(mu_ctrl[5,i],48),rep(mu_ctrl[6,i],48),rep(mu_ctrl[7,i],48),rep(mu_ctrl[8,i],48),
                 rep(mu_ctrl[9,i],48),rep(mu_ctrl[10,i],48),rep(mu_ctrl[11,i],48),rep(mu_ctrl[12,i],48),
                 rep(mu_ctrl[13,i],48),rep(mu_ctrl[14,i],48),rep(mu_ctrl[15,i],48),rep(mu_ctrl[16,i],48),
                 rep(mu_ctrl[17,i],3072),rep(S_2[,i], each=2))
  model_2_Y_HP[,i]<- signal[,i]+error[,i]
}


model_3_Y_HP<- matrix(NA, nrow=34610,ncol=5)

S_3<- matrix(NA,15385,5)
mu_ctrl<- matrix(NA,17,5)
tau_active<- rep(NA,5)
mu_active<- rep(NA,5)
sigma_3<- rep(NA,5)
error<- matrix(NA, 34610,5)
signal<- matrix(NA,34610,5)

for(i in 1:5){
  sigma_3[i]<- runif(1,0,10^3)
  error[,i]<- rnorm(34610,0,sigma_3[i])
  mu_active[i]<- runif(1,-10^2,10^2)
  tau_active[i]<- runif(1,0,10^3)
  mu_ctrl[,i]<- runif(17,-10^2,10^2)
  S_3[,i]<- rnorm(15385,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],48),rep(mu_ctrl[2,i],48),rep(mu_ctrl[3,i],48),rep(mu_ctrl[4,i],48),
                 rep(mu_ctrl[5,i],48),rep(mu_ctrl[6,i],48),rep(mu_ctrl[7,i],48),rep(mu_ctrl[8,i],48),
                 rep(mu_ctrl[9,i],48),rep(mu_ctrl[10,i],48),rep(mu_ctrl[11,i],48),rep(mu_ctrl[12,i],48),
                 rep(mu_ctrl[13,i],48),rep(mu_ctrl[14,i],48),rep(mu_ctrl[15,i],48),rep(mu_ctrl[16,i],48),
                 rep(mu_ctrl[17,i],3072),rep(S_3[,i], each=2))
  model_3_Y_HP[,i]<- signal[,i]+error[,i]
}


write.csv(S_1,"S_1_huprot.csv")
write.csv(S_2,"S_2_huprot.csv")
write.csv(S_3,"S_3_huprot.csv")

I_h<-34610
K_h<-15402
ones_h<-rep(1,I_h)

type_H<- c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48),rep(10,48),
           rep(11,48),rep(12,48),rep(13,48),rep(14,48),rep(15,48),rep(16,48),rep(17,3072),rep(18:15402,each=2))



get_huprot_data_list_alpha25<- function(r,Y,type,beta,c,I,K,ones){
  df_r_H<- as.data.frame(cbind(Y,type))
  colnames(df_r_H)<- c("Y", "Type_1")
  df_r_H<- df_r_H[order(df_r_H$Type_1),c("Y", "Type_1")]
  df_r_H$Y<- as.numeric(df_r_H$Y)
  df_r_H$Type_1<- as.numeric(df_r_H$Type_1)
  
  return(list('Y'=df_r_H$Y,"I"=I, "K"=K,"type_1"=df_r_H$Type_1,
              'ones'=ones,'beta_cons'=beta,'c_alpha'=c))
}


get_huprot_data_list_alpha1<- function(r,Y,type,I,K){
  df_r_H<- as.data.frame(cbind(Y,type))
  colnames(df_r_H)<- c("Y", "Type_1")
  df_r_H<- df_r_H[order(df_r_H$Type_1),]
  df_r_H$Y<- as.numeric(df_r_H$Y)
  df_r_H$Type_1<- as.numeric(df_r_H$Type_1)
  
  return(list('Y'=df_r_H$Y,"I"=I, "K"=K,"type_1"=df_r_H$Type_1))
}

get_huprot_data_list_hatSigma<- function(r,Y,type,I,K,sigma){
  df_r_H<- as.data.frame(cbind(Y,type))
  colnames(df_r_H)<- c("Y", "Type_1")
  df_r_H<- df_r_H[order(df_r_H$Type_1),]
  df_r_H$Y<- as.numeric(df_r_H$Y)
  df_r_H$Type_1<- as.numeric(df_r_H$Type_1)
  
  return(list('Y'=df_r_H$Y,"I"=I, "K"=K,"type_1"=df_r_H$Type_1, 'sigma'=sigma))
}


########
#Models
########


model_alpha_25<- "
model{
C<-10^5
for(i in 1:I){
L[i]<- (c_alpha/(sigma*beta_cons))*
((pnorm((Y[i]-S[type_1[i]])*c_alpha/sigma,0,1))^(1.5))*
((1-pnorm((Y[i]-S[type_1[i]])*c_alpha/sigma,0,1))^(1.5))*
((dnorm((Y[i]-S[type_1[i]])*c_alpha/sigma,0,1))^(1.5))
p[i]=L[i]/C
ones[i]~dbern(p[i])
}

for(j in 1:17){
S[j]~dunif(-10^2,10^2)
}

for(k in 18:K){
S[k]~dnorm(mu_s,sigma_s)
}


mu_s~dunif(-10^2,10^2)

sigma_s=1/(tau_s*tau_s)
tau_s~dunif(0,10^3)
sigma~dunif(0,10^3)
}
"


model_alpha_1<- "
model{
C<-10^5
for(i in 1:I){
Y[i]~dnorm(S[type_1[i]],1/(sigma*sigma))
}

for(j in 1:17){
S[j]~dunif(-10^2,10^2)
}

for(k in 18:K){
S[k]~dnorm(mu_s,sigma_s)
}


mu_s~dunif(-10^2,10^2)

sigma_s=1/(tau_s*tau_s)
tau_s~dunif(0,10^3)
sigma~dunif(0,10^3)
}
"

model_hatSigma<- "
model{
C<-10^5
for(i in 1:I){
Y[i]~dnorm(S[type_1[i]],1/(sigma*sigma))
}

for(j in 1:17){
S[j]~dunif(-10^2,10^2)
}

for(k in 18:K){
S[k]~dnorm(mu_s,sigma_s)
}


mu_s~dunif(-10^2,10^2)

sigma_s=1/(tau_s*tau_s)
tau_s~dunif(0,10^3)

}
"

get_percentile_residuals_huprot<- function(posterior, r, df_rY,type){
  percentile<- ecdf(c(as.vector(as.matrix(posterior[[1]][,18:15402])),as.vector(as.matrix(posterior[[2]][,18:15402])),as.vector(as.matrix(posterior[[3]][,18:15402]))))
  
  df_Y<- as.data.frame(cbind(df_rY,
                      type))
  colnames(df_Y)<- c("Y_vals","type")
  df_Y$Y_vals<- as.numeric(df_Y$Y_vals)
  df_Y$type<- as.factor(df_Y$type)
  
  Y_avg<- df_Y%>%
    group_by(type)%>%
    summarize(meanY=mean(Y_vals))
  
  percentile_residuals_r<- rep(NA,15385)
  for(l in 1:15385){
    percentile_residuals_r[l]<- percentile(Y_avg$meanY[l+17])
  }
  percentile_residuals_r_transformed<- qnorm(percentile_residuals_r)
  return(percentile_residuals_r_transformed)
}


get_posteriors<- function(model, data_list, updating_params,r){
  model_1 <-textConnection(model)
  model_initial<- jags.model(model_1,data=data_list, n.chains=3,
                             n.adapt=1000,quiet=FALSE)
  
  update(model_initial,5000)
  
  model_coda_r <- coda.samples(model_initial,updating_params,n.iter=10000)
  
  return(model_coda_r)
}

##extract posteriors with analysis model 1 data generating D

for(p in 1:5){
  huprot_mcmc_data_p_alpha25<- get_huprot_data_list_alpha25(r=p,Y=model_1_Y_HP[,p],type=type_H,
                                                            beta=beta_cons,c=c_alpha,I=34610,K=15402,ones = ones_h)
  huprot_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                               data_list = huprot_mcmc_data_p_alpha25,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha25<- as.data.frame(rbind(huprot_posteriors_p_alpha25[[1]],
                                              huprot_posteriors_p_alpha25[[2]],
                                              huprot_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha25.csv")
  write.csv(huprot_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha25,r=p,df_rY = huprot_mcmc_data_p_alpha25[[1]],type=huprot_mcmc_data_p_alpha25[[4]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha25.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}

##extract posteriors with analysis model 1 data generating E
for(p in 1:5){
  huprot_mcmc_data_p_alpha25<- get_huprot_data_list_alpha25(r=p,Y=model_2_Y_HP[,p],type=type_H,
                                                            beta=beta_cons,c=c_alpha,I=34610,K=15402,ones = ones_h)
  huprot_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                               data_list = huprot_mcmc_data_p_alpha25,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha25<- as.data.frame(rbind(huprot_posteriors_p_alpha25[[1]],
                                              huprot_posteriors_p_alpha25[[2]],
                                              huprot_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha25Model2.csv")
  write.csv(huprot_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha25,r=p,df_rY = huprot_mcmc_data_p_alpha25[[1]],type=huprot_mcmc_data_p_alpha25[[4]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha25Model2.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}

##extract posteriors with analysis model 1 data generating F

for(p in 1:5){
  huprot_mcmc_data_p_alpha25<- get_huprot_data_list_alpha25(r=p,Y=model_3_Y_HP[,p],type=type_H,
                                                            beta=beta_cons,c=c_alpha,I=34610,K=15402,ones = ones_h)
  huprot_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                               data_list = huprot_mcmc_data_p_alpha25,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha25<- as.data.frame(rbind(huprot_posteriors_p_alpha25[[1]],
                                              huprot_posteriors_p_alpha25[[2]],
                                              huprot_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha25Model3.csv")
  write.csv(huprot_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha25,r=p,df_rY = huprot_mcmc_data_p_alpha25[[1]],type=huprot_mcmc_data_p_alpha25[[4]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha25Model3.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}


##extract posteriors with analysis model 2 data generating D
for(p in 1:5){
  huprot_mcmc_data_p_alpha1<- get_huprot_data_list_alpha1(r=p,Y=model_1_Y_HP[,p],
                                                          type=type_H,I=I_h,K=K_h)
  huprot_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                              data_list = huprot_mcmc_data_p_alpha1,
                                              updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha1<- as.data.frame(rbind(huprot_posteriors_p_alpha1[[1]],
                                             huprot_posteriors_p_alpha1[[2]],
                                             huprot_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha1.csv")
  write.csv(huprot_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha1,r=p,df_rY = huprot_mcmc_data_p_alpha1[[1]],type=huprot_mcmc_data_p_alpha1[[4]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha1.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}
##extract posteriors with analysis model 2 data generating E

for(p in 1:5){
  huprot_mcmc_data_p_alpha1<- get_huprot_data_list_alpha1(r=p,Y=model_2_Y_HP[,p],
                                                          type=type_H,I=I_h,K=K_h)
  huprot_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                              data_list = huprot_mcmc_data_p_alpha1,
                                              updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha1<- as.data.frame(rbind(huprot_posteriors_p_alpha1[[1]],
                                             huprot_posteriors_p_alpha1[[2]],
                                             huprot_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha1Model2.csv")
  write.csv(huprot_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha1,r=p,df_rY = huprot_mcmc_data_p_alpha1[[1]],type=huprot_mcmc_data_p_alpha1[[4]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha1Model2.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}

##extract posteriors with analysis model 2 data generating F
for(p in 1:5){
  huprot_mcmc_data_p_alpha1<- get_huprot_data_list_alpha1(r=p,Y=model_3_Y_HP[,p],
                                                          type=type_H,I=I_h,K=K_h)
  huprot_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                              data_list = huprot_mcmc_data_p_alpha1,
                                              updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha1<- as.data.frame(rbind(huprot_posteriors_p_alpha1[[1]],
                                             huprot_posteriors_p_alpha1[[2]],
                                             huprot_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha1Model3.csv")
  write.csv(huprot_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha1,r=p,df_rY = huprot_mcmc_data_p_alpha1[[1]],type=huprot_mcmc_data_p_alpha1[[4]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha1Model3.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}

#####################################
###Get SEL ranks from posteriors####
###################################
library(data.table)


read_and_process_arrays_25_1<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_alpha25.csv"), select=c(19:15403))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_25_2<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_alpha25Model2.csv"), select=c(19:15403))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_25_3<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_alpha25Model3.csv"), select=c(19:15403))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_1_1<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_alpha1.csv"), select=c(19:15403))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_1_2<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_alpha1Model2.csv"), select=c(19:15403))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_1_3<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_alpha1Model3.csv"), select=c(19:15403))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}


get_huprot_array_ranks<- function(posterior_matrix){
  ranks_posterior<- apply(posterior_matrix,1,rank)
  array_ranks_bar<- apply(t(ranks_posterior),2,mean)
  
  array_ranks_hat<- rank(array_ranks_bar)
  
  return(list(array_ranks_bar, array_ranks_hat, t(ranks_posterior)))
}


huprot_array_ranks_alpha_1_1_bar_df<- matrix(NA,nrow=5,ncol=15385)
huprot_array_ranks_alpha_1_2_bar_df<- matrix(NA,nrow=5,ncol=15385)
huprot_array_ranks_alpha_1_3_bar_df<- matrix(NA,nrow=5,ncol=15385)

huprot_array_ranks_alpha_25_1_bar_df<- matrix(NA,nrow=5,ncol=15385)
huprot_array_ranks_alpha_25_2_bar_df<- matrix(NA,nrow=5,ncol=15385)
huprot_array_ranks_alpha_25_3_bar_df<- matrix(NA,nrow=5,ncol=15385)

for(i in 1:5){
  global_huprot_MCMC_Matrix<- read_and_process_arrays_1_1(i)
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix)
  huprot_array_ranks_alpha_1_1_bar_df[i,]<- huprot_array_ranks[[1]]
}
for(i in 1:5){
  global_huprot_MCMC_Matrix<- read_and_process_arrays_1_2(i)
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix)
  huprot_array_ranks_alpha_1_2_bar_df[i,]<- huprot_array_ranks[[1]]
}

for(i in 1:5){
  global_huprot_MCMC_Matrix<- read_and_process_arrays_1_3(i)
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix)
  huprot_array_ranks_alpha_1_3_bar_df[i,]<- huprot_array_ranks[[1]]
}

for(i in 1:5){
  global_huprot_MCMC_Matrix<- read_and_process_arrays_25_1(i)
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix)
  huprot_array_ranks_alpha_25_1_bar_df[i,]<- huprot_array_ranks[[1]]
}

for(i in 1:5){
  global_huprot_MCMC_Matrix<- read_and_process_arrays_25_2(i)
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix)
  huprot_array_ranks_alpha_25_2_bar_df[i,]<- huprot_array_ranks[[1]]
}

for(i in 1:5){
  global_huprot_MCMC_Matrix<- read_and_process_arrays_25_3(i)
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix)
  huprot_array_ranks_alpha_25_3_bar_df[i,]<- huprot_array_ranks[[1]]
}




write.csv(huprot_array_ranks_alpha_1_1_bar_df,"huprot_array_ranks_alpha_1_1_bar_df.csv")
write.csv(huprot_array_ranks_alpha_1_2_bar_df,"huprot_array_ranks_alpha_1_2_bar_df.csv")
write.csv(huprot_array_ranks_alpha_1_3_bar_df,"huprot_array_ranks_alpha_1_3_bar_df.csv")

write.csv(huprot_array_ranks_alpha_25_1_bar_df,"huprot_array_ranks_alpha_25_1_bar_df.csv")
write.csv(huprot_array_ranks_alpha_25_2_bar_df,"huprot_array_ranks_alpha_25_2_bar_df.csv")
write.csv(huprot_array_ranks_alpha_25_3_bar_df,"huprot_array_ranks_alpha_25_3_bar_df.csv")