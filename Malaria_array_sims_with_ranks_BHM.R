library(R2jags)
library(rjags)
library(dplyr)



beta_cons<- beta(2.5,2.5)
c_alpha<- sd(qnorm(qbeta(1:10000/10001,2.5,2.5)))


###########################################
########Malaria Generating################
##########################################

#generate (10 arrays): 

model_1_Y_Malaria<- matrix(NA, nrow=1156,ncol=10)

S_1<- matrix(NA,1038,10)
mu_ctrl<- matrix(NA,17,10)
tau_active<- rep(NA,10)
mu_active<- rep(NA,10)
sigma<- rep(NA,10)
error<- matrix(NA, 1156,10)
signal<- matrix(NA,1156,10)

for(i in 1:10){
  sigma[i]<- runif(1,0,10)
  error[,i]<- sigma[i]*(1/c_alpha)*qnorm(rbeta(1156,2.5,2.5))
  mu_active[i]<- runif(1,-10^2,10^2)
  tau_active[i]<- runif(1,0,10^3)
  mu_ctrl[,i]<- runif(17,-10^2,10^2)
  S_1[,i]<- rnorm(1038,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],5),rep(mu_ctrl[2,i],4),rep(mu_ctrl[3,i],5),rep(mu_ctrl[4,i],24),
                 rep(mu_ctrl[5,i],32),rep(mu_ctrl[6,i],4),rep(mu_ctrl[7,i],4),rep(mu_ctrl[8,i],4),
                 rep(mu_ctrl[9,i],4),rep(mu_ctrl[10,i],4),rep(mu_ctrl[11,i],4),rep(mu_ctrl[12,i],4),
                 rep(mu_ctrl[13,i],4),rep(mu_ctrl[14,i],4),rep(mu_ctrl[15,i],4),rep(mu_ctrl[16,i],4),
                 rep(mu_ctrl[17,i],4),S_1[,i])
  model_1_Y_Malaria[,i]<- signal[,i]+error[,i]
}


model_2_Y_Malaria<- matrix(NA, nrow=1156,ncol=10)

S_2<- matrix(NA,1038,10)
mu_ctrl<- matrix(NA,17,10)
tau_active<- rep(NA,10)
mu_active<- rep(NA,10)
sigma<- rep(NA,10)
error<- matrix(NA, 1156,10)
signal<- matrix(NA,1156,10)

for(i in 1:10){
  sigma[i]<- runif(1,0,10)
  error[,i]<- rnorm(1156,0,sigma[i])
  mu_active[i]<- runif(1,-10^2,10^2)
  tau_active[i]<- runif(1,0,10^3)
  mu_ctrl[,i]<- runif(17,-10^2,10^2)
  S_2[,i]<- rnorm(1038,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],5),rep(mu_ctrl[2,i],4),rep(mu_ctrl[3,i],5),rep(mu_ctrl[4,i],24),
                 rep(mu_ctrl[5,i],32),rep(mu_ctrl[6,i],4),rep(mu_ctrl[7,i],4),rep(mu_ctrl[8,i],4),
                 rep(mu_ctrl[9,i],4),rep(mu_ctrl[10,i],4),rep(mu_ctrl[11,i],4),rep(mu_ctrl[12,i],4),
                 rep(mu_ctrl[13,i],4),rep(mu_ctrl[14,i],4),rep(mu_ctrl[15,i],4),rep(mu_ctrl[16,i],4),
                 rep(mu_ctrl[17,i],4),S_2[,i])
  model_2_Y_Malaria[,i]<- signal[,i]+error[,i]
}

model_3_Y_Malaria<- matrix(NA, nrow=1156,ncol=10)

S_3<- matrix(NA,1038,10)
mu_ctrl<- matrix(NA,17,10)
tau_active<- rep(NA,10)
mu_active<- rep(NA,10)
sigma<- rep(NA,10)
error<- matrix(NA, 1156,10)
signal<- matrix(NA,1156,10)

for(i in 1:10){
  sigma[i]<- runif(1,0,10^3)
  error[,i]<- rnorm(1156,0,sigma[i])
  mu_active[i]<- runif(1,-10^2,10^2)
  tau_active[i]<- runif(1,0,10^3)
  mu_ctrl[,i]<- runif(17,-10^2,10^2)
  S_3[,i]<- rnorm(1038,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],5),rep(mu_ctrl[2,i],4),rep(mu_ctrl[3,i],5),rep(mu_ctrl[4,i],24),
                 rep(mu_ctrl[5,i],32),rep(mu_ctrl[6,i],4),rep(mu_ctrl[7,i],4),rep(mu_ctrl[8,i],4),
                 rep(mu_ctrl[9,i],4),rep(mu_ctrl[10,i],4),rep(mu_ctrl[11,i],4),rep(mu_ctrl[12,i],4),
                 rep(mu_ctrl[13,i],4),rep(mu_ctrl[14,i],4),rep(mu_ctrl[15,i],4),rep(mu_ctrl[16,i],4),
                 rep(mu_ctrl[17,i],4),S_3[,i])
  model_3_Y_Malaria[,i]<- signal[,i]+error[,i]
}




I_m<- 1156
K_m<- 1055
ones<- rep(1,I_m)

type_1_M<- c(rep(1,5),rep(2,4),rep(3,5),rep(4,24),rep(5,32),rep(6,4),rep(7,4),rep(8,4),rep(9,4),
             rep(10,4),rep(11,4),rep(12,4),rep(13,4),rep(14,4),rep(15,4),rep(16,4),rep(17,4),
             rep(18:1055))

get_malaria_data_list_alpha25<- function(r,Y,type,beta,c,I,K,ones){
  df_r_M<- as.data.frame(cbind(Y,type))
  colnames(df_r_M)<- c("Y", "Type_1")
  df_r_M<- df_r_M[order(df_r_M$Type_1),c("Y", "Type_1")]
  df_r_M$Y<- as.numeric(df_r_M$Y)
  df_r_M$Type_1<- as.numeric(df_r_M$Type_1)
  
  return(list('Y'=df_r_M$Y,"I"=I, "K"=K,"type_1"=df_r_M$Type_1,
              'ones'=ones,'beta_cons'=beta,'c_alpha'=c))
}




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
tau_s~dunif(0,10^6)
sigma~dunif(0,10^3)
}
"



get_percentile_residuals_malaria<- function(posterior, r, df_rY){
  percentile<- ecdf(c(c(posterior[[1]][,18:1055]),c(posterior[[2]][,18:1055]),c(posterior[[3]][,18:1055])))
  percentile_residuals_r<- rep(NA,1038)
  for(l in 1:1038){
    percentile_residuals_r[l]<- percentile(df_rY[l+118])
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


#extract postierors when analysis model is Model 1- data generating is A


for(p in 1:10){
  malaria_mcmc_data_p_alpha25<- get_malaria_data_list_alpha25(r=p,Y=model_1_Y_Malaria[,p],type=type_1_M,
                                                              beta=beta_cons,c=c_alpha,I=1156,K=1055,ones = ones)
  malaria_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                                data_list = malaria_mcmc_data_p_alpha25,
                                                updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=1)
  
  malaria_post_p_alpha25<- as.data.frame(rbind(malaria_posteriors_p_alpha25[[1]],
                                               malaria_posteriors_p_alpha25[[2]],
                                               malaria_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha25.csv")
  write.csv(malaria_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha25,r=p,df_rY = malaria_mcmc_data_p_alpha25[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha25.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}

#extract posteiors when analysis model is model 1 data generating is B

for(p in 1:10){
  malaria_mcmc_data_p_alpha25<- get_malaria_data_list_alpha25(r=p,Y=model_2_Y_Malaria[,p],type=type_1_M,
                                                              beta=beta_cons,c=c_alpha,I=1156,K=1055,ones = ones)
  malaria_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                                data_list = malaria_mcmc_data_p_alpha25,
                                                updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  malaria_post_p_alpha25<- as.data.frame(rbind(malaria_posteriors_p_alpha25[[1]],
                                               malaria_posteriors_p_alpha25[[2]],
                                               malaria_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha25Model2.csv")
  write.csv(malaria_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha25,r=p,df_rY = malaria_mcmc_data_p_alpha25[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha25Model2.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}

#extract posteiors when analysis model is model 1 data generating is C

for(p in 1:10){
  malaria_mcmc_data_p_alpha25<- get_malaria_data_list_alpha25(r=p,Y=model_3_Y_Malaria[,p],type=type_1_M,
                                                              beta=beta_cons,c=c_alpha,I=1156,K=1055,ones = ones)
  malaria_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                                data_list = malaria_mcmc_data_p_alpha25,
                                                updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  malaria_post_p_alpha25<- as.data.frame(rbind(malaria_posteriors_p_alpha25[[1]],
                                               malaria_posteriors_p_alpha25[[2]],
                                               malaria_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha25Model3.csv")
  write.csv(malaria_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha25,r=p,df_rY = malaria_mcmc_data_p_alpha25[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha25Model3.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}



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
tau_s~dunif(0,10^6)
sigma~dunif(0,10^3)
}
"

get_malaria_data_list_alpha1<- function(r,Y,type,I,K){
  df_r_M<- as.data.frame(cbind(Y,type))
  colnames(df_r_M)<- c("Y", "Type_1")
  df_r_M<- df_r_M[order(df_r_M$Type_1),]
  df_r_M$Y<- as.numeric(df_r_M$Y)
  df_r_M$Type_1<- as.numeric(df_r_M$Type_1)
  
  return(list('Y'=df_r_M$Y,"I"=I, "K"=K,"type_1"=df_r_M$Type_1))
}
########################################
##########Extract Posteriors############
#######################################

#extract posteriors when analysis model is model 2 data generating is model A
for(p in 1:10){
  malaria_mcmc_data_p_alpha1<- get_malaria_data_list_alpha1(r=p,Y=model_1_Y_Malaria[,p],
                                                            type=type_1_M,I=1156,K=1055)
  malaria_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                               data_list = malaria_mcmc_data_p_alpha1,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  malaria_post_p_alpha1<- as.data.frame(rbind(malaria_posteriors_p_alpha1[[1]],
                                              malaria_posteriors_p_alpha1[[2]],
                                              malaria_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha1.csv")
  write.csv(malaria_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha1,r=p,df_rY = malaria_mcmc_data_p_alpha1[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha1.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}


#extract posteriors when analysis model is model 2 data generating is model B
for(p in 1:10){
  malaria_mcmc_data_p_alpha1<- get_malaria_data_list_alpha1(r=p,Y=model_2_Y_Malaria[,p],
                                                            type=type_1_M,I=1156,K=1055)
  malaria_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                               data_list = malaria_mcmc_data_p_alpha1,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  malaria_post_p_alpha1<- as.data.frame(rbind(malaria_posteriors_p_alpha1[[1]],
                                              malaria_posteriors_p_alpha1[[2]],
                                              malaria_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha1Model2.csv")
  write.csv(malaria_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha1,r=p,df_rY = malaria_mcmc_data_p_alpha1[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha1Model2.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}


#extract posteriors when analysis model is model 2 data generating is model C

for(p in 1:10){
  malaria_mcmc_data_p_alpha1<- get_malaria_data_list_alpha1(r=p,Y=model_3_Y_Malaria[,p],
                                                            type=type_1_M,I=1156,K=1055)
  malaria_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                               data_list = malaria_mcmc_data_p_alpha1,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  malaria_post_p_alpha1<- as.data.frame(rbind(malaria_posteriors_p_alpha1[[1]],
                                              malaria_posteriors_p_alpha1[[2]],
                                              malaria_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha1Model3.csv")
  write.csv(malaria_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha1,r=p,df_rY = malaria_mcmc_data_p_alpha1[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha1Model3.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}
#############################################
######Get SEL ranks using full posteriors####
#############################################
library(data.table)

read_and_process_arrays_25_1<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_alpha25.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_25_2<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_alpha25Model2.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_25_3<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_alpha25Model3.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_1_1<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_alpha1.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_1_2<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_alpha1Model2.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

read_and_process_arrays_1_3<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_alpha1Model3.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}


global_Malaria_MCMC_Matrix_25_1<- list(read_and_process_arrays_25_1(1),read_and_process_arrays_25_1(2),
                                       read_and_process_arrays_25_1(3),read_and_process_arrays_25_1(4),
                                       read_and_process_arrays_25_1(5),read_and_process_arrays_25_1(6),
                                       read_and_process_arrays_25_1(7),read_and_process_arrays_25_1(8),
                                       read_and_process_arrays_25_1(9),read_and_process_arrays_25_1(10))

global_Malaria_MCMC_Matrix_25_2<- list(read_and_process_arrays_25_2(1),read_and_process_arrays_25_2(2),
                                       read_and_process_arrays_25_2(3),read_and_process_arrays_25_2(4),
                                       read_and_process_arrays_25_2(5),read_and_process_arrays_25_2(6),
                                       read_and_process_arrays_25_2(7),read_and_process_arrays_25_2(8),
                                       read_and_process_arrays_25_2(9),read_and_process_arrays_25_2(10))

global_Malaria_MCMC_Matrix_25_3<- list(read_and_process_arrays_25_3(1),read_and_process_arrays_25_3(2),
                                       read_and_process_arrays_25_3(3),read_and_process_arrays_25_3(4),
                                       read_and_process_arrays_25_3(5),read_and_process_arrays_25_3(6),
                                       read_and_process_arrays_25_3(7),read_and_process_arrays_25_3(8),
                                       read_and_process_arrays_25_3(9),read_and_process_arrays_25_3(10))


############

get_malaria_array_ranks<- function(posterior_matrix){
  ranks_posterior<- apply(posterior_matrix,1,rank)
  array_ranks_bar<- apply(t(ranks_posterior),2,mean)
  
  array_ranks_hat<- rank(array_ranks_bar)
  
  return(list(array_ranks_bar, array_ranks_hat, t(ranks_posterior)))
}

malaria_array_ranks_alpha_25_1_bar_df<- matrix(NA,nrow=10,ncol=1038)
malaria_array_ranks_alpha_25_2_bar_df<- matrix(NA,nrow=10,ncol=1038)
malaria_array_ranks_alpha_25_3_bar_df<- matrix(NA,nrow=10,ncol=1038)

for(i in 1:10){
  malaria_array_ranks_alpha_25_1<- get_malaria_array_ranks(global_Malaria_MCMC_Matrix_25_1[[i]])
  malaria_array_ranks_alpha_25_2<- get_malaria_array_ranks(global_Malaria_MCMC_Matrix_25_2[[i]])
  malaria_array_ranks_alpha_25_3<- get_malaria_array_ranks(global_Malaria_MCMC_Matrix_25_3[[i]])
  malaria_array_ranks_alpha_25_1_bar_df[i,]<- malaria_array_ranks_alpha_25_1[[1]]
  malaria_array_ranks_alpha_25_2_bar_df[i,]<- malaria_array_ranks_alpha_25_2[[1]]
  malaria_array_ranks_alpha_25_3_bar_df[i,]<- malaria_array_ranks_alpha_25_3[[1]]
}

write.csv(malaria_array_ranks_alpha_25_1_bar_df,"Malaria_array_ranks_alpha_25_1_bar_df.csv")
write.csv(malaria_array_ranks_alpha_25_2_bar_df,"Malaria_array_ranks_alpha_25_2_bar_df.csv")
write.csv(malaria_array_ranks_alpha_25_3_bar_df,"Malaria_array_ranks_alpha_25_3_bar_df.csv")


write.csv(S_1,"S_1.csv")
write.csv(S_2,"S_2.csv")
write.csv(S_3,"S_3.csv")
################################

Model_1_ranks<- rank(apply(apply(S_1,2,rank),1,mean))
Model_2_ranks<- rank(apply(apply(S_2,2,rank),1,mean))
Model_3_ranks<- rank(apply(apply(S_3,2,rank),1,mean))



alpha_25_1_ranks<- rank(apply(apply(malaria_array_ranks_alpha_25_1_bar_df,1,rank),1,mean))
alpha_25_2_ranks<- rank(apply(apply(malaria_array_ranks_alpha_25_2_bar_df,1,rank),1,mean))
alpha_25_3_ranks<- rank(apply(apply(malaria_array_ranks_alpha_25_3_bar_df,1,rank),1,mean))


alpha_1_1_ranks<- rank(apply(apply(malaria_array_ranks_alpha_1_1_bar_df,1,rank),1,mean))
alpha_1_2_ranks<- rank(apply(apply(malaria_array_ranks_alpha_1_2_bar_df,1,rank),1,mean))
alpha_1_3_ranks<- rank(apply(apply(malaria_array_ranks_alpha_1_3_bar_df,1,rank),1,mean))

