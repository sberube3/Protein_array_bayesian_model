library(R2jags)
library(runjags)
library(MCMCpack)
library(rjags)
library(coda)
library(bayesplot)
library(distr)
library(mcmcplots)
library(dplyr)
library(tidyr)


#constants

beta_cons<- beta(2.5,2.5)
c_alpha<- sd(qnorm(qbeta(1:10000/10001,2.5,2.5)))


###########################################
#####Malaria Generating=Model 2############
##########################################

#generate (25 arrays): 

model_2_Y_Malaria<- matrix(NA, nrow=1156,ncol=10)

S<- matrix(NA,1038,10)
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
  S[,i]<- rnorm(1038,mu_active[i],tau_active[i])
  signal[,i]<- c(rep(mu_ctrl[1,i],5),rep(mu_ctrl[2,i],4),rep(mu_ctrl[3,i],5),rep(mu_ctrl[4,i],24),
             rep(mu_ctrl[5,i],32),rep(mu_ctrl[6,i],4),rep(mu_ctrl[7,i],4),rep(mu_ctrl[8,i],4),
             rep(mu_ctrl[9,i],4),rep(mu_ctrl[10,i],4),rep(mu_ctrl[11,i],4),rep(mu_ctrl[12,i],4),
             rep(mu_ctrl[13,i],4),rep(mu_ctrl[14,i],4),rep(mu_ctrl[15,i],4),rep(mu_ctrl[16,i],4),
             rep(mu_ctrl[17,i],4),S[,i])
  model_2_Y_Malaria[,i]<- signal[,i]+error[,i]
}





#analysis model is Model 2


I_m<- 1156
K_m<- 1055
ones<- rep(1,I_m)

type<- c(rep(1,5),rep(2,4),rep(3,5),rep(4,24),rep(5,32),rep(6,4),rep(7,4),rep(8,4),rep(9,4),
         rep(10,4),rep(11,4),rep(12,4),rep(13,4),rep(14,4),rep(15,4),rep(16,4),rep(17,4),
         rep(18:1055))

get_malaria_data_list_alpha1<- function(r,Y,type,I,K){
  df_r_M<- as.data.frame(cbind(Y,type))
  colnames(df_r_M)<- c("Y", "Type_1")
  df_r_M<- df_r_M[order(df_r_M$Type_1),]
  df_r_M$Y<- as.numeric(df_r_M$Y)
  df_r_M$Type_1<- as.numeric(df_r_M$Type_1)
  
  return(list('Y'=df_r_M$Y,"I"=I, "K"=K,"type_1"=df_r_M$Type_1))
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
tau_s~dunif(0,10^3)
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


for(p in 1:10){
  malaria_mcmc_data_p_alpha1<- get_malaria_data_list_alpha1(r=p,Y=model_2_Y_Malaria[,p],
                                                            type=type,I=1156,K=1055)
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

percentiles_1<- read.csv("Malaria_percentiles_array_1_alpha1.csv")
percentiles_2<- read.csv("Malaria_percentiles_array_2_alpha1.csv")
percentiles_3<- read.csv("Malaria_percentiles_array_3_alpha1.csv")
percentiles_4<- read.csv("Malaria_percentiles_array_4_alpha1.csv")
percentiles_5<- read.csv("Malaria_percentiles_array_5_alpha1.csv")
percentiles_6<- read.csv("Malaria_percentiles_array_6_alpha1.csv")
percentiles_7<- read.csv("Malaria_percentiles_array_7_alpha1.csv")
percentiles_8<- read.csv("Malaria_percentiles_array_8_alpha1.csv")
percentiles_9<- read.csv("Malaria_percentiles_array_9_alpha1.csv")
percentiles_10<- read.csv("Malaria_percentiles_array_10_alpha1.csv")


malaria_percentiles_1_1<- c(percentiles_1$x,
                            percentiles_2$x,
                            percentiles_3$x,
                            percentiles_4$x,
                            percentiles_5$x,
                            percentiles_6$x,
                            percentiles_7$x,
                            percentiles_8$x,
                            percentiles_9$x,
                            percentiles_10$x)


Malaria_perc_res_df<- as.data.frame(cbind(c(malaria_percentiles_1_1),
                                          rep(c("Model 1"),
                                              length(malaria_percentiles_1_1))))

colnames(Malaria_perc_res_df)<- c("P_double_dagger","Model")

Malaria_perc_res_df$P_double_dagger<- as.numeric(Malaria_perc_res_df$P_double_dagger)

p_res_malaria<- ggplot(Malaria_perc_res_df,aes(sample=P_double_dagger))+
  stat_qq(size=0.75)+
  geom_abline(slope=1,intercept=0,color="red")+
  facet_grid(~Model)+
  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                   legend.title=element_text(size=5),legend.text=element_text(size=5))+
  ggtitle("Malaria Arrays")+xlab("Theoretical N(0,1)")+ylab("Percentile Residuals")+
  scale_x_continuous(breaks = round(seq(-5,5, by = 1),1))


################################