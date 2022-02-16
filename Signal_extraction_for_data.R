
library(R2jags)
library(dplyr)



####################
#Constants and Data
###################

beta_cons<- beta(2.5,2.5)
c_alpha<- sd(qnorm(qbeta(1:10000/10001,2.5,2.5)))

I_m<- 1156
K_m<- 1055
ones<- rep(1,I_m)

I_h<-48384
K_h<-15402
ones_h<-rep(1,I_h)

std_ratio_matrix_M<- read.csv("malaria_std_ratio_matrix.csv",header = F)
global_sd_M<- read.csv("Malaria_sd_vector.csv")

std_ratio_matrix_H<- read.csv("HuProt_std_ratio_matrix.csv")
global_sd_H<- read.csv("HuProt_sd_vector.csv")



index_neg_ctrl_empty_M<- c(574,863,866,1152,1155)
index_neg_ctrl_blank_M<- c(289,578,867,1156)
index_neg_ctrl_blank1_M<- c(285,349,695,1144,1147)
index_neg_ctrl_noDNA_M<- c(17,18,19,20,21,22,306,307,308,309,310,311,
                           595,596,597,598,599,600,884,885,886,887,888,889)
index_neg_ctrl_ttbs_M<- c(13,14,15,16,29,30,31,32,302,303,304,305,
                          318,319,320,321,591,592,593,594,607,608,609,610,
                          880,881,882,883,896,897,898,899)
  
nc_index_M<- c(index_neg_ctrl_empty_M,index_neg_ctrl_blank_M,
               index_neg_ctrl_blank1_M,index_neg_ctrl_noDNA_M,index_neg_ctrl_ttbs_M)


HuProt_name_indices<-read.csv("HuProt_name_indices.csv")

index_neg_ctrl_gst10<-which(HuProt_name_indices$Name=="GST10n")
index_neg_ctrl_gst50<-which(HuProt_name_indices$Name=="GST50n")
index_neg_ctrl_gst100<-which(HuProt_name_indices$Name=="GST100n")
index_neg_ctrl_gst200<-which(HuProt_name_indices$Name=="GST200n")
index_neg_ctrl_mouse<-which(HuProt_name_indices$Name=="Mouse-anti-biotin")
index_neg_ctrl_rabbit<-which(HuProt_name_indices$Name=="Rabbit-anti-biotin")
index_neg_ctrl_bsa<-which(HuProt_name_indices$Name=="BSA")
index_neg_ctrl_empty_H<-which(HuProt_name_indices$Name=="empty")
index_neg_ctrl_buffer_H<-which(HuProt_name_indices$Name=="buffer"|HuProt_name_indices$Name=="Buffer")


nc_index_H<- c(index_neg_ctrl_gst10,index_neg_ctrl_gst50,index_neg_ctrl_gst100,
                   index_neg_ctrl_gst200,index_neg_ctrl_mouse,index_neg_ctrl_rabbit,
                   index_neg_ctrl_bsa,index_neg_ctrl_empty_H,index_neg_ctrl_buffer_H)


index_Igg003<-c(11,300,589,878)
index_Igg03<- c(9,298,587,876)
index_Igg3<- c(7,296,585,874)
index_MIgg003<-c(5,294,583,872)
index_MIgg03<- c(3,292,581,870)
index_MIgg3<- c(1,290,579,868)
index_Igg001<- c(12,301,590,879)
index_Igg01<- c(10,299,588,877)
index_Igg1<- c(8,297,586,875)
index_MIgg001<- c(6,295,584,873)
index_MIgg01<- c(4,293,582,871)
index_MIgg1<- c(2,291,580,869)

pc_index_M<- c(index_Igg003,index_Igg03,index_Igg3,index_MIgg003,
               index_MIgg03,index_MIgg3,index_Igg001,index_Igg01,
               index_Igg1,index_MIgg001,index_MIgg01,index_MIgg1)


index_h1<-which(HuProt_name_indices$Name=="H1")
index_h2<- which(HuProt_name_indices$Name=="H2(A+B)")
index_h3<- which(HuProt_name_indices$Name=="H3")
index_h4<-which(HuProt_name_indices$Name=="H4")
index_alexaIgg<- which(HuProt_name_indices$Name=="IgG488/594")
index_RhodaAlexa<- which(HuProt_name_indices$Name=="Rhodamine+IgG647")
index_BiotinBSA<- which(HuProt_name_indices$Name=="Biotin-BSA")
index_MouseIGM<- which(HuProt_name_indices$Name=="Mouse-IgM\xe1")

pc_index_H<- c(index_h1,index_h2,index_h3,index_h4,index_alexaIgg,index_RhodaAlexa,
                   index_BiotinBSA,index_MouseIGM)


type_1_M<- rep(NA,1156)

type_1_M[index_neg_ctrl_blank_M]<-1
type_1_M[index_neg_ctrl_blank1_M]<-2
type_1_M[index_neg_ctrl_empty_M]<-3
type_1_M[index_neg_ctrl_noDNA_M]<-4
type_1_M[index_neg_ctrl_ttbs_M]<-5

type_1_M[index_Igg003]<-6
type_1_M[index_Igg03]<- 7
type_1_M[index_Igg3]<- 8
type_1_M[index_Igg001]<- 9
type_1_M[index_Igg01]<-10
type_1_M[index_Igg1]<-11
type_1_M[index_MIgg003]<-12
type_1_M[index_MIgg03]<-13
type_1_M[index_MIgg3]<-14
type_1_M[index_MIgg001]<-15
type_1_M[index_MIgg01]<-16
type_1_M[index_MIgg1]<-17

for(i in 1:1038){
  type_1_M[-c(pc_index_M,nc_index_M)][i]<- (i+17)
}

type_1_H<- rep(NA,48384)

type_1_H[index_neg_ctrl_gst10]<-1
type_1_H[index_neg_ctrl_gst50]<-2
type_1_H[index_neg_ctrl_gst100]<-3
type_1_H[index_neg_ctrl_gst200]<-4
type_1_H[index_neg_ctrl_mouse]<-5
type_1_H[index_neg_ctrl_rabbit]<-6
type_1_H[index_neg_ctrl_bsa]<-7
type_1_H[index_neg_ctrl_buffer_H]<-8
type_1_H[index_neg_ctrl_empty_H]<-9

type_1_H[index_h1]<-10
type_1_H[index_h2]<- 11
type_1_H[index_h3]<- 12
type_1_H[index_h4]<- 13
type_1_H[index_alexaIgg]<-14
type_1_H[index_RhodaAlexa]<-15
type_1_H[index_BiotinBSA]<-16
type_1_H[index_MouseIGM]<-17

unique_names_noctrl<- unique(HuProt_name_indices$ID[-c(pc_index_H,nc_index_H)])
unique_names_ID_num<- rep((1+17):(length(unique_names_noctrl)+17))

active_spot_IDs<- as.data.frame(cbind(unique_names_noctrl,unique_names_ID_num))
colnames(active_spot_IDs)<- c("Name", "Number")

active_spot_IDs$Number<- as.numeric(active_spot_IDs$Number)

for(i in 1:length(unique_names_noctrl)){
  type_1_H[which(HuProt_name_indices$Name==active_spot_IDs$Name[i])]<- active_spot_IDs$Number[i]
}


###############################################
#Functions to get data list for Malaria Arrays
##############################################


get_malaria_data_list_alpha25<- function(r,Y,type,beta,c,I,K,ones){
  df_r_M<- as.data.frame(cbind(Y,type))
  colnames(df_r_M)<- c("Y", "Type_1")
  df_r_M<- df_r_M[order(df_r_M$Type_1),c("Y", "Type_1")]
  df_r_M$Y<- as.numeric(df_r_M$Y)
  df_r_M$Type_1<- as.numeric(df_r_M$Type_1)
  
  return(list('Y'=df_r_M$Y,"I"=I, "K"=K,"type_1"=df_r_M$Type_1,
              'ones'=ones,'beta_cons'=beta,'c_alpha'=c))
}


get_malaria_data_list_alpha1<- function(r,Y,type,I,K){
  df_r_M<- as.data.frame(cbind(Y,type))
  colnames(df_r_M)<- c("Y", "Type_1")
  df_r_M<- df_r_M[order(df_r_M$Type_1),]
  df_r_M$Y<- as.numeric(df_r_M$Y)
  df_r_M$Type_1<- as.numeric(df_r_M$Type_1)
  
  return(list('Y'=df_r_M$Y,"I"=I, "K"=K,"type_1"=df_r_M$Type_1))
}

get_malaria_data_list_hatSigma<- function(r,Y,type,I,K,sigma){
  df_r_M<- as.data.frame(cbind(Y,type))
  colnames(df_r_M)<- c("Y", "Type_1")
  df_r_M<- df_r_M[order(df_r_M$Type_1),]
  df_r_M$Y<- as.numeric(df_r_M$Y)
  df_r_M$Type_1<- as.numeric(df_r_M$Type_1)
  
  return(list('Y'=df_r_M$Y,"I"=I, "K"=K,"type_1"=df_r_M$Type_1, 'sigma'=sigma))
}


###############################################
#Functions to get data list for HuProt Arrays
##############################################


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


#############################
#Function to get posteriors
############################


get_posteriors<- function(model, data_list, updating_params,r){
  model_1 <-textConnection(model)
  model_initial<- jags.model(model_1,data=data_list, n.chains=3,
                             n.adapt=1000,quiet=FALSE)
  
  update(model_initial,5000)
  
  model_coda_r <- coda.samples(model_initial,updating_params,n.iter=10000)
  
  return(model_coda_r)
}


################################
#Percentile Residuals Functions
###############################
get_percentile_residuals_malaria<- function(posterior, r, df_rY){
  percentile<- ecdf(c(c(posterior[[1]][,18:1055]),c(posterior[[2]][,18:1055]),c(posterior[[3]][,18:1055])))
  percentile_residuals_r<- rep(NA,1038)
  for(l in 1:1038){
    percentile_residuals_r[l]<- percentile(df_rY[l+118])
  }
  percentile_residuals_r_transformed<- qnorm(percentile_residuals_r)
  return(percentile_residuals_r_transformed)
}



get_percentile_residuals_huprot<- function(posterior, r, df_rY){
  percentile<- ecdf(c(as.vector(as.matrix(posterior[[1]][,18:15402])),as.vector(as.matrix(posterior[[2]][,18:15402])),as.vector(as.matrix(posterior[[3]][,18:15402]))))
  percentile_residuals_r<- rep(NA,15385)
  for(l in 1:15385){
    percentile_residuals_r[l]<- percentile(df_rY[l+4608])
  }
  percentile_residuals_r_transformed<- qnorm(percentile_residuals_r)
  return(percentile_residuals_r_transformed)
}

###################
#Malaria For Loops
##################

#alpha=2.5
for(p in 1:503){
  malaria_mcmc_data_p_alpha25<- get_malaria_data_list_alpha25(r=p,Y=std_ratio_matrix_M[,p+1],type=type_1_M,
                                                      beta=beta_cons,c=c_alpha,I=1156,K=1055,ones = ones)
  malaria_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                        data_list = malaria_mcmc_data_p_alpha25,
                                        updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  malaria_post_p_alpha25<- as.data.frame(rbind(malaria_posteriors_p_alpha25[[1]],
                                               malaria_posteriors_p_alpha25[[2]],
                                              malaria_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_alpha25.csv")
  write.csv(malaria_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_malaria(malaria_posteriors_p_alpha25,r=p,df_rY = malaria_mcmc_data_p_alpha25[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_alpha25.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}

#alpha=1

for(p in 1:503){
  malaria_mcmc_data_p_alpha1<- get_malaria_data_list_alpha1(r=p,Y=std_ratio_matrix_M[,p+1],
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

#hatSigma

for(p in 1:503){
  malaria_mcmc_data_p_hatSigma<- get_malaria_data_list_hatSigma(r=p,Y=std_ratio_matrix_M[,p+1],
                                                            type=type_1_M,I=1156,K=1055,sigma = global_sd_M$x[p])
  malaria_posteriors_p_hatSigma<- get_posteriors(model = model_hatSigma,
                                               data_list = malaria_mcmc_data_p_hatSigma,
                                               updating_params = c("S", "mu_s",'tau_s'), r=p)
  
  malaria_post_p_hatSigma<- as.data.frame(rbind(malaria_posteriors_p_hatSigma[[1]],
                                              malaria_posteriors_p_hatSigma[[2]],
                                              malaria_posteriors_p_hatSigma[[3]]))
  file_name_posteriors<- paste0("Malaria_Posteriors_Array_",p,"_hatSigma.csv")
  write.csv(malaria_post_p_hatSigma,file_name_posteriors)
  
  percentiles_Array_p_hatSigma<-get_percentile_residuals_malaria(malaria_posteriors_p_hatSigma,r=p,df_rY = malaria_mcmc_data_p_hatSigma[[1]])
  
  file_name_percentiles<- paste0("Malaria_Percentiles_Array_",p,"_hatSigma.csv")
  write.csv(percentiles_Array_p_hatSigma,file_name_percentiles)
}



###################
#HuProt For Loops
##################

#alpha=2.5
for(p in c(3,23,43,57,83)){
  huprot_mcmc_data_p_alpha25<- get_huprot_data_list_alpha25(r=p,Y=std_ratio_matrix_H[,p+1],type=type_1_H,
                                                              beta=beta_cons,c=c_alpha,I=48384,K=15402,ones = ones_H)
  huprot_posteriors_p_alpha25<- get_posteriors(model = model_alpha_25,
                                                data_list = huprot_mcmc_data_p_alpha25,
                                                updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha25<- as.data.frame(rbind(huprot_posteriors_p_alpha25[[1]],
                                               huprot_posteriors_p_alpha25[[2]],
                                               huprot_posteriors_p_alpha25[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha25.csv")
  write.csv(huprot_post_p_alpha25,file_name_posteriors)
  
  percentiles_Array_p_alpha25<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha25,r=p,df_rY = huprot_mcmc_data_p_alpha25[[1]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha25.csv")
  write.csv(percentiles_Array_p_alpha25,file_name_percentiles)
}

#alpha=1

posterior<- fread("HuProt_Posteriors_Array_43_alpha1.csv")
posterior1<- as.data.frame(posterior)
huprot_mcmc_data_43_alpha1<- get_huprot_data_list_alpha1(r=43,Y=std_ratio_matrix_H[,44],
                                                         type=type_1_H,I=I_h,K=K_h)

get_percentile_residuals_huprot<- function(posterior, r, df_rY){
  percentile<- ecdf(as.vector(as.matrix(posterior1[,19:15403])))
  percentile_residuals_r<- rep(NA,15385)
  for(l in 1:15385){
    percentile_residuals_r[l]<- percentile(huprot_mcmc_data_43_alpha1[[1]][l+4608])
  }
  percentile_residuals_r_transformed<- qnorm(percentile_residuals_r)
  return(percentile_residuals_r_transformed)
}
write.csv(percentile_residuals_r_transformed,"HuProt_Percentiles_Array_43_alpha1.csv")


for(p in c(3,23,43,57,81)){
  huprot_mcmc_data_p_alpha1<- get_huprot_data_list_alpha1(r=p,Y=std_ratio_matrix_H[,p+1],
                                                            type=type_1_H,I=I_h,K=K_h)
  huprot_posteriors_p_alpha1<- get_posteriors(model = model_alpha_1,
                                               data_list = huprot_mcmc_data_p_alpha1,
                                               updating_params = c("S", "mu_s",'tau_s', 'sigma'), r=p)
  
  huprot_post_p_alpha1<- as.data.frame(rbind(huprot_posteriors_p_alpha1[[1]],
                                              huprot_posteriors_p_alpha1[[2]],
                                              huprot_posteriors_p_alpha1[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_alpha1.csv")
  write.csv(huprot_post_p_alpha1,file_name_posteriors)
  
  percentiles_Array_p_alpha1<-get_percentile_residuals_huprot(huprot_posteriors_p_alpha1,r=p,df_rY = huprot_mcmc_data_p_alpha1[[1]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_alpha1.csv")
  write.csv(percentiles_Array_p_alpha1,file_name_percentiles)
}

#hatSigma

for(p in c(3,23,43,57,83)){
  huprot_mcmc_data_p_hatSigma<- get_huprot_data_list_hatSigma(r=p,Y=std_ratio_matrix_H[,p+1],
                                                                type=type_1_H,I=I_h,K=K_h,sigma = global_sd_H$x[p])
  huprot_posteriors_p_hatSigma<- get_posteriors(model = model_hatSigma,
                                                 data_list = huprot_mcmc_data_p_hatSigma,
                                                 updating_params = c("S", "mu_s",'tau_s'), r=p)
  
  huprot_post_p_hatSigma<- as.data.frame(rbind(huprot_posteriors_p_hatSigma[[1]],
                                                huprot_posteriors_p_hatSigma[[2]],
                                                huprot_posteriors_p_hatSigma[[3]]))
  file_name_posteriors<- paste0("HuProt_Posteriors_Array_",p,"_hatSigma.csv")
  write.csv(huprot_post_p_hatSigma,file_name_posteriors)
  
  percentiles_Array_p_hatSigma<-get_percentile_residuals_huprot(huprot_posteriors_p_hatSigma,r=p,df_rY = huprot_mcmc_data_p_hatSigma[[1]])
  
  file_name_percentiles<- paste0("HuProt_Percentiles_Array_",p,"_hatSigma.csv")
  write.csv(percentiles_Array_p_hatSigma,file_name_percentiles)
}



