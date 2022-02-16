library(data.table)
library(abind)
library(dplyr)
########################################
#Reading in  malaria data
########################################

read_and_process_arrays<- function(index){
  array<- fread(paste0("Malaria_Posteriors_Array_",index,"_hatSigma.csv"), select=c(19:1056))
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}




global_Malaria_MCMC_Matrix<- list(read_and_process_arrays(1),read_and_process_arrays(2),
                                  read_and_process_arrays(3),read_and_process_arrays(4),
                                  read_and_process_arrays(5),read_and_process_arrays(6),
                                  read_and_process_arrays(7),read_and_process_arrays(8),
                                  read_and_process_arrays(9),read_and_process_arrays(10))


quantile_new<- quantile(as.vector(unlist(global_Malaria_MCMC_Matrix)),c(0.9,0.95,0.97))



print("done reading malaria data")
#################################
#getting array bar and hat ranks
#################################
get_malaria_array_ranks<- function(posterior_matrix){
  ranks_posterior<- apply(posterior_matrix,1,rank)
  array_ranks_bar<- apply(t(ranks_posterior),2,mean)
  
  array_ranks_hat<- rank(array_ranks_bar)
  
  return(list(array_ranks_bar, array_ranks_hat, t(ranks_posterior)))
}

malaria_array_ranks_bar_df<- matrix(NA,nrow=503,ncol=1038)
malaria_array_ranks_hat_df<- matrix(NA,nrow=503,ncol=1038)
for(i in 1:503){
  malaria_array_ranks<- get_malaria_array_ranks(global_Malaria_MCMC_Matrix[,,i])
  malaria_array_ranks_bar_df[i,]<- malaria_array_ranks[[1]]
  malaria_array_ranks_hat_df[i,]<-malaria_array_ranks[[2]]
}

write.csv(malaria_array_ranks_bar_df,"Malaria_array_ranks_bar_df.csv")
write.csv(malaria_array_ranks_hat_df,"Malaria_array_ranks_hat_df.csv")

#######################################
#Global ranks, max and min ranks
######################################


global_ranks_bar<-apply(malaria_array_ranks_bar_df,2,mean)
global_ranks_hat<-apply(malaria_array_ranks_hat_df,2,mean)



write.csv(global_ranks_bar,"Malaria_global_ranks_bar.csv")
write.csv(global_ranks_hat,"Malaria_global_ranks_hat.csv")
print("done doing global ranks")
#############
#gamma ranks
#############

malaria_percentiles<- quantile(as.vector(unlist(global_Malaria_MCMC_Matrix)),c(0.97,0.95,0.90))


malaria_gamma_mean_array_top97<- matrix(NA,1038,10)

for(i in 1:10){
  for(j in 1:1038){
    
    malaria_gamma_mean_array_top97[j,i]<-length(which(as.matrix(global_Malaria_MCMC_Matrix[[i]])[,j]>malaria_percentiles[1]))/10000
  }
  
}



malaria_gamma_mean_array_top_95<- matrix(NA,1038,10)

for(i in 1:10){
  for(j in 1:1038){
    
    malaria_gamma_mean_array_top95[j,i]<-length(which(as.matrix(global_Malaria_MCMC_Matrix[[i]])[,j]>malaria_percentiles[2]))/10000
  }
  
}


malaria_gamma_mean_array_top_90<- matrix(NA,1038,10)

for(i in 1:10){
  for(j in 1:1038){
    
    malaria_gamma_mean_array_top90[j,i]<-length(which(as.matrix(global_Malaria_MCMC_Matrix[[i]])[,j]>malaria_percentiles[3]))/10000
  }
  
}



write.csv(malaria_gamma_mean_array_top97,"Malaria_gamma_array_top97.csv")
write.csv(malaria_gamma_mean_array_top95,"Malaria_gamma_array_top95.csv")
write.csv(malaria_gamma_mean_array_top90,"Malaria_gamma_array_top90.csv")

malaria_gamma_mean_top97<- apply(malaria_gamma_mean_array_top97,1,mean)
malaria_gamma_mean_top95<- apply(malaria_gamma_mean_array_top95,1,mean)
malaria_gamma_mean_top90<- apply(malaria_gamma_mean_array_top90,1,mean)

write.csv(malaria_gamma_mean_top97,"Malaria_gamma_top97.csv")
write.csv(malaria_gamma_mean_top95,"Malaria_gamma_top95.csv")
write.csv(malaria_gamma_mean_top90,"Malaria_gamma_top90.csv")


########################################
############Read in Huprot data#########
########################################
read_and_process_arrays<- function(index){
  array<- fread(paste0("HuProt_Posteriors_Array_",index,"_hatSigma.csv"), select=c(19:15403),)
  array_clean<- array[sample(1:30000,size=10000),]
  return(array_clean)
}

global_huprot_MCMC_Matrix<- array(as.numeric(NA), c(10000,15385,100))

for(i in 1:100){
  global_huprot_MCMC_Matrix[,,i]<- read_and_process_arrays(i)
}


print("done reading huprot data")
#################################
#getting array bar and hat ranks
#################################
get_huprot_array_ranks<- function(posterior_matrix){
  ranks_posterior<- apply(posterior_matrix,1,rank)
  array_ranks_bar<- apply(t(ranks_posterior),2,mean)
  
  array_ranks_hat<- rank(array_ranks_bar)
  
  return(list(array_ranks_bar, array_ranks_hat, t(ranks_posterior)))
}

huprot_array_ranks_bar_df<- matrix(NA,nrow=100,ncol=15385)
huprot_array_ranks_hat_df<- matrix(NA,nrow=100,ncol=15385)
for(i in 1:100){
  huprot_array_ranks<- get_huprot_array_ranks(global_huprot_MCMC_Matrix[,,i])
  huprot_array_ranks_bar_df[i,]<- huprot_array_ranks[[1]]
  huprot_array_ranks_hat_df[i,]<-huprot_array_ranks[[2]]
}

write.csv(huprot_array_ranks_bar_df,"HuProt_array_ranks_bar_df.csv")
write.csv(huprot_array_ranks_hat_df,"HuProt_array_ranks_hat_df.csv")
print("done doing huprot array ranks")
#######################################
#Global ranks, max and min ranks
######################################


global_ranks_bar<-apply(huprot_array_ranks_bar_df,2,mean)
global_ranks_hat<-apply(huprot_array_ranks_hat_df,2,mean)


write.csv(global_ranks_bar,"HuProt_global_ranks_bar.csv")
write.csv(global_ranks_hat,"HuProt_global_ranks_hat.csv")

print("done doing huprot global ranks")
#############
#gamma ranks
#############

huprot_percentiles<- quantile(as.vector(global_huprot_MCMC_Matrix),c(0.97,0.95,0.90))


huprot_gamma_mean_array_top97<- matrix(NA,15385,100)

for(i in 1:100){
  for(j in 1:15385){
    
    huprot_gamma_mean_array_top97[j,i]<- length(which(global_huprot_MCMC_Matrix[,j,i]>huprot_percentiles[1]))/10000
  }
  
}


huprot_gamma_mean_array_top95<- matrix(NA,15385,100)

for(i in 1:100){
  for(j in 1:15385){
    
    huprot_gamma_mean_array_top95[j,i]<- length(which(global_huprot_MCMC_Matrix[,j,i]>huprot_percentiles[2]))/10000
  }
  
}

huprot_gamma_mean_array_top90<- matrix(NA,15385,100)

for(i in 1:100){
  for(j in 1:15385){
    
    huprot_gamma_mean_array_top90[j,i]<- length(which(global_huprot_MCMC_Matrix[,j,i]>huprot_percentiles[3]))/10000
  }
  
}





write.csv(huprot_gamma_mean_array_top97,"HuProt_gamma_array_top97.csv")
write.csv(huprot_gamma_mean_array_top95,"HuProt_gamma_array_top95.csv")
write.csv(huprot_gamma_mean_array_top90,"HuProt_gamma_array_top90.csv")

huprot_gamma_mean_top97<- apply(huprot_gamma_mean_array_top97,1,mean)
huprot_gamma_mean_top95<- apply(huprot_gamma_mean_array_top95,1,mean)
huprot_gamma_mean_top90<- apply(huprot_gamma_mean_array_top90,1,mean)

write.csv(huprot_gamma_mean_top97,"HuProt_gamma_top97.csv")
write.csv(huprot_gamma_mean_top95,"HuProt_gamma_top95.csv")
write.csv(huprot_gamma_mean_top90,"HuProt_gamma_top90.csv")

