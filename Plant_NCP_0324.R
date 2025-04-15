#########Library############
library(bipartite)
library(igraph)
library(rnetcarto)
library(devtools)
library(BayesLogit)
library(withr)
library(knitr)
library(Hmsc)
library(fields)
library(MASS)
library(parallel)
library(doParallel)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(car)
library(effects)
library(emmeans)
library(reshape)
library(corrplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
#install.packages("fmsb")
library(fmsb)

######Part 1##############
#upload data set for the four different glaciers and the persistence prob
#clean the data set and perform simulation

##Load general data set
Prob_df<- read.csv("Simone_NCP/tabone.csv")##persisence probability data set. Use p_m (mean prob)
NCP<- read.csv("Simone_NCP/Plants-NCP_n.csv", sep = ";")##Plants, their stage and the NCP database
#NCP<- NCP[c(1:99),]## keep only the plant with an attibuted stage
NCP_cat<-read.csv("Simone_NCP/Plants-NCP_cat.csv", sep = ";")#caetgories of NCP

##Load glacier dataset
amola<- read.csv("Simone_NCP/amola_matb.csv")
amola_env<- read.csv("Simone_NCP/amola_env.csv")

cedec<- read.csv("Simone_NCP/cedec_matb.csv")
cedec_env<-read.csv("Simone_NCP/cedec_env.csv")

rutor<-read.csv("Simone_NCP/rutor_matb.csv")
rutor_env<-read.csv("Simone_NCP/rutor_env.csv")

trobio<-read.csv("Simone_NCP/trobio_matb.csv")
trobio_env<-read.csv("Simone_NCP/trobio_env.csv")

glacierbrutlist<- list(amola, cedec, rutor,trobio)
glacierenvlist<-list(amola_env, cedec_env, rutor_env, trobio_env)

###Prepare data set
###Simulation of plant persistnce for each glacier bases on persistance probabilites. 


listofpred<- list()##store proba for each plant on each glacier
listofproj<-list()##store the projecction for each simulation (per plant and per glacier)
name<-list()##store the name pof plant species 

###This loop prepare an empty data frame to store the simulations
##based on the subset of plant for which we have proba data. 
for (i in 1:length(glacierbrutlist)){   
  SP<-as.data.frame(colnames(glacierbrutlist[[i]][,-1]))
  colnames(SP)<- "sp"
  Pred<-subset(Prob_df, sp %in% SP$sp)
  listofpred[[i]]<-Pred
  sp_name = listofpred[[i]]$sp
  name[[i]]<-sp_name
  proj = data.frame(matrix(data=NA,
                           nrow= 100,
                           ncol=length(name[[i]])))
  
  
  listofproj[[i]]<-proj
  colnames(listofproj[[i]]) = name[[i]]
}


for(j in 1:length(listofproj)){ ##calculate the persistence (presence or absence) of each plant based on the proba.
  for(i in 1:ncol(listofproj[[j]])){
    listofproj[[j]][,i]=
      rbinom(n =100, size = 1, prob = listofpred[[j]]$p_m[i])
    
  }
}

rm(amola, cedec, rutor, trobio, Pred, proj,SP)

#######Step 2#######
#Calculate the stage of each plant species based on simualtion. 
#Each plant can be found on different year following glacier retreat, that is why we want to generalise whether a species
#is pioneer, intermediate or late. 


##sum each plant species, for each year following glacier retreat. 

N_sp<-list()
matbip<-list()


for (i in 1:4){
  
  glacierbrutlist[[i]]<-glacierbrutlist[[i]]%>% remove_rownames %>% column_to_rownames(var="X")
  N_sp[[i]]<- ncol(glacierbrutlist[[i]])
  glacierenvlist[[i]]$years<-as.factor(glacierenvlist[[i]]$years)
  mat = matrix(0, length(levels(glacierenvlist[[i]]$years)), N_sp[[i]])
  matbip[[i]]=mat
  
  for(j in 1:length(levels(glacierenvlist[[i]]$years))){ 
    matbip[[i]][j,] = colSums(glacierbrutlist[[i]][which(glacierenvlist[[i]]$years == levels(glacierenvlist[[i]]$years)[j]),])
    colnames(matbip[[i]]) = colnames(glacierbrutlist[[i]])
    rownames(matbip[[i]]) = levels(glacierenvlist[[i]]$years)
  }
}


##based on the year, generalise the stage of each plant. 

## network matrix
distrnetlist<- list()
commlist<- list()
comm_glacier<- list()
sp_predy_list<-list()
Sitelist<-c("Amola", "Cedec", "Rutor", "Trobio")

for (i in 1:length(matbip)){
  distrnet = graph_from_incidence_matrix(matbip[[i]], weighted = TRUE)
  distrnetlist[[i]] <- distrnet 
  commlist[[i]] = cluster_fast_greedy(distrnetlist[[i]])
  comm_glacier[[i]] = membership(commlist[[i]])[-(1:length(levels(glacierenvlist[[i]]$years)))]
  sp_predy_list[[i]] = data.frame(species = colnames(glacierbrutlist[[i]]),
                                  comm = comm_glacier[[i]],
                                  Site = Sitelist[i])
}

rm(commlist, distrnet, distrnetlist, mat,sp_name)
rm(amola_env, cedec_env, trobio_env, rutor_env)

##########STEP 3 #########
##Link NCP with each glacier

Prob<- Prob_df[,c(1,15,17,19)]
colnames(Prob)[1]<- "species"
Problist<- list()

NCPlist<- list()###Curent NCP list


for (i in 1:4){
  NCPlist[[i]]<- NCP[,-c(1,6,7,12,21:26)]
  colnames(NCPlist[[i]])[1]<- "species"
  NCPlist[[i]]$species<- gsub(" ", ".",NCPlist[[i]]$species)
  NCPlist[[i]]<-merge(NCPlist[[i]], Prob)
  
}

for(i in 1:4){
  Problist[[i]]<- merge(Prob, sp_predy_list[[i]])
  NCPlist[[i]]<- merge(NCPlist[[i]],Problist[[i]])
}

##Merge each simulation with NCPs. 


# Initialize a list to store the results for each data frame
Simulationlist <- list()

# Loop over the four different data frames
for (j in 1:4) {
  # Initialize a list to store the results for the current data frame
  listofdf <- list()
  
  # Loop over the rows (simulations)
  for (i in 1:100) {
    # Select columns where the value is 1
    n <- colnames(listofproj[[j]])[which(listofproj[[j]][i,] == 1)]
    
    # Subset NCP based on selected species
    ncp_proj <- NCPlist[[j]][NCPlist[[j]]$species %in% n,]
    
    # Convert to row names
    ncp_proj <- ncp_proj %>% remove_rownames() %>% column_to_rownames(var = "species")
    
    # Store the result in listofdf
    listofdf[[i]] <- ncp_proj
  }  
  # Store the list of data frames for the current data frame in listoflists
  Simulationlist[[j]] <- listofdf
}



# Iterate over each list in listoflists
for (j in 1:length(Simulationlist)) {
  # Iterate over each data frame in the current list
  for (i in 1:length(Simulationlist[[j]])) {
    # Remove the first 3 columns from the current data frame
    Simulationlist[[j]][[i]] <- Simulationlist[[j]][[i]][, -(1:3)]
  }
}


rm(ncp_proj, NCP, NCP_cat)
rm(Prob, Prob_df, matbip, name)

####Figure NCP early stages comparison with current situation. 

dft<- Simulationlist[[1]][[1]]

Sum1list<- list()


for (i in 1:4){
  
  dfsum1<-data.frame(matrix(nrow = 100, ncol = 15))
  colnames(dfsum1)<-colnames(dft[,c(1:15)])
  
  for (k in 1:100){
    subset_df <- Simulationlist[[i]][[k]]
    column_sums <- colSums(subset_df[subset_df$comm == 1, c(1:15)])
    dfsum1[k,]<-column_sums  
    
    
  }
  Sum1list[[i]]<-dfsum1
}

NCP1list<-list()
for (i in 1:4){
  
  subset_df <- NCPlist[[i]]
  total <- colSums(subset_df[subset_df$comm == 1, c(5:19)])
  
  
  datancp <- data.frame(
    Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
    NCP= c(rep(colnames(Sum1list[[i]]),2)),
    value=as.numeric(cbind(c(colMeans(Sum1list[[i]])), total)),
    max=as.numeric(cbind(c(sapply(Sum1list[[i]], max), total))),
    min=as.numeric(cbind(c(sapply(Sum1list[[i]], min), total)))
  )
  
  
  NCP1list[[i]]<-datancp
}

P_earlyNCP<-list()

for (i in 1:4){
  
  P_Early<-ggplot(data=NCP1list[[i]], aes(x=NCP, y=value, fill=Time)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme( axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
    ggtitle(Sitelist[[i]])+ scale_fill_manual(values=c("skyblue","orange"))+
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9))
  P_earlyNCP[[i]]<-P_Early
  
}

library(ggpubr)

figure_NCP_early<- ggarrange(P_earlyNCP[[1]],P_earlyNCP[[2]], P_earlyNCP[[3]],P_earlyNCP[[4]],
                             labels = letters[1:4],
                             ncol = 2, nrow = 2,
                             common.legend = TRUE, legend = "bottom")  

figure_NCP_early  




####Figure NCP Intermediate stages comparison with current situation. 

dft<- Simulationlist[[1]][[1]]

Sum2list<- list()


for (i in 1:4){
  
  dfsum2<-data.frame(matrix(nrow = 100, ncol = 15)) #create en empty data frame to store the sum of the contribution of plants for each simulation
  colnames(dfsum2)<-colnames(dft[,c(1:15)])
  
  for (k in 1:100){
    subset_df <- Simulationlist[[i]][[k]]
    column_sums <- colSums(subset_df[subset_df$comm == 2, c(1:15)])
    dfsum2[k,]<-column_sums  
    
    
  }
  Sum2list[[i]]<-dfsum2
}

NCP2list<-list()
for (i in 1:4){
  
  subset_df <- NCPlist[[i]]
  total <- colSums(subset_df[subset_df$comm == 2, c(5:19)])
  
  
  datancp <- data.frame(
    Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
    NCP= c(rep(colnames(Sum2list[[i]]),2)),
    value=as.numeric(cbind(c(colMeans(Sum2list[[i]])), total)),
    max=as.numeric(cbind(c(sapply(Sum2list[[i]], max), total))),
    min=as.numeric(cbind(c(sapply(Sum2list[[i]], min), total)))
  )
  
  
  NCP2list[[i]]<-datancp
}

P_IntermediateNCP<-list()

for (i in 1:4){
  
  P_inter<-ggplot(data=NCP2list[[i]], aes(x=NCP, y=value, fill=Time)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme( axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
    ggtitle(Sitelist[[i]])+ scale_fill_manual(values=c("skyblue","orange"))+
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9))
  P_IntermediateNCP[[i]]<-P_inter
  
}



figure_NCP_Inter<- ggarrange(P_IntermediateNCP[[1]],P_IntermediateNCP[[2]], P_IntermediateNCP[[3]],P_IntermediateNCP[[4]],
                             labels = letters[1:4],
                             ncol = 2, nrow = 2,
                             common.legend = TRUE, legend = "bottom")  

figure_NCP_Inter  



####Figure NCP Late stages comparison with current situation. 

dft<- Simulationlist[[1]][[1]]

Sum3list<- list()


for (i in 1:4){
  
  dfsum3<-data.frame(matrix(nrow = 100, ncol = 15)) #create en empty data frame to store the sum of the contribution of plants for each simulation
  colnames(dfsum3)<-colnames(dft[,c(1:15)])
  
  for (k in 1:100){
    subset_df <- Simulationlist[[i]][[k]]
    column_sums <- colSums(subset_df[subset_df$comm == 3, c(1:15)])
    dfsum3[k,]<-column_sums  
    
    
  }
  Sum3list[[i]]<-dfsum3
}

NCP3list<-list()
for (i in 1:4){
  
  subset_df <- NCPlist[[i]]
  total <- colSums(subset_df[subset_df$comm == 1, c(5:19)])
  
  
  datancp <- data.frame(
    Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
    NCP= c(rep(colnames(Sum3list[[i]]),2)),
    value=as.numeric(cbind(c(colMeans(Sum3list[[i]])), total)),
    max=as.numeric(cbind(c(sapply(Sum3list[[i]], max), total))),
    min=as.numeric(cbind(c(sapply(Sum3list[[i]], min), total)))
  )
  
  
  NCP3list[[i]]<-datancp
}

P_LateNCP<-list()

for (i in 1:4){
  
  P_late<-ggplot(data=NCP3list[[i]], aes(x=NCP, y=value, fill=Time)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme( axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
    ggtitle(Sitelist[[i]])+ scale_fill_manual(values=c("skyblue","orange"))+
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9))
  P_LateNCP[[i]]<-P_late
  
}

figure_NCP_Late<- ggarrange(P_LateNCP[[1]],P_LateNCP[[2]], P_LateNCP[[3]],P_LateNCP[[4]],
                            labels = letters[1:4],
                            ncol = 2, nrow = 2,
                            common.legend = TRUE, legend = "bottom")  

figure_NCP_Late  


rm(comm_glacier, datancp)

#######Zscore calculation
dft<- Simulationlist[[1]][[1]]

###Loop for stage 1

Zearly<-list()
for (i in 1:4){
  #create empty dataframe
  z<- data.frame(matrix(nrow = 4, ncol = 15)) 
  rownames(z)<-c("Mean", "SD", "Observed", "Z-Score")
  colnames(z)<- colnames(dft[,c(1:15)])
  
  subset_df <- NCPlist[[i]]
  
  z[1,]<-colMeans(Sum1list[[i]])
  z[2,]<-sqrt(sapply(Sum1list[[i]], var))
  z[3,]<- colSums(subset_df[subset_df$comm == 1, c(5:19)]) 
  z[4,]<- -(z[3,]-z[1,])/z[2,]
  
  Zearly[[i]]<-z
  
}

######zscore intermediate

Zinter<-list()
for (i in 1:4){
  #create empty dataframe
  z<- data.frame(matrix(nrow = 4, ncol = 15)) 
  rownames(z)<-c("Mean", "SD", "Observed", "Z-Score")
  colnames(z)<- colnames(dft[,c(1:15)])
  subset_df <- NCPlist[[i]]
  
  z[1,]<-colMeans(Sum2list[[i]])
  z[2,]<-sqrt(sapply(Sum2list[[i]], var))
  z[3,]<- colSums(subset_df[subset_df$comm == 1, c(5:19)]) 
  z[4,]<- -(z[3,]-z[1,])/z[2,]
  
  Zinter[[i]]<-z
  
}


######zscore Late

Zlate<-list()
for (i in 1:4){
  #create empty dataframe
  z<- data.frame(matrix(nrow = 4, ncol = 15)) 
  rownames(z)<-c("Mean", "SD", "Observed", "Z-Score")
  colnames(z)<- colnames(dft[,c(1:15)])
  subset_df <- NCPlist[[i]]
  
  z[1,]<-colMeans(Sum3list[[i]])
  z[2,]<-sqrt(sapply(Sum3list[[i]], var))
  z[3,]<- colSums(subset_df[subset_df$comm == 1, c(5:19)]) 
  z[4,]<- -(z[3,]-z[1,])/z[2,]
  
  Zlate[[i]]<-z
  
}


rm (z, subset_df)
####p-value calculation
##Early stage

N <- data.frame(matrix(0L, nrow = 101, ncol = 15))

for (k in 1:4) {
  subset_df <- NCPlist[[k]]
  somme <- colSums(subset_df[subset_df$comm == 1, 5:19])
  
  for (j in 1:15) {
    for (i in 1:100) {
      if (Sum1list[[k]][i, j] < somme[j]) 
        N[i, j] <- 1
    }
    # Calculate p-value for each column after updating N
    s <- colSums(N)
    Zearly[[k]][5,j] <- 1 - (s[j] / 101) 
    rownames(Zearly[[k]])[5]<- ("P-value")
    
  }
}


##Intermediate stage

N <- data.frame(matrix(0L, nrow = 101, ncol = 15))

for (k in 1:4) {
  subset_df <- NCPlist[[k]]
  somme <- colSums(subset_df[subset_df$comm == 2, 5:19])
  
  for (j in 1:15) {
    for (i in 1:100) {
      if (Sum2list[[k]][i, j] < somme[j]) 
        N[i, j] <- 1
    }
    # Calculate p-value for each column after updating N
    s <- colSums(N)
    Zinter[[k]][5,j] <- 1 - (s[j] / 101) 
    rownames(Zinter[[k]])[5]<- ("P-value")
    
  }
}

##Late stage

N <- data.frame(matrix(0L, nrow = 101, ncol = 15))
# Calculate the column sums for comm == 1 in subset_df
somme <- colSums(subset_df[subset_df$comm == 3, 5:19])

for (k in 1:4) {
  subset_df <- NCPlist[[k]]
  somme <- colSums(subset_df[subset_df$comm == 1, 5:19])
  
  for (j in 1:15) {
    for (i in 1:100) {
      if (Sum3list[[k]][i, j] < somme[j]) 
        N[i, j] <- 1
    }
    # Calculate p-value for each column after updating N
    s <- colSums(N)
    Zlate[[k]][5,j] <- 1 - (s[j] / 101) 
    rownames(Zlate[[k]])[5]<- ("P-value")
    
  }
}






#########Figure##########

##Figure plant species
library(tidyverse)

sim1list<-list() 
sim2list<-list()
sim3list<-list()

for (i in 1:4){
  
  sim1<-listofproj[[i]] %>% select(any_of(sp_predy_list[[i]]$species[sp_predy_list[[i]]$comm ==1]))
  sim1list[[i]]<-sim1
  
  sim2<-listofproj[[i]] %>% select(any_of(sp_predy_list[[i]]$species[sp_predy_list[[i]]$comm ==2]))
  sim2list[[i]]<-sim2
  
  sim3<-listofproj[[i]] %>% select(any_of(sp_predy_list[[i]]$species[sp_predy_list[[i]]$comm ==3]))
  sim3list[[i]]<-sim3
  
  
}

Speciescomp<- list()

for (i in 1:4){
  Comp <- data.frame(
    Time=c(rep("Future - no glacier" , 4) , rep("Current - 2020" , 4) ),
    Stage=c(rep(c("Early", "Intermediate", "Late", "Total"), 2)),
    value=c(mean(rowSums(sim1list[[i]])),mean(rowSums(sim2list[[i]])),mean(rowSums(sim3list[[i]])),mean(rowSums(listofproj[[i]])), length(sim1list[[i]]), length(sim2list[[i]]), length(sim3list[[i]]), nrow(NCPlist[[i]])),
    max=c(max(rowSums(sim1list[[i]])), max(rowSums(sim2list[[i]])), max(rowSums(sim3list[[i]])), max(rowSums(listofproj[[i]])),length(sim1list[[i]]), length(sim2list[[i]]), length(sim3list[[i]]), nrow(NCPlist[[i]])),
    min=c( min(rowSums(sim1list[[i]])), min(rowSums(sim2list[[i]])), min(rowSums(sim3list[[i]])),min(rowSums(listofproj[[i]])), length(sim1list[[i]]), length(sim2list[[i]]), length(sim3list[[i]]), nrow(NCPlist[[i]])),
    name = letters[1:8]
  )
  Speciescomp[[i]]<-Comp
} 

Plotspecies <- list()

for (i in 1:4) {
  print(paste("Iteration:", i))
  print(Speciescomp[[i]])  # Print out the data for debugging
  
  P_species <- ggplot(data = Speciescomp[[i]], aes(x = Stage, y = value, fill = Time)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    theme(panel.background = element_rect(fill = "white")) +
    xlab("Stage") + ylab("Number of Species") +
    ggtitle(Sitelist[i]) + 
    scale_fill_manual(values = c("skyblue", "orange")) +
    geom_errorbar(aes(ymin = min, ymax = max), width = 0.2, position = position_dodge(0.9))
  
  Plotspecies[[i]] <- P_species
  
  print(Plotspecies[[i]])  # Print the plot for debugging
}

library(ggpubr)

figure_Plant<- ggarrange(Plotspecies[[1]],Plotspecies[[2]], Plotspecies[[3]],Plotspecies[[4]],
                         labels = letters[1:4],
                         ncol = 2, nrow = 2,
                         common.legend = TRUE, legend = "bottom")  

figure_Plant  

rm(Comp, sim1, sim2, sim3, n)


#####Figure Zscore

##early 

datazearly<- list()

for (i in 1:4){
  datazearly[[i]]<- data.frame(
    NCP=colnames(Zearly[[i]]),
    Zscore=as.numeric(Zearly[[i]][4,]),
    Pvalue=as.numeric(Zearly[[i]][5,])
  ) 
}

Pzearly<-list()

for (i in 1:4){
  
  Pzearly[[i]]<-ggplot(datazearly[[i]], aes(x=NCP, y=Zscore)) +
    geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Zscore")+ ggtitle(Sitelist[i])
}
print(Pzearly)


figure_Pzearly<- ggarrange(Pzearly[[1]],Pzearly[[2]], Pzearly[[3]],Pzearly[[4]],
                           labels = letters[1:4],
                           ncol = 2, nrow = 2,
                           common.legend = TRUE, legend = "bottom")  

figure_Pzearly 


###intermediate
datazinter<- list()

for (i in 1:4){
  datazinter[[i]]<- data.frame(
    NCP=colnames(Zinter[[i]]),
    Zscore=as.numeric(Zinter[[i]][4,]),
    Pvalue=as.numeric(Zinter[[i]][5,])
  ) 
}

Pzinter<-list()

for (i in 1:4){
  
  Pzinter[[i]]<-ggplot(datazinter[[i]], aes(x=NCP, y=Zscore)) +
    geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Zscore")+ ggtitle(Sitelist[i])
}
print(Pzinter)


figure_Pzinter<- ggarrange(Pzinter[[1]],Pzinter[[2]], Pzinter[[3]],Pzinter[[4]],
                           labels = letters[1:4],
                           ncol = 2, nrow = 2,
                           common.legend = TRUE, legend = "bottom")  

figure_Pzinter 


###Late
datazlate<- list()

for (i in 1:4){
  datazlate[[i]]<- data.frame(
    NCP=colnames(Zinter[[i]]),
    Zscore=as.numeric(Zlate[[i]][4,]),
    Pvalue=as.numeric(Zlate[[i]][5,])
  ) 
}

Pzlate<-list()

for (i in 1:4){
  
  Pzlate[[i]]<-ggplot(datazlate[[i]], aes(x=NCP, y=Zscore)) +
    geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Zscore")+ ggtitle(Sitelist[i])
}
print(Pzlate)


figure_Pzlate<- ggarrange(Pzlate[[1]],Pzlate[[2]], Pzlate[[3]],Pzlate[[4]],
                           labels = letters[1:4],
                           ncol = 2, nrow = 2,
                           common.legend = TRUE, legend = "bottom")  

figure_Pzlate 



############Summary #################

NCP1<- NCPlist[[1]]
NCP2<- NCPlist[[2]]
NCP3<- NCPlist[[3]]
NCP4<- NCPlist[[4]]

###Stage 1 : 
Early_ncp_future<-rbind(Sum1list[[1]], Sum1list[[2]], Sum1list[[3]],Sum1list[[4]])
Early_ncp_current<-rbind(colSums(NCP1[NCP1$comm == 1, c(5:19)]),colSums(NCP2[NCP2$comm == 1, c(5:19)]), colSums(NCP3[NCP3$comm == 1, c(5:19)]), colSums(NCP4[NCP4$comm == 1, c(5:19)]))
Early_ncp_current<-as.data.frame(Early_ncp_current)
Early_ncp_current[5,]<- sapply(Early_ncp_current, mean)   
Early_ncp_current[6,]<- sapply(Early_ncp_current[1:4,], sd)


summaryncp_1 <- data.frame(
  Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
  NCP= c(rep(colnames(Early_ncp_future),2)),
  value=as.numeric(cbind(c(colMeans(Early_ncp_future), Early_ncp_current[5,]))),
  max=as.numeric(cbind(c(sapply(Early_ncp_future, max), Early_ncp_current[5,]+Early_ncp_current[6,]))),
  min=as.numeric(cbind(c(sapply(Early_ncp_future, min), Early_ncp_current[5,]-Early_ncp_current[6,])))
)


P_summary_early<-ggplot(data=summaryncp_1, aes(x=NCP, y=value, fill=Time)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme( axis.text.x = element_text(hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
  ggtitle("NCP - Early summary")+ scale_fill_manual(values=c("skyblue","orange"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9)) +  scale_x_discrete(labels=labels) +
  coord_flip()

P_summary_early

###Stage 2 : 
Inter_ncp_future<-rbind(Sum2list[[1]], Sum2list[[2]], Sum2list[[3]],Sum2list[[4]])
Inter_ncp_current<-rbind(colSums(NCP1[NCP1$comm == 2, c(5:19)]),colSums(NCP2[NCP2$comm == 2, c(5:19)]), colSums(NCP3[NCP3$comm == 2, c(5:19)]), colSums(NCP4[NCP4$comm == 2, c(5:19)]))
Inter_ncp_current<-as.data.frame(Inter_ncp_current)
Inter_ncp_current[5,]<- sapply(Inter_ncp_current, mean)   
Inter_ncp_current[6,]<- sapply(Inter_ncp_current[1:4,], sd)


summaryncp_2 <- data.frame(
  Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
  NCP= c(rep(colnames(Inter_ncp_future),2)),
  value=as.numeric(cbind(c(colMeans(Inter_ncp_future), Inter_ncp_current[5,]))),
  max=as.numeric(cbind(c(sapply(Inter_ncp_future, max), Inter_ncp_current[5,]+Inter_ncp_current[6,]))),
  min=as.numeric(cbind(c(sapply(Inter_ncp_future, min), Inter_ncp_current[5,]-Inter_ncp_current[6,])))
)


P_summary_Inter<-ggplot(data=summaryncp_2, aes(x=NCP, y=value, fill=Time)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme( axis.text.x = element_text( hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
  ggtitle("NCP - Intermediate summary")+ scale_fill_manual(values=c("skyblue","orange"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9))+  scale_x_discrete(labels=labels) + 
  coord_flip()

P_summary_Inter


###Stage 3 : 
Late_ncp_future<-rbind(Sum3list[[1]], Sum3list[[2]], Sum3list[[3]],Sum3list[[4]])
Late_ncp_current<-rbind(colSums(NCP1[NCP1$comm == 3, c(5:19)]),colSums(NCP2[NCP2$comm == 3, c(5:19)]), colSums(NCP3[NCP3$comm == 3, c(5:19)]), colSums(NCP4[NCP4$comm == 3, c(5:19)]))
Late_ncp_current<-as.data.frame(Late_ncp_current)
Late_ncp_current[5,]<- sapply(Late_ncp_current, mean)   
Late_ncp_current[6,]<- sapply(Late_ncp_current[1:4,], sd)


summaryncp_3 <- data.frame(
  Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
  NCP= c(rep(colnames(Inter_ncp_future),2)),
  value=as.numeric(cbind(c(colMeans(Late_ncp_future), Late_ncp_current[5,]))),
  max=as.numeric(cbind(c(sapply(Late_ncp_future, max), Late_ncp_current[5,]+Late_ncp_current[6,]))),
  min=as.numeric(cbind(c(sapply(Late_ncp_future, min), Late_ncp_current[5,]-Late_ncp_current[6,])))
)


P_summary_late<-ggplot(data=summaryncp_3, aes(x=NCP, y=value, fill=Time)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme( axis.text.x = element_text(hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
  ggtitle("NCP - Late summary")+ scale_fill_manual(values=c("skyblue","orange"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9))+  scale_x_discrete(labels=labels)+
  coord_flip()

P_summary_late


figure_NCP_summary<- ggarrange(P_summary_early,P_summary_Inter,P_summary_late,
                            labels = letters[1:3],
                            ncol = 1, nrow = 3,
                            common.legend = TRUE, legend = "bottom")  

figure_NCP_summary 


####zscore 



##1. Calculate Mean of simulation


#z=current-meanfutur/sdfut

z_early<--(Early_ncp_current[5,] - sapply(Early_ncp_future, function(x) mean(x, na.rm = TRUE)) /
             sapply(Early_ncp_future, function(x) sd(x, na.rm = TRUE)))

z_inter<- -(Inter_ncp_current[5,]-sapply(Inter_ncp_future, function(x) mean(x, na.rm = TRUE))/
              sapply(Inter_ncp_future, function(x) sd(x, na.rm = TRUE)))
z_late <- -(Late_ncp_current[5,]-sapply(Late_ncp_future, function(x) mean(x, na.rm = TRUE))/
              sapply(Late_ncp_future, function(x) sd(x, na.rm = TRUE)))

z_late[] <- lapply(z_late, function(x) {
  x[is.nan(x)] <- 0
  return(x)
})

str(z_late)
####Pvalue 
N<- data.frame(matrix(0L,nrow = 400, ncol = 15))

Early_Ncp.2<- rbind(Early_ncp_current[5,], Early_ncp_future) 

#Stage 1

# Initialiser la matrice N
N <- matrix(0, nrow = 401, ncol = 15)

# Boucle pour parcourir les colonnes
for (j in 1:15) {
  # Boucle pour parcourir les lignes
  for (i in 2:401) {
    # Vérifier si la valeur est non-NA
    if (!is.na(Early_Ncp.2[i, j])) {
      # Vérifier la condition
      if (Early_Ncp.2[i, j] < Early_Ncp.2[1, j]) {
        N[i, j] <- 1
      }
    }
  }
}

# Calculer la somme des colonnes de N
s <- colSums(N)

# Calculer z_early pour la deuxième ligne
z_early[2, ] <- 1 - (s / 401)
rownames(z_early)[2]<- ("P-value")
str(z_early)






# Initialiser la matrice N avec des zéros
N <- matrix(0, nrow = 401, ncol = 15)

# Initialiser z_early avec des zéros (ou une matrice de la taille appropriée)


# Boucle pour parcourir les colonnes
for (j in 1:15) {
  # Boucle pour parcourir les lignes
  for (i in 2:401) {
    # Vérifier si la valeur est non-NA
    if (!is.na(Early_Ncp.2[i, j]) && !is.na(Early_Ncp.2[1, j])) {
      # Vérifier la condition
      if (Early_Ncp.2[i, j] < Early_Ncp.2[1, j]) {
        N[i, j] <- 1
      }
    }
  }
  # Calculer la somme des colonnes de N après avoir parcouru toutes les lignes
  s <- colSums(N)
  # Calculer z_early pour la deuxième ligne
  z_early[2, j] <- 1 - (s[j] / 401)
}

# Afficher le résultat
print(z_early)





#Stage 2

Inter_Ncp.2<- rbind(Inter_ncp_current[5,], Inter_ncp_future) 

for (j in 1:15){
  for (i in 2:401){
    if (Inter_Ncp.2[i,j]< Inter_Ncp.2[1,j]) N[i,j]<-1
    s<-colSums(N)
    z_inter[2,j]<- 1-(s[j]/401)
  }
}
rownames(z_inter)[2]<- ("P-value")

#Stage 3

Late_Ncp.2<- rbind(Late_ncp_current[5,], Late_ncp_future) 

for (j in 1:15){
  for (i in 2:401){
    if (Late_Ncp.2[i,j]< Late_Ncp.2[1,j]) N[i,j]<-1
    s<-colSums(N)
    z_late[2,j]<- 1-(s[j]/401)
  }
}
rownames(z_late)[2]<- ("P-value")


####Zscore plot 

###Early
z<- z_early

dataz<- data.frame(
  NCP=colnames(z),
  Zscore=as.numeric(z[1,]),
  Pvalue=as.numeric(z[2,])
)  


P_zscore_early<-ggplot(dataz, aes(x=NCP, y=Zscore)) +
  geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Zscore")+ ggtitle("Zscore - Early stage")+  scale_x_discrete(labels=labels)

P_zscore_early

###Intermediate 

z<- z_inter

dataz<- data.frame(
  NCP=colnames(z),
  Zscore=as.numeric(z[1,]),
  Pvalue=as.numeric(z[2,])
)  


P_zscore_inter<-ggplot(dataz, aes(x=NCP, y=Zscore)) +
  geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Zscore")+ ggtitle("Zscore - Intermediate stage")+  scale_x_discrete(labels=labels)

P_zscore_inter
###Late
z<- z_late

dataz<- data.frame(
  NCP=colnames(z),
  Zscore=as.numeric(z[1,]),
  Pvalue=as.numeric(z[2,])
)  


P_zscore_late<-ggplot(dataz, aes(x=NCP, y=Zscore)) +
  geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Zscore")+ ggtitle("Zscore - Late stage")+  scale_x_discrete(labels=labels)

P_zscore_late


figure_zscore_summary<- ggarrange(P_zscore_early,P_zscore_inter,P_zscore_late,
                               labels = letters[1:3],
                               ncol = 3, nrow = 1,
                               common.legend = TRUE, legend = "bottom")  

figure_zscore_summary 



#####################RATIO calculation and alyses###########
#Here instead of talking the raw number of plants we use the NCP / by the nzumber of plant species
#we cans ee ft the contribution each NCP changes or not. 


nspecies <- rep(NA, 100)


###Early  

listofnspecies1<-list()##calculate the number of species per simulation
for (i in 1:4){
  for (k in 1:100){
    subset_df <- Simulationlist[[i]][[k]]
    nsp <- length(subset_df[subset_df$comm == 1, 1])
    nspecies[k]<-nsp  
  }
  listofnspecies1[[i]]<- nspecies
}

ratioearly<-list() ##calculate the tation
dfratio1 <- as.data.frame(matrix(NA, nrow = 101, ncol = 15)) #df combinuing the current ratio and the simulation
colnames(dfratio1)<-colnames(dft[,c(1:15)])
rownames(dfratio1)[1]<- "now"

for (i in 1:4){
  NCP<- NCPlist[[i]]
  dfratio1[1, ]<-colSums(NCP[NCP$comm == 1, c(5:19)])/length(((NCP[NCP$comm == 1, 1])  )) ##current ratio
  dfratio1[2:101, ]<-Sum1list[[i]]/listofnspecies1[[i]]
  ratioearly[[i]]<- dfratio1
}


##intermediate

listofnspecies2<-list()
for (i in 1:4){
  for (k in 1:100){
    subset_df <- Simulationlist[[i]][[k]]
    nsp <- length(subset_df[subset_df$comm == 2, 1])
    nspecies[k]<-nsp  
  }
  listofnspecies2[[i]]<- nspecies
}


ratiointer<<-list()
dfratio2 <- as.data.frame(matrix(NA, nrow = 101, ncol = 15))
colnames(dfratio2)<-colnames(dft[,c(1:15)])
rownames(dfratio2)[1]<- "now"

for (i in 1:4){
  NCP<- NCPlist[[i]]
  dfratio2[1, ]<-colSums(NCP[NCP$comm == 2, c(5:19)])/length(((NCP[NCP$comm == 2, 1])  )) ##current ratio
  dfratio2[2:101, ]<-Sum2list[[i]]/listofnspecies2[[i]]
  ratiointer[[i]]<- dfratio2
}


##late

listofnspecies3<-list()
for (i in 1:4){
  for (k in 1:100){
    subset_df <- Simulationlist[[i]][[k]]
    nsp <- length(subset_df[subset_df$comm == 3, 1])
    nspecies[k]<-nsp  
  }
  listofnspecies3[[i]]<- nspecies
}


ratiolate<-list()
dfratio3 <- as.data.frame(matrix(NA, nrow = 101, ncol = 15))
colnames(dfratio3)<-colnames(dft[,c(1:15)])
rownames(dfratio3)[1]<- "now"

for (i in 1:4){
  NCP<- NCPlist[[i]]
  dfratio3[1, ]<-colSums(NCP[NCP$comm == 3, c(5:19)])/length(((NCP[NCP$comm == 3, 1])  )) ##current ratio
  dfratio3[2:101, ]<-Sum3list[[i]]/listofnspecies3[[i]]
  ratiolate[[i]]<- dfratio3
}


####figure#####

SummaryEarly_ratio_list<-list()

for (i in 1:4){
  NCP_val<-ratioearly[[i]]
  
  datancp <- data.frame(
    Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
    NCP= c(rep(colnames(NCP_val),2)),
    value=as.numeric(cbind(c(colMeans(NCP_val[-1,], na.rm =TRUE), NCP_val[1,]))),
    max=as.numeric(cbind(c(sapply(NCP_val[-1,], max, na.rm =TRUE), NCP_val[1,]))),
    min=as.numeric(cbind(c(sapply(NCP_val[-1,], min,na.rm =TRUE), NCP_val[1,])))
  )
  SummaryEarly_ratio_list[[i]]<-datancp
}


P_NCP_Early_ratio<-list()
for (i in 1:4){
  P_NCP_Early_ratio[[i]]<-ggplot(data=SummaryEarly_ratio_list[[i]], aes(x=NCP, y=value, fill=Time)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme( axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
    ggtitle(Sitelist[i])+ scale_fill_manual(values=c("skyblue","orange"))+
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9))
  
}

figure_P_NCP_Early_ratio<- ggarrange(P_NCP_Early_ratio[[1]],P_NCP_Early_ratio[[2]], P_NCP_Early_ratio[[3]],P_NCP_Early_ratio[[4]],
                                     labels = letters[1:4],
                                     ncol = 2, nrow = 2,
                                     common.legend = TRUE, legend = "bottom")  
figure_P_NCP_Early_ratio



###intermediate
Summaryinter_ratio_list<-list()


for (i in 1:4){
  NCP_val<-ratiointer[[i]]
  
  datancp <- data.frame(
    Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
    NCP= c(rep(colnames(NCP_val),2)),
    value=as.numeric(cbind(c(colMeans(NCP_val[-1,], na.rm =TRUE), NCP_val[1,]))),
    max=as.numeric(cbind(c(sapply(NCP_val[-1,], max, na.rm =TRUE), NCP_val[1,]))),
    min=as.numeric(cbind(c(sapply(NCP_val[-1,], min,na.rm =TRUE), NCP_val[1,])))
  )
  Summaryinter_ratio_list[[i]]<-datancp
}


P_NCP_inter_ratio<-list()
for (i in 1:4){
  P_NCP_inter_ratio[[i]]<-ggplot(data=Summaryinter_ratio_list[[i]], aes(x=NCP, y=value, fill=Time)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme( axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
    ggtitle(Sitelist[i])+ scale_fill_manual(values=c("skyblue","orange"))+
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9))
  
}


figure_P_NCP_inter_ratio<- ggarrange(P_NCP_inter_ratio[[1]],P_NCP_inter_ratio[[2]], P_NCP_inter_ratio[[3]],P_NCP_inter_ratio[[4]],
                                     labels = letters[1:4],
                                     ncol = 2, nrow = 2,
                                     common.legend = TRUE, legend = "bottom")  
figure_P_NCP_inter_ratio


###late

Summarylate_ratio_list<-list()

for (i in 1:4){
  NCP_val<-ratiolate[[i]]
  
  datancp <- data.frame(
    Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
    NCP= c(rep(colnames(NCP_val),2)),
    value=as.numeric(cbind(c(colMeans(NCP_val[-1,], na.rm =TRUE), NCP_val[1,]))),
    max=as.numeric(cbind(c(sapply(NCP_val[-1,], max, na.rm =TRUE), NCP_val[1,]))),
    min=as.numeric(cbind(c(sapply(NCP_val[-1,], min,na.rm =TRUE), NCP_val[1,])))
  )
  Summarylate_ratio_list[[i]]<-datancp
}


P_NCP_late_ratio<-list()
for (i in 1:4){
  P_NCP_late_ratio[[i]]<-ggplot(data=Summarylate_ratio_list[[i]], aes(x=NCP, y=value, fill=Time)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme( axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
    ggtitle(Sitelist[i])+ scale_fill_manual(values=c("skyblue","orange"))+
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9))
  
}

figure_P_NCP_late_ratio<- ggarrange(P_NCP_late_ratio[[1]],P_NCP_late_ratio[[2]], P_NCP_late_ratio[[3]],P_NCP_late_ratio[[4]],
                                    labels = letters[1:4],
                                    ncol = 2, nrow = 2,
                                    common.legend = TRUE, legend = "bottom")  
figure_P_NCP_late_ratio


#######Zscore calculation


###Loop for stage 1

Zearlyratio<-list()
for (i in 1:4){
  #create empty dataframe
  z<- data.frame(matrix(nrow = 4, ncol = 15)) 
  rownames(z)<-c("Mean", "SD", "Observed", "Z-Score")
  colnames(z)<- colnames(dft[,c(1:15)])
  
  z[1,]<-colMeans(ratioearly[[i]], na.rm = TRUE)
  z[2,]<-sqrt(sapply(ratioearly[[i]], var, na.rm = TRUE))
  z[3,]<- ratioearly[[1]][1,c(1:15)]
  z[4,]<- -(z[3,]-z[1,])/z[2,]
  
  Zearlyratio[[i]]<-z
  
}


#######Zscore calculation


###Loop for stage 2

Zintereratio<-list()
for (i in 1:4){
  #create empty dataframe
  z<- data.frame(matrix(nrow = 4, ncol = 15)) 
  rownames(z)<-c("Mean", "SD", "Observed", "Z-Score")
  colnames(z)<- colnames(dft[,c(1:15)])
  
  z[1,]<-colMeans(ratiointer[[i]], na.rm = TRUE)
  z[2,]<-sqrt(sapply(ratiointer[[i]], var, na.rm = TRUE))
  z[3,]<- ratiointer[[1]][1,c(1:15)]
  z[4,]<- -(z[3,]-z[1,])/z[2,]
  
  Zintereratio[[i]]<-z
  
}



zlateratio<-list()
for (i in 1:4){
  #create empty dataframe
  z<- data.frame(matrix(nrow = 4, ncol = 15)) 
  rownames(z)<-c("Mean", "SD", "Observed", "Z-Score")
  colnames(z)<- colnames(dft[,c(1:15)])
  
  z[1,]<-colMeans(ratiolate[[i]], na.rm = TRUE)
  z[2,]<-sqrt(sapply(ratiolate[[i]], var, na.rm = TRUE))
  z[3,]<- ratiolate[[1]][1,c(1:15)]
  z[4,]<- -(z[3,]-z[1,])/z[2,]
  
  zlateratio[[i]]<-z
  
}





N <- data.frame(matrix(0L, nrow = 101, ncol = 15))
# Calculate the column sums for comm == 1 in subset_df
somme <- ratioearly[[k]][1,]



for (k in 1:4) {
  somme <- ratioearly[[k]][1,]
  for (j in 1:15) {
    for (i in 2:101) {
      if (!is.na(ratioearly[[k]][i, j]) && !is.na(somme[j]) && ratioearly[[k]][i, j] < somme[j]) {
        N[i, j] <- 1
      }
    }
    # Calculate p-value for each column after updating N
    s <- colSums(N)
    Zearlyratio[[k]][5,j] <- 1 - (s[j] / 101) 
    rownames(Zearlyratio[[k]])[5]<- ("P-value")
    
  }
}


# Calculate the column sums for comm == 1 in subset_df


for (k in 1:4) {
  somme <- ratiointer[[k]][1,]
  for (j in 1:15) {
    for (i in 2:101) {
      if (!is.na(ratiointer[[k]][i, j]) && !is.na(somme[j]) && ratiointer[[k]][i, j] < somme[j]) {
        N[i, j] <- 1
      }
    }
    # Calculate p-value for each column after updating N
    s <- colSums(N)
    Zintereratio[[k]][5,j] <- 1 - (s[j] / 101) 
    rownames(Zintereratio[[k]])[5]<- ("P-value")
  }
}



for (k in 1:4) {
  somme <- ratiolate[[k]][1,]
  for (j in 1:15) {
    for (i in 2:101) {
      if (!is.na(ratiolate[[k]][i, j]) && !is.na(somme[j]) && ratiolate[[k]][i, j] < somme[j]) {
        N[i, j] <- 1
      }
    }
    # Calculate p-value for each column after updating N
    s <- colSums(N)
    zlateratio[[k]][5,j] <- 1 - (s[j] / 101) 
    rownames(zlateratio[[k]])[5]<- ("P-value")
  }
}

#####Figure Zscore

##early 

datazearlyratio<- list()

for (i in 1:4){
  datazearlyratio[[i]]<- data.frame(
    NCP=colnames(Zearlyratio[[i]]),
    Zscore=as.numeric(Zearlyratio[[i]][4,]),
    Pvalue=as.numeric(Zearlyratio[[i]][5,])
  ) 
}

Pzearlyratio<-list()

for (i in 1:4){
  
  Pzearlyratio[[i]]<-ggplot(datazearlyratio[[i]], aes(x=NCP, y=Zscore)) +
    geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Zscore")+ ggtitle(Sitelist[i])
}
print(Pzearlyratio)


figure_Pzearlyratio<- ggarrange(Pzearlyratio[[1]],Pzearlyratio[[2]], Pzearlyratio[[3]],Pzearlyratio[[4]],
                                labels = letters[1:4],
                                ncol = 2, nrow = 2,
                                common.legend = TRUE, legend = "bottom")  

figure_Pzearlyratio 

##Intermediate 

datazinterratio<- list()

for (i in 1:4){
  datazinterratio[[i]]<- data.frame(
    NCP=colnames(Zintereratio[[i]]),
    Zscore=as.numeric(Zintereratio[[i]][4,]),
    Pvalue=as.numeric(Zintereratio[[i]][5,])
  ) 
}

Pzinterratio<-list()

for (i in 1:4){
  
  Pzinterratio[[i]]<-ggplot(datazinterratio[[i]], aes(x=NCP, y=Zscore)) +
    geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Zscore")+ ggtitle(Sitelist[i])
}
print(Pzinterratio)


figure_Pzinterratio<- ggarrange(Pzinterratio[[1]],Pzinterratio[[2]], Pzinterratio[[3]],Pzinterratio[[4]],
                                labels = letters[1:4],
                                ncol = 2, nrow = 2,
                                common.legend = TRUE, legend = "bottom")  

figure_Pzinterratio 


##late 

datazlateratio<- list()

for (i in 1:4){
  datazlateratio[[i]]<- data.frame(
    NCP=colnames(zlateratio[[i]]),
    Zscore=as.numeric(zlateratio[[i]][4,]),
    Pvalue=as.numeric(zlateratio[[i]][5,])
  ) 
}

Pzlateratio<-list()

for (i in 1:4){
  
  Pzlateratio[[i]]<-ggplot(datazlateratio[[i]], aes(x=NCP, y=Zscore)) +
    geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("") +
    ylab("Zscore")+ ggtitle(Sitelist[i])
}
print(Pzlateratio)


figure_Pzlateratio<- ggarrange(Pzlateratio[[1]],Pzlateratio[[2]], Pzlateratio[[3]],Pzlateratio[[4]],
                               labels = letters[1:4],
                               ncol = 2, nrow = 2,
                               common.legend = TRUE, legend = "bottom")  

figure_Pzlateratio



#############Summary ratio##############

####figure NCP ratio



###Stage 1 : 
Early_ncp_futureratio<-rbind(ratioearly[[1]][c(2:101),], ratioearly[[2]][c(2:101),], ratioearly[[3]][c(2:101),],ratioearly[[4]][c(2:101),])
Early_ncp_currentratio<-rbind(ratioearly[[1]][1,], ratioearly[[2]][1,], ratioearly[[3]][1,],ratioearly[[4]][1,])
Early_ncp_currentratio<-as.data.frame(Early_ncp_currentratio)
Early_ncp_currentratio[5,]<- sapply(Early_ncp_currentratio, mean)   
Early_ncp_currentratio[6,]<- sapply(Early_ncp_currentratio[1:4,], sd)


summaryncpratio_1 <- data.frame(
  Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
  NCP= c(rep(colnames(Early_ncp_futureratio),2)),
  value=as.numeric(cbind(c(colMeans(Early_ncp_futureratio, na.rm = TRUE), Early_ncp_currentratio[5,]))),
  max=as.numeric(cbind(c(sapply(Early_ncp_futureratio, max, na.rm = TRUE), Early_ncp_currentratio[5,]+Early_ncp_currentratio[6,]))),
  min=as.numeric(cbind(c(sapply(Early_ncp_futureratio, min, na.rm = TRUE), Early_ncp_currentratio[5,]-Early_ncp_currentratio[6,])))
)


P_summary_early_ratio<-ggplot(data=summaryncpratio_1, aes(x=NCP, y=value, fill=Time)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme( axis.text.x = element_text( hjust = 0.5),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
  ggtitle("NCP - Early summary - Ratio")+ scale_fill_manual(values=c("skyblue","orange"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9))+  scale_x_discrete(labels=labels)+
  coord_flip()

P_summary_early_ratio




###Stage 2 : 
Inter_ncp_futureratio<-rbind(ratiointer[[1]][c(2:101),], ratiointer[[2]][c(2:101),], ratiointer[[3]][c(2:101),],ratiointer[[4]][c(2:101),])
Inter_ncp_currentratio<-rbind(ratiointer[[1]][1,], ratiointer[[2]][1,], ratiointer[[3]][1,],ratiointer[[4]][1,])
Inter_ncp_currentratio<-as.data.frame(Inter_ncp_currentratio)
Inter_ncp_currentratio[5,]<- sapply(Inter_ncp_currentratio, mean)   
Inter_ncp_currentratio[6,]<- sapply(Inter_ncp_currentratio[1:4,], sd)


summaryncpratio_2 <- data.frame(
  Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
  NCP= c(rep(colnames(Inter_ncp_futureratio),2)),
  value=as.numeric(cbind(c(colMeans(Inter_ncp_futureratio, na.rm = TRUE), Inter_ncp_currentratio[5,]))),
  max=as.numeric(cbind(c(sapply(Inter_ncp_futureratio, max, na.rm = TRUE), Inter_ncp_currentratio[5,]+Inter_ncp_currentratio[6,]))),
  min=as.numeric(cbind(c(sapply(Inter_ncp_futureratio, min, na.rm = TRUE), Inter_ncp_currentratio[5,]-Inter_ncp_currentratio[6,])))
)


P_summary_inter_ratio<-ggplot(data=summaryncpratio_2, aes(x=NCP, y=value, fill=Time)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme( axis.text.x = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
  ggtitle("NCP - Intermediate summary - Ratio")+ scale_fill_manual(values=c("skyblue","orange"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9))+  scale_x_discrete(labels=labels)+
  coord_flip()

P_summary_inter_ratio


###Stage 3 : 
Late_ncp_futureratio<-rbind(ratiolate[[1]][c(2:101),], ratiolate[[2]][c(2:101),], ratiolate[[3]][c(2:101),],ratiolate[[4]][c(2:101),])
Late_ncp_currentratio<-rbind(ratiolate[[1]][1,], ratiolate[[2]][1,], ratiolate[[3]][1,],ratiolate[[4]][1,])
Late_ncp_currentratio<-as.data.frame(Late_ncp_currentratio)
Late_ncp_currentratio[5,]<- sapply(Late_ncp_currentratio, mean)   
Late_ncp_currentratio[6,]<- sapply(Late_ncp_currentratio[1:4,], sd)


summaryncpratio_3 <- data.frame(
  Time=c(rep("Future - No glacier" , 15) , rep("Current - 2020" , 15) ),
  NCP= c(rep(colnames(Late_ncp_futureratio),2)),
  value=as.numeric(cbind(c(colMeans(Late_ncp_futureratio, na.rm = TRUE), Late_ncp_currentratio[5,]))),
  max=as.numeric(cbind(c(sapply(Late_ncp_futureratio, max, na.rm = TRUE), Late_ncp_currentratio[5,]+Late_ncp_currentratio[6,]))),
  min=as.numeric(cbind(c(sapply(Late_ncp_futureratio, min, na.rm = TRUE), Late_ncp_currentratio[5,]-Late_ncp_currentratio[6,])))
)


P_summary_late_ratio<-ggplot(data=summaryncpratio_3, aes(x=NCP, y=value, fill=Time)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme( axis.text.x = element_text( hjust = 0.5),panel.background = element_rect(fill = "white"))+xlab("NCP")+ylab("Number of Species")+
  ggtitle("NCP - Late summary - Ratio")+ scale_fill_manual(values=c("skyblue","orange"))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9))+  scale_x_discrete(labels=labels)+
  coord_flip()

P_summary_late_ratio



figure_ratio_summary<- ggarrange(P_summary_early_ratio, P_summary_inter_ratio,P_summary_late_ratio,
                               labels = letters[1:3],
                               ncol = 1, nrow = 3,
                               common.legend = TRUE, legend = "bottom")  

figure_ratio_summary





####zscore 


#z=current-meanfutur/sdfut

z_earlyratiosummary<--(Early_ncp_currentratio[5,]-sapply(Early_ncp_futureratio, mean, na.rm =TRUE)/
                         sapply(Early_ncp_futureratio, sd, na.rm =TRUE))

z_interratiosummary<- -(Inter_ncp_currentratio[5,]-sapply(Inter_ncp_futureratio, mean, na.rm =TRUE)/
                          sapply(Inter_ncp_futureratio, sd, na.rm =TRUE))
zlateratiosummary <- -(Late_ncp_currentratio[5,]-sapply(Late_ncp_futureratio, mean, na.rm =TRUE)/
                         sapply(Late_ncp_futureratio, sd, na.rm =TRUE))



####Pvalue 
N<- data.frame(matrix(0L,nrow = 400, ncol = 15))

Early_Ncp.2ratio<- rbind(Early_ncp_currentratio[5,], Early_ncp_futureratio) 

#Stage 1

# Initialize N matrix
N <- matrix(0, nrow = 401, ncol = 15)

# Loop through each column
for (j in 1:15) {
  # Loop through each row starting from the second row
  for (i in 2:401) {
    # Check if both values are not missing before comparison
    if (!is.na(Early_Ncp.2ratio[i, j]) && !is.na(Early_Ncp.2ratio[1, j]) && Early_Ncp.2ratio[i, j] < Early_Ncp.2ratio[1, j]) {
      N[i, j] <- 1
    }
  }
  # Calculate p-value for each column after updating N
  s <- colSums(N)
  z_earlyratiosummary[2, j] <- 1 - (s[j] / 401)
}

# Set rownames
rownames(z_earlyratiosummary) <- c("Value", "P-value")



#Stage 2

N<- data.frame(matrix(0L,nrow = 400, ncol = 15))

Inter_Ncp.2ratio<- rbind(Inter_ncp_currentratio[5,], Inter_ncp_futureratio) 

#Stage 2

# Initialize N matrix
N <- matrix(0, nrow = 401, ncol = 15)

# Loop through each column
for (j in 1:15) {
  # Loop through each row starting from the second row
  for (i in 2:401) {
    # Check if both values are not missing before comparison
    if (!is.na(Inter_Ncp.2ratio[i, j]) && !is.na(Inter_Ncp.2ratio[1, j]) && Inter_Ncp.2ratio[i, j] < Inter_Ncp.2ratio[1, j]) {
      N[i, j] <- 1
    }
  }
  # Calculate p-value for each column after updating N
  s <- colSums(N)
  z_interratiosummary[2, j] <- 1 - (s[j] / 401)
}

# Set rownames
rownames(z_interratiosummary) <- c("Value", "P-value")

#Stage 2

# Initialize N matrix
N <- matrix(0, nrow = 401, ncol = 15)

Late_Ncp.2ratio<- rbind(Late_ncp_currentratio[5,], Late_ncp_futureratio) 
# Loop through each column
for (j in 1:15) {
  # Loop through each row starting from the second row
  for (i in 2:401) {
    # Check if both values are not missing before comparison
    if (!is.na(Late_Ncp.2ratio[i, j]) && !is.na(Late_Ncp.2ratio[1, j]) && Late_Ncp.2ratio[i, j] < Late_Ncp.2ratio[1, j]) {
      N[i, j] <- 1
    }
  }
  # Calculate p-value for each column after updating N
  s <- colSums(N)
  zlateratiosummary[2, j] <- 1 - (s[j] / 401)
}

# Set rownames
rownames(zlateratiosummary) <- c("Value", "P-value")
zlateratiosummary[] <- lapply(zlateratiosummary, function(x) {
  x[is.nan(x)] <- 0
  return(x)
})
####Zscore plot 

###Early
z<- z_earlyratiosummary

dataz<- data.frame(
  NCP=colnames(z),
  Zscore=as.numeric(z[1,]),
  Pvalue=as.numeric(z[2,])
)  

P_zscore_earlyratio<-ggplot(dataz, aes(x=NCP, y=Zscore)) +
  geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Zscore")+ ggtitle("Zscore -Early stage Summary")+  scale_x_discrete(labels=labels)

P_zscore_earlyratio

##intermediate

z<- z_interratiosummary

dataz<- data.frame(
  NCP=colnames(z),
  Zscore=as.numeric(z[1,]),
  Pvalue=as.numeric(z[2,])
)  


P_zscore_interratio<-ggplot(dataz, aes(x=NCP, y=Zscore)) +
  geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Zscore")+ ggtitle("Zscore -Intermediate stage Summary")+  scale_x_discrete(labels=labels)

P_zscore_interratio


##Late

z<- zlateratiosummary

dataz<- data.frame(
  NCP=colnames(z),
  Zscore=as.numeric(z[1,]),
  Pvalue=as.numeric(z[2,])
)  


P_zscore_lateratio<-ggplot(dataz, aes(x=NCP, y=Zscore)) +
  geom_segment( aes(x=NCP, xend=NCP, y=0, yend=Zscore), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Zscore")+ ggtitle("Zscore -Late stage Summary")+  scale_x_discrete(labels=labels)

P_zscore_lateratio


figure_zsummary<- ggarrange(P_zscore_earlyratio,P_zscore_interratio,P_zscore_lateratio,
                            labels = letters[1:3],
                            ncol = 3, nrow = 1,
                            common.legend = TRUE, legend = "bottom")  

figure_zsummary 



labels<- c("Air quality", "Feed", "Food", "Habitat", "Medical Resources", "Natural Hazards", "Nutrient regulation", 
           "Ornamental Resources", "Pest and Disease", "Experiences", "Pollination", "Raw Materials","Learning and Inspiration", "Soil", "Supporting Identities")


write.csv(z_earlyratiosummary, "Zearlyratiosummary.csv", sep = ",")

write.csv(z_interratiosummary, "Zinterratiosummary.csv", sep = ",")

write.csv(zlateratiosummary, "Zlateatiosummary.csv", sep = ",")


write.csv(z_early, "z_early_tot.csv")
write.csv(z_inter, "z_inter_tot.csv")
write.csv(z_late, "z_late_tot.csv")



