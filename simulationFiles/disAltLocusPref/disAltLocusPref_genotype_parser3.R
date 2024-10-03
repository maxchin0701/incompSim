#### PARSING THROUGH GENOTYPES ####
#### PACKAGES ####
library(hierfstat)
#store command line arguements as variables
arrArgs <- commandArgs(trailingOnly = TRUE)
parameters<-as.character(arrArgs[1])
prefCoeff <- as.numeric(arrArgs[2])/10
if(prefCoeff == 0.9){
  prefCoeff <- prefCoeff+0.01
}

# parameters_list = strsplit(parameters,"_")
# project<- parameters_list[[1]][1]
# esize<- parameters_list[[1]][2]
# ss<- parameters_list[[1]][3]
# ns<- parameters_list[[1]][4]

#### SETTING UP DATA STRUCTURES AND INDICIES ####
#create a dataframe to store the mean/sd/cv phenotype for each iteration/replicate
df <- data.frame(NULL)
df1 <- data.frame(NULL)
df2 <- data.frame(NULL)
df3 <- data.frame(NULL)

#Create temporary lists to store
df_list <- vector("list",5000)
df1_list <- vector("list",5000)
df2_list <- vector("list",5000)
df3_list <- vector("list",5000)

#Create index variables
index <- 1
index1 <- 1

#### LOOP THROUGH ITERATIONS ####
#loop through directories/iterations
for(i in 1:100){
  setwd(paste0("./",parameters,"_",i))
  # f = 'hybrids101'
  # if(file.exists(f)==FALSE){
  #   num_failed_sims = num_failed_sims + 1
  #   next
  # }
  
  #Create list for genes files
  hg <- list()
  
  #Create list for phenotypes files
  hp <- list()
  
  #Create vector for male dataframe names
  hv <- c()
  
  #### LOAD GENOTYPES FOR GENERATIONS/POPS ####
  for(k in 1:20){
    ##Create dataframe names
    pop1_geno <- paste0("h",k,"l_geno",sep = "")
    pop2_geno <- paste0("h",k,"m_geno",sep = "")
    pop3_geno <- paste0("h",k,"h_geno",sep = "")
    pop1_pheno <- paste0("h",k,"l_pheno",sep = "")
    pop2_pheno <- paste0("h",k,"m_pheno",sep = "")
    pop3_pheno <- paste0("h",k,"h_pheno",sep = "")
    
    #load in dataframe
    curGeno <- read.delim(file = paste0("Gen",k,"_genes.txt"), header = F)
    curPheno <- read.delim(file = paste0("Gen",k,"_phenotypes.txt"), header = T)
    
    ##Assign dataframes to names (add male data frames to list)
    hg[[length(hg) + 1 ]] <- assign(pop1_geno,curGeno[which(curGeno[,2] == 3),])
    hg[[length(hg) + 1 ]] <- assign(pop2_geno,curGeno[which(curGeno[,2] == 4),])
    hg[[length(hg) + 1 ]] <- assign(pop3_geno,curGeno[which(curGeno[,2] == 5),])
    hp[[length(hp) + 1 ]] <- assign(pop1_pheno,curPheno[which(curPheno[,2] == 3),])
    hp[[length(hp) + 1 ]] <- assign(pop2_pheno,curPheno[which(curPheno[,2] == 4),])
    hp[[length(hp) + 1 ]] <- assign(pop3_pheno,curPheno[which(curPheno[,2] == 5),])

    ##Add dataframe names to vector
    hv <- c(hv,pop1_geno,pop2_geno,pop3_geno)
  }
  
  #### LOOP THROUGH GENERATION/POPULATION DFS ####
  for(j in 1:length(hg)){
    ##Get generation and site from data name
    if(nchar(hv[j]) == 9){
      gen <- substr(hv[j],2,3)
      
      ##Get site from data name
      if(substr(hv[j],4,4) == "l"){
        site <- "STL"
      } else if (substr(hv[j],4,4) == "m"){
        site <- "STM"
      } else if (substr(hv[j],4,4) == "h"){
        site <- "STH"
      }
    } else if(nchar(hv[j]) == 8){
      gen <- substr(hv[j],2,2)
      
      ##Get site from data name
      if(substr(hv[j],3,3) == "l"){
        site <- "STL"
      } else if (substr(hv[j],3,3) == "m"){
        site <- "STM"
      } else if (substr(hv[j],3,3) == "h"){
        site <- "STH"
      }
    }
    
    #### SUMMARIZE GENOTYPES FOR INDIVIDUALS ####

    #get individual IDs
    indID <- levels(factor(hg[[j]][,1],levels=unique(hg[[j]][,1])))
    
    #set up list for individuals
    indGenoList <- vector("list",1000)
    indSex <- vector("list",100)
    
    #loop through individuals
    for(l in indID){
      #get corresponding rows
      curIndHaps <- hg[[j]][which(hg[[j]][,1] == l),]
      
      #set locus counter and column
      locusN <- 1
      locusColumn <- 23
      curIndGeno <- vector("list",100)
      curIndGeno[[1]] <- curIndHaps[1,2]
      
      #loop through loci
      while(locusN <= 1){
        #skip if SDL
        if(curIndHaps[1,locusColumn-1] == "Sm" || 
           curIndHaps[1,locusColumn-1] == "sb" ||
           curIndHaps[1,locusColumn-1] == "sm"){
          locusColumn <- locusColumn + 2
          next
        }
        #append genotype to list
        curIndGeno[[locusN+1]] <- paste0(curIndHaps[1,locusColumn],
                                         curIndHaps[2,locusColumn])
        
        #increment locusN and locusColumn
        locusN <- locusN + 1
        locusColumn <- locusColumn + 3
      }
      
      #bind and append to ind list
      indGenoList[[which(indID==l)]] <- as.numeric(do.call(cbind,curIndGeno))
      
      #get sex
      SDLGeno <- paste0(curIndHaps[1,70],
                        curIndHaps[2,70])
      if(SDLGeno=="10"){
        indSex[[which(indID==l)]] <- 1
      } else {
        indSex[[which(indID==l)]] <- 0
      }
    }
    
    #bind list 
    indGenoDf <- as.data.frame(do.call(rbind,indGenoList))
    
    #set colnames
    colnames(indGenoDf) <- c("pop",paste0("locus",1))

    #### CALCULATE FIS AND H0 ####
    curPopFis <- basic.stats(data=indGenoDf,diploid=T)$Fis[1,1]
    curPopHo <- basic.stats(data=indGenoDf,diploid=T)$Ho[1,1]
    
    #### APPEND TO LIST ####
    #append to df
    df_list[[index]] <- cbind(i,
                              indID,
                              do.call(rbind,indSex),
                              site,
                              gen,
                              prefCoeff,
                              indGenoDf[,-1])
    
    #append to df1
    df1_list[[index1]] <- c(i,
                               site,
                               gen,
                               prefCoeff,
                               curPopFis)
    
    #append to df3
    df2_list[[index1]] <- c(i,
                            site,
                            gen,
                            prefCoeff,
                            curPopHo)
    
    #append to df3
    df3_list[[index]] <- cbind(i,
                              indID,
                              do.call(rbind,indSex),
                              site,
                              gen,
                              prefCoeff,
                              hp[[j]]$perHet,
                              hp[[j]]$perMal)
    
    #Increment indices
    index <- index + 1
    
    index1 <- index1 + 1
    
  }
  setwd("..")
}

#### PREPARE FINAL DF FOR SAVING ####
#Bind all elements in lists
df <- rbind(df,do.call(rbind,df_list))

df1 <- rbind(df1,do.call(rbind,df1_list))

df2 <- rbind(df2,do.call(rbind,df2_list))

df3 <- rbind(df3,do.call(rbind,df3_list))

#Set column names
colnames(df)[1:7] = c("replicate","ID","sex","pop","gen","prefS","locus1")

colnames(df1)[1:5] = c("replicate","pop","gen","prefS","Fis")

colnames(df2)[1:5] = c("replicate","pop","gen","prefS","H0")

colnames(df3)[1:8] = c("replicate","ID","sex","pop","gen","prefS","perHet","perMal")

#add all parameters to df dataframe for writing
df$selModel <- parameters

#add all parameters to df1 dataframe for writing
df1$selModel <- parameters

#add all parameters to df1 dataframe for writing
df2$selModel <- parameters

df3$selModel <- parameters

#### SAVE ####
write.csv(df, paste0("../../outputs/",parameters,"/",parameters,"_RawGeno_",prefCoeff,".csv"), quote = FALSE, row.names = FALSE)
write.csv(df1, paste0("../../outputs/",parameters,"/",parameters,"_Fis_",prefCoeff,".csv"), quote = FALSE, row.names = FALSE)
write.csv(df2, paste0("../../outputs/",parameters,"/",parameters,"_Ho_",prefCoeff,".csv"), quote = FALSE, row.names = FALSE)
write.csv(df3, paste0("../../outputs/",parameters,"/",parameters,"_Phenotypes_",prefCoeff,".csv"), quote = FALSE, row.names = FALSE)
