### PACKAGES ####
library(ggplot2)

#### LOAD IN DATA ####
#load in parameters
arrArgs <- commandArgs(trailingOnly = TRUE)
parameters <- as.character(arrArgs[1])

#load in csvs
rawGeno <- read.csv(paste0("../../outputs/",parameters,"/",parameters,"_RawGeno_All.csv"))
pheno <- read.csv(paste0("../../outputs/",parameters,"/",parameters,"_Phenotypes_All.csv"))
selModel <- rawGeno[1,8]
prefS <- as.numeric(levels(factor(rawGeno$prefS)))
pops <- levels(factor(rawGeno$pop,levels=unique(rawGeno$pop)))

#set up df for storing af
df_list <- vector("list",50000)

#### LOOP THROUGH AND CALCULATE AF + HETEROZYGOTE FREQUENCT####
rowIndex<-1
print("parsing genotypes and phenotypes")
for(n in prefS){
  for(i in 1:100){
    for(j in pops){
      for(k in 1:20){
        curRepPopGen <- rawGeno[rawGeno$prefS==n & 
                                  rawGeno$replicate==i & 
                                  rawGeno$pop==j & 
                                  rawGeno$gen==k,]
        
        curRepPopGenPheno <- pheno[pheno$prefS==n & 
                                       pheno$replicate==i & 
                                       pheno$pop==j & 
                                       pheno$gen==k,]
        
        #get het freq
        hetFreq <- (sum(curRepPopGen$locus1==14) + sum(curRepPopGen$locus1==41))/nrow(curRepPopGen)
        hetN <- sum(curRepPopGen$locus1==14) + sum(curRepPopGen$locus1==41)
        
        #calculate allele freq
        for(l in 1:1){
          curPopAlleles <- strsplit(as.character(curRepPopGen[,6+l]),"")
          nInteracting <- 0
          for(m in 1:length(curPopAlleles)){
            if(curPopAlleles[[m]][1] == "1" &&
               curPopAlleles[[m]][2] == "1"){
              nInteracting <- nInteracting + 2
            } else if (curPopAlleles[[m]][1] == "1" ||
                       curPopAlleles[[m]][2] == "1"){
              nInteracting <- nInteracting + 1
            }
          }
          interactAF <- nInteracting/(length(curPopAlleles)*2)
          interactN <- nInteracting
        }
        df_list[[rowIndex]] <- c(i,j,k,n,
                                 interactAF,interactN,
                                 hetFreq,hetN,
                                 mean(curRepPopGenPheno$perHet),
                                 mean(curRepPopGenPheno$perMal))
        rowIndex <- rowIndex+1
      }
    }
  }
}
print("finished parsing")

df <- as.data.frame(do.call(rbind,df_list))
colnames(df) <- c("replicate","pop","gen","prefS",
                  "locus1","nBirLocus1",
                  "hetFreq","hetN","perHet","perMal")

#### PLOT ####

#density plot theme
theme_dens <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background= element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.text.y = element_text(size = 10),
                    axis.title = element_text(face = "bold",
                                              size = 10),
                    legend.title = element_text(face = "bold",
                                                size = 10,
                                                hjust=0.5),
                    legend.key=element_blank(),
                    plot.title = element_text(face = "bold",
                                              size = 17,
                                              hjust=0.5))

#density plot for af at incomp locus (10 generations)
locusAF <- ggplot(data=subset(df,gen=="10"),aes(x=as.numeric(locus1),
                                                               fill=as.factor(prefS)))+
  geom_density(alpha=0.5)+
  scale_fill_viridis_d()+
  #xlim(0,0.6)+
  labs(fill="Strength of preference")+
  xlab("Birchmanni allele frequency")+
  ylab("Density")+
  theme_dens

#density plot for genome wide heterozygosity (10 generations)
perHet <- ggplot(data=subset(df,gen=="10"),aes(x=as.numeric(perHet)/100,
                                                               fill=as.factor(prefS)))+
  geom_density(alpha=0.5)+
  scale_fill_viridis_d()+
  #xlim(0.4,0.6)+
  labs(fill="Strength of preference")+
  xlab("Genome-wide heterozygosity")+
  ylab("Density")+
  theme_dens

#density plot for genome wide mal ancestry (10 generations)
perMal <- ggplot(data=subset(df,gen=="10"),aes(x=as.numeric(perMal)/100,
                                                     fill=as.factor(prefS)))+
  geom_density(alpha=0.5)+
  scale_fill_viridis_d()+
  #xlim(0,0.6)+
  labs(fill="Strength of preference")+
  xlab("Genome-wide malinche ancestry")+
  ylab("Density")+
  theme_dens

#### SAVE ####
ggsave(locusAF,
       filename = paste0(selModel,"_locusAF.pdf"),
       width = 5.31,
       height = 4.60,
       units = "in")

ggsave(perHet,
       filename = paste0(selModel,"_perHet.pdf"),
       width = 5.31,
       height = 4.60,
       units = "in")

ggsave(perMal,
       filename = paste0(selModel,"_perMal.pdf"),
       width = 5.31,
       height = 4.60,
       units = "in")




