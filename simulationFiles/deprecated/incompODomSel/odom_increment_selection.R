#store command line arguements as variables
arrArgs <- commandArgs(trailingOnly = TRUE)
natsel_file<-as.character(arrArgs[1])
odomCoeff <- as.numeric(arrArgs[3])/10
if(odomCoeff <= 0.6){
  increment_value<-as.numeric(arrArgs[2])/10
} else {
  increment_value<-as.numeric(arrArgs[2])/10+0.01
}

my_func = function(increment_value,odomCoeff){
  natsel = read.delim(paste(natsel_file), sep = "\t")
  natsel$Selection[3:5] = paste0("s1=0.91,s2=",odomCoeff+increment_value,
                                 ",if(incomp==2,1-s1,if(incomp==8,1-s2,1))")
  write.table(natsel, paste(natsel_file), sep = '\t', quote = FALSE, row.names = FALSE)
}

my_func(increment_value,odomCoeff)
