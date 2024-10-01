#store command line arguements as variables
arrArgs <- commandArgs(trailingOnly = TRUE)
sexsel_file<-as.character(arrArgs[1])
disCoeff <- as.numeric(arrArgs[3])/10
if(disCoeff <= 0.6){
  increment_value<-as.numeric(arrArgs[2])/10
} else {
  increment_value<-as.numeric(arrArgs[2])/10+0.01
}

my_func = function(increment_value,disCoeff){
  sexsel = read.delim(paste(sexsel_file), sep = "\t")
  sexsel$Selection[3:5] = paste0("s=",disCoeff+increment_value,
                                 ",1-((1-(Courter_perHet/100))*s)")
  write.table(sexsel, paste(sexsel_file), sep = '\t', quote = FALSE, row.names = FALSE)
}

my_func(increment_value,disCoeff)
