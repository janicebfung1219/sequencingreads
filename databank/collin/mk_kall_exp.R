library(dplyr)
#setwd("databank/collin")
###col names are ""target_id"  "length"     "eff_length" "est_counts" "tpm" "

##pattern + group
groups <- c("nissle", "mutans")
pattern <- c("_N", "_M")
schema <- data.frame(groups, pattern)
print(schema)
  
##function to extract out each fastq file and grab expression(tpm)
sample_exp <- function(file, sample_name){
  file <- paste0("kall_out/",sample_name,"/abundance.tsv")
  df <- read.delim(file)
  return(df[,c(1,5)])
}


##collect all files in directory "kall_out"
sample_name <- c(list.dirs(path = "kall_out/", full.names = FALSE, recursive = FALSE))
print(sample_name)


##extract tpm from each file in directory 
exp_table <- data.frame()
for (name in sample_name){
  print(name)
  file <- paste0("kall_out/",name,"/abundance.tsv")
  #print(file)

  df_tmp <- sample_exp(file, name)
  colnames(df_tmp)[2] <- name;

  if( dim(exp_table)[1]==0){ exp_table = df_tmp}
  else { exp_table <- full_join(exp_table,df_tmp,by = "target_id")}
  #print(colnames(exp_table))

}

for (i in 1:nrow(schema)){
  name <- paste0(schema[i,"groups"], "_table")
  pattern <- schema[i,"pattern"]
  table <- cbind(exp_table$target_id, exp_table %>% select(matches(pattern)))
  filename <- paste0("expout/", schema[i,"groups"], "_table.csv")
  
  write.csv(table, file=filename, row.names = FALSE)
  print(filename)
}
