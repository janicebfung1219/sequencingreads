##remove all 0 expressions and NA

file="mutans_table.csv"

#echo ${file}
cat expout/${file} | sed '/NA,NA,NA,NA,NA,NA,NA/d'|sed '/0,0,0,0,0,0,0,0/d' > expout/clean_${file}
wc -l expout/clean_${file}
