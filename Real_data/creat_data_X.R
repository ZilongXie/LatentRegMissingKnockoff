# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

##############################
# Extract data from sas file #
##############################
library(haven)
library(tidyverse)
df = read_sas('./PISA_2015/cy6_ms_cmb_stu_qqq.sas7bdat')
write.csv(df[which(df$CNT == 'USA'),],
          file = './PISA_2015/2015_USA_STU.csv',
          row.names = F)

#################
# Extract items #
################# 
data <- read.csv('./PISA_2015/2015_USA_STU.csv')

# Discard country and school id
data <- data[,-which(names(data) %in% c('CNTRYID','CNT','CNTSCHID'))]
# Discard plausible values
data <- data[,-which(names(data) %in% grep(pattern = 'PV', x = names(data), value = T))]
# Discard weights
data <- data[,-which(names(data) %in% grep(pattern = 'W_FST', x = names(data), value = T))]
data <- data[,-which(names(data) %in% c('WVARSTRR','SENWT'))]
# Discard coding
data <- data[,-which(names(data) %in% c('SUBNATIO','STRATUM','OECD','NatCen','CYC','PROGN'))]
data <- data[,-which(names(data) %in% c('COBN_M','COBN_F','COBN_S'))]
data <- data[,-which(names(data) %in% c('OCOD1','OCOD2','OCOD3'))]
# Discard RANDOM NUMBER
data <- data[,-which(names(data) %in% c('UNIT'))]
data <- data[,-which(names(data) %in% c('CBASCI'))]
# Discard date
data <- data[,-which(names(data) %in% c('VER_DAT'))]
# Discard response mode(Paper/Computer)
data <- data[,-which(names(data) %in% c('ADMINMODE'))]
# Discard option
data <- data[,-which(names(data) %in% grep(names(data), pattern = 'Option*', value = T))]
# Discard LANGTEST
data <- data[,-which(names(data) %in% grep(names(data), pattern = 'LANGTEST*', value = T))]
# Discard region
data <- data[,-which(names(data) %in% c('Region'))] 
# Discard BOOKID
data <- data[,-which(names(data) %in% c('BOOKID'))]
# Discard variables with only missing values
data <- data[,which(colSums(is.na(data)) < dim(data)[1])]

save(data, file = '.Real_data/2015_USA_STU.RData')