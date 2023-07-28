# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

##############################
# Extract data from sas file #
##############################
library(haven)
data <- read_sas('./PISA_2015/cy6_ms_cmb_stu_cog.sas7bdat')
write.csv(data[which(data[,'CNT'] == 'USA'),],
          file = './PISA_2015/2015_USA_COG.csv', row.names = F)

#################
# Extract items #
#################
cog <- read.csv('./PISA_2015/2015_USA_COG.csv')

# unique(substr(names(cog)[which(colSums(is.na(cog)) < dim(cog)[1])], 1, 2))
CS <- grep(names(cog)[which(colSums(is.na(cog)) < dim(cog)[1])], pattern = 'CS', value = T)
DS <- grep(names(cog)[which(colSums(is.na(cog)) < dim(cog)[1])], pattern = 'DS', value = T)

CS.scored <- grep(CS, pattern = 'S$', value = T) # match the end
DS.scored <- c(grep(DS, pattern = '[0-9]C$', value = T),
               grep(DS, pattern = '[a-z]C$', value = T)) # match the end

cog <- cog[,c('CNTSTUID',CS.scored, DS.scored)]
cog.max <- apply(cog, 2, max, na.rm=T)

# recode some items
dic.recode <- cog[,which((cog.max > 2) & (cog.max <= 12))]
poly.recode <- cog[,which((cog.max > 12))[-1]]
dic.recode[(dic.recode >= 1) & (dic.recode < 11)] <- 0
dic.recode[(dic.recode >= 10) & (dic.recode <= 12)] <- 1
poly.recode[(poly.recode >= 1) & (poly.recode < 11)] <- 0
poly.recode[(poly.recode >= 11) & (poly.recode <= 12)] <- 1
poly.recode[(poly.recode >= 21)] <- 2
poly.recode <- as.data.frame(poly.recode)
cog[,which((cog.max > 2) & (cog.max <= 12))] <- dic.recode
cog[,which((cog.max > 12))[-1]] <- poly.recode


########################
# Save item parameters #
########################
# Table-Extraction_params.csv and Table-Extraction_ref.csv are extracted from
# Annex A of 'PISA 2015 Technical Report'.
ref <- read.csv('./PISA_2015/Table-Extraction_ref.csv', header = F)
names(ref) <- c('Generic ID', 'CBA item ID')
ref$`Generic ID`[which(ref$`CBA item ID` == 'DS304Q03aC')] = 'S304Q03a'
ref$`Generic ID`[which(ref$`CBA item ID` == 'DS304Q03bC')] = 'S304Q03b'
params <- read.csv('./PISA_2015/Table-Extraction_params.csv', header = F)
names(params) <- c('Generic ID', 'Slope','Difficulty','Step1','Step2')
params$`Generic ID`[which(params$`Generic ID` == 'S304Q03')] <- c('S304Q03a', 'S304Q03b')
ref_params <- merge(ref, params, by = 'Generic ID')[,-1]
ref_params$Step1[is.na(ref_params$Step1)] <- 0

a.vec <- c() 
for (item in names(cog)[-1]) {
  if (item %in% ref_params$`CBA item ID`){
    a.vec <- c(a.vec, ref_params$Slope[which(ref_params$`CBA item ID` == item)])
  }
  else{
    print(item)
  }
}

# Note that 'step 1' + 'step 2' = 0
d1.vec <- c()
b.vec <- c()
for (item in names(cog)[-1]) {
  if (item %in% ref_params$`CBA item ID`){
    a <- ref_params$Slope[which(ref_params$`CBA item ID` == item)]
    b <- ref_params$Difficulty[which(ref_params$`CBA item ID` == item)]
    b.vec <- c(b.vec , b)
    d1 <- ref_params$Step1[which(ref_params$`CBA item ID` == item)]
    d1.vec <- c(d1.vec, a*(d1 - b))
  }
  else{
    print(item)
  }
}

d2.vec <- c()
for (item in names(cog)[-1]) {
  if (item %in% ref_params$`CBA item ID`){
    a <- ref_params$Slope[which(ref_params$`CBA item ID` == item)]
    b <- ref_params$Difficulty[which(ref_params$`CBA item ID` == item)]
    d2 <- ref_params$Step2[which(ref_params$`CBA item ID` == item)]
    d2.vec <- c(d2.vec, a*(d2 - b))
  }
  else{
    print(item)
  }
}
################################################################################
save(cog, a.vec, d1.vec, d2.vec, file = '.Real_data/2015_USA_COG.RData')
