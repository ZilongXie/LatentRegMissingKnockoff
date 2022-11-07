load('./Real_data/2015_USA_STU.RData')
load('./Real_data/2015_USA_COG.RData')

##############################################
# Merge questionaire data amd cognitive data #
##############################################
ALL <- merge(data, cog, by = 'CNTSTUID')
Y.names <- names(cog)[-1]
X.names <- names(data)[-1]
X <- ALL[,X.names]
Y <- ALL[,Y.names]
save(a.vec, d1.vec, d2.vec, file = './Real_data/Item_parameters_2015.RData')

##########################
# Common functions       #
##########################
rowmax <- function(x, y){
  temp <- rep(0, length(x))
  z <- cbind(x, y)
  temp[rowSums(is.na(z)) < 2] <- apply(z[which(rowSums(is.na(z)) < 2),], 1, max, na.rm = T)
  temp[rowSums(is.na(z)) == 2] <- NA
  temp
}

count.factors <- function(x){
  y <- as.factor(x)
  length(levels(y))
}

poly.trans <- function(x){
  as.numeric(as.factor(x))
}

min_prop <- function(x){t <- table(x);min(t)/sum(t)}
merge_level <- function(x, t = 0.05){
  tx <- table(x)
  v <- sort(unique(x))
  if(length(tx) == 2){
    return(x)
  }
  minlevel <- min(tx)
  sumlevel <- sum(tx)
  prop <- tx/sumlevel
  if(min(prop) > t){
    return(x)
  }
  loc <- which(tx == minlevel)
  if(loc == 1){
    x[which(x == v[1])] <- v[2]
    print(paste0(1, ' to ', 2))
  }
  if(loc == length(tx)){
    x[which(x >= v[length(tx) - 1])] <- v[length(tx) - 1]
    print(paste0(loc, ' to ', loc - 1))
  }
  if((loc > 1) & (loc < length(tx))){
    if(prop[loc-1] <  prop[loc+1]){
      x[which(x == loc)] <- v[loc - 1]
      print(paste0(loc, ' to ', loc - 1))
    }
    if(prop[loc-1] >  prop[loc+1]){
      x[which(x == loc)] <- v[loc + 1]
      print(paste0(loc, ' to ', loc + 1))
    }
  }
  return(poly.trans(x))
}

###################
# select ST items #
###################
# Select items that are not used to construct index
ST <- grep(X.names, pattern = 'ST[0-9]', value = T)
ST.items <- c('ST004D01T', # Gender
              # 'ST019AQ01T', # Country of Birth International-Self, correlated with ST022Q01TA
              # 'ST019BQ01T', # Country of Birth International-Mother, correlated with ST022Q01TA
              # 'ST019CQ01T', # Country of Birth International-Father, correlated with ST022Q01TA
              # 'ST021Q01TA', # Age of Arrived
              'ST022Q01TA', # International Language at Home
              'ST111Q01TA', # Level that Expect to Complete 
              # 'ST121Q01NA', 
              # 'ST121Q02NA',
              # 'ST121Q03NA',
              'ST062Q01TA', # Skip school day
              'ST062Q02TA', # Skip classes
              'ST062Q03TA', # Arrive late
              'ST031Q01NA', # Moderate Exercises Days per week
              'ST032Q01NA', # Hard Exercises Days per week
              # 'ST063Q01NA', # Physics
              # 'ST063Q01NB', # Physics
              # 'ST063Q02NA', # Chemistry
              # 'ST063Q02NB', # Chemistry
              # 'ST063Q03NA', # Biology
              # 'ST063Q03NB', # Biology
              # 'ST063Q04NA', # Earth
              # 'ST063Q04NB', # Earth
              # 'ST063Q05NA', # Applied
              # 'ST063Q05NB', # Applied
              # 'ST063Q06NA', # General
              # 'ST063Q06NB', # General
              'ST064Q01NA', # Choose Courses
              'ST064Q02NA', # Choose Difficulty
              'ST064Q03NA' # Choose Number
              # 'ST071Q01NA', # Out of School Time: Science
              # 'ST071Q02NA', # Out of School Time: Mathematics
              # 'ST071Q03NA', # Out of School Time: Test Language
              # 'ST071Q04NA', # Out of School Time: Foreign Language
              # 'ST071Q05NA', # Out of School Time: Other
              # 'ST076Q01NA', # Before School: Breakfast
              # 'ST076Q02NA', # Before School: Study
              # 'ST076Q03NA', # Before School: TV
              # 'ST076Q04NA', # Before School: Book
              # 'ST076Q05NA', # Before School: Social
              # 'ST076Q06NA', # Before School: Game
              # 'ST076Q07NA', # Before School: Friend
              # 'ST076Q08NA', # Before School: Parents
              # 'ST076Q09NA', # Before School: Household
              # 'ST076Q10NA', # Before School: Work
              # 'ST076Q11NA', # Before School: Exercise
              # 'ST078Q01NA', # After School: Dinner
              # 'ST078Q02NA', # After School: Study
              # 'ST078Q03NA', # After School: TV
              # 'ST078Q04NA', # After School: Book
              # 'ST078Q05NA', # After School: Social
              # 'ST078Q06NA', # After School: Game
              # 'ST078Q07NA', # After School: Friend
              # 'ST078Q08NA', # After School: Parents
              # 'ST078Q09NA', # After School: Household
              # 'ST078Q10NA', # After School: Work
              # 'ST078Q11NA' # After School: Exercise
              # 'ST065Class' # Student Coded Science Class
)


# Combine ST063: Science classes
courses <- NULL
for(i in 1:6){
  courses <- cbind(courses, rowmax(X[,paste0('ST063Q0',i,'NA')],
                                   X[,paste0('ST063Q0',i,'NB')]))
}
courses <- as.data.frame(courses)
names(courses) <- c('SCI.PHY','SCI.CHE','SCI.BIO',
                    'SCI.EAR','SCI.APP','SCI.GEN')

# Combine ST076 and ST078: Out of school activities
activities <- NULL
for(i in 1:11){
  if(i <= 9){
    activities <- cbind(activities, rowmax(X[,paste0('ST076Q0',i,'NA')],
                                           X[,paste0('ST078Q0',i,'NA')]))  
  }
  if(i >= 10){
    activities <- cbind(activities, rowmax(X[,paste0('ST076Q',i,'NA')],
                                           X[,paste0('ST078Q',i,'NA')]))  
  }
}
activities <- as.data.frame(activities)
names(activities) <- c('OUT.MEAL','OUT.STUDY','OUT.VEDIO',
                       'OUT.READ','OUT.NET','Out.GAME',
                       'OUT.FRD', 'OUT.PAR','OUT.HOUSE',
                       'OUT.JOB', 'OUT.SPORT')

######################
##  select indexes  ##
######################
indexes <- X.names[which(!(X.names %in% ST))]
indexes <- c(
  # 'AGE',
  'GRADE',
  # 'ISCEDL', # Correlated with REPEAT
  # 'ISCEDD', # All Same
  # 'ISCEDO' # All Same
  'DISCLISCI', 
  'TEACHSUP',
  'IBTEACH',
  'TDTEACH',
  'ENVAWARE',
  'ENVOPT',
  'JOYSCIE',
  'INTBRSCI',
  'INSTSCIE',
  'SCIEEFF',
  'EPIST',
  'SCIEACT',
  'BSMJ',
  # 'IMMIG', # ST019
  'MISCED',
  'FISCED',
  # 'HISCED',
  'BMMJ1',
  'BFMJ2',
  # 'hisei', 
  'REPEAT',
  'DURECEC', 
  'OUTHOURS',
  # 'MMINS', # TMINS
  # 'LMINS', # TMINS
  # 'SMINS', # TMINS
  'TMINS',
  'BELONG',
  'ANXTEST',
  'MOTIVAT',
  'COOPERATE',
  'CPSVALUE',
  'EMOSUPS',
  'PERFEED',
  'ADINST',
  'unfairteacher', 
  # 'PARED',
  # 'LANGN',
  'CULTPOSS',
  'HEDRES',
  # 'HOMEPOS',
  # 'ICTRES',
  'WEALTH'
  # 'ESCS'
)

X <- cbind(X[,c(ST.items, indexes)], courses, activities)

################
# Drop NA rows #
################
dropped <- which((rowSums(is.na(X)) <= 60) & (rowSums(is.na(Y)) <= 180))
Y <- Y[dropped,]
X <- X[dropped,]

######################
# Log-transformation #
######################
X[which(X[,'OUTHOURS'] == 0), 'OUTHOURS'] = NA
X[,'OUTHOURS'] <- log(X[,'OUTHOURS'])

#################
# Creat indices #
#################
num.facs <- (apply(X, 2, count.factors))
bin.ind <- which(num.facs == 2)
ord.ind <- which(num.facs >= 3 & num.facs <= 10)
con.ind <- which(num.facs >= 20)

####################
# Merge ord-levels #
####################
while(1){
  print(which(apply(X[,ord.ind], 2, min_prop) < 0.05))
  X[,ord.ind] <- apply(X[,ord.ind], 2, merge_level)
  if(min(apply(X[,ord.ind], 2, min_prop)) >= 0.05){
    break()
  }
}

###############
# New indices #
###############
num.facs <- (apply(X, 2, count.factors))
con.ind <- which(num.facs >= 9)
ord.ind <- which(num.facs >= 3 & num.facs <= 8)
bin.ind <- which(num.facs == 2)

# Reorder
X <- X[,c(bin.ind, ord.ind, con.ind)]
num.facs <- (apply(X, 2, count.factors))
con.ind <- which(num.facs >= 9)
ord.ind <- which(num.facs >= 3 & num.facs <= 8)
bin.ind <- which(num.facs == 2)

###################
# Discrete coding #
###################
X[,c(bin.ind, ord.ind)] <- apply(X[,c(bin.ind, ord.ind)], 2, poly.trans) - 1


#####################
# Check correlation #
#####################
require(psych)
res <- psych::mixedCor(as.data.frame(X),
                       c = con.ind,
                       d = bin.ind,
                       p = ord.ind,
                       correct = T,
                       smooth = T)
rho = res$rho
rho[upper.tri(rho, diag = T)] = 0
rho[which(abs(rho)>0.7, arr.ind = T)]

#############
# Save data #
#############
write.csv(X, file = './Real_data/USA_X.csv', row.names = F)
write.csv(Y, file = './Real_data/USA_Y.csv', row.names = F)

#############################
# Histogram of missing data #
#############################
# Plot missing proportion
library(ggplot2)
missing.data <- data.frame(Missing = colMeans(is.na(X)))
ggplot(data = missing.data, aes(x = Missing)) + 
  geom_histogram(binwidth = 0.05, color="black", fill='gray') + 
  xlab("Missing propotion") + 
  ylab("Count")  + 
  theme(panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
max(colMeans(is.na(X)))

