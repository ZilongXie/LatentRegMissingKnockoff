# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

source('./R_code/gci_estimate.R')
source('./R_code/knockoff_construct.R')
source('./R_code/knockoff_select.R')
source('./R_code/stem.R')

X <- read.csv('./Real_data/USA_X.csv')
Y <- read.csv('./Real_data/USA_Y.csv')
load('./Real_data/Item_parameters_2015.RData')

bin.ind <- 1:20
ord.ind <- 21:33
con.ind <- 34:62

set.seed(123456789, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind='default')
estimate.res <- GCI.estimate(X = X, 
                             bin.ind = bin.ind,
                             ord.ind = ord.ind,
                             con.ind = con.ind,
                             max.iter = 10000,
                             burn.in = 5000,
                             print.iter = T,
                             is.params.trace = T,
                             is.return.under = T,
)
con.params.es <- list()
con.params.es$mean <- estimate.res$params.avg$con.params.mean
con.params.es$sd <- estimate.res$params.avg$con.params.sd
bin.params.es <- estimate.res$params.avg$bin.params
ord.params.es <- estimate.res$params.avg$ord.params
Sigma.es <- estimate.res$params.avg$Sigma
Omega.es <- solve(Sigma.es)
X.under0 <- estimate.res$X.under
S <- MVR.knockoffs(Sigma.es)
print(min(diag(S)))

latreg.es <- StEM(X=X, X.under=X.under0, X.knock=NULL, X.knock.under=NULL, Y=Y, 
                  Sigma=Sigma.es, S=NULL, a.vec=a.vec, d1.vec=d1.vec, d2.vec=d2.vec, 
                  bin.ind=bin.ind, ord.ind=ord.ind, con.ind=con.ind, 
                  bin.params=bin.params.es, ord.params=ord.params.es, con.params=con.params.es, 
                  max_iter=2000, burn_in=1000, is.beta.trace=T, is.print.iter=T, 
                  is.knockoff=F, lambda=0) # Linear regression

beta.es <- latreg.es$beta
X.under <- latreg.es$X.under
theta.interp.es <- latreg.es$interp
theta.var.es <- latreg.es$var
XX <- latreg.es$XX
XX.ind <- latreg.es$XX.ind

save(bin.ind, ord.ind, con.ind, bin.params.es, ord.params.es, con.params.es, 
     X, X.under, XX, XX.ind, Sigma.es, Omega.es, S, beta.es,  theta.interp.es, theta.var.es, latreg.es,
     file = './Real_data/Real_data_estimate.RData')
