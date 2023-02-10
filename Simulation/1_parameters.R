# Set working directory as the folder downloaded from 
# https://github.com/ZilongXie/LatentRegMissingKnockoff
setwd('') 

################################################################################
########################
## Correlation matrix ##  
########################
J <- 60
p <- 100
set.seed(123456789)
# Block-wise underlying correlation matrix
block1 <- 1:20
block2 <- 21:40
block3 <- 41:60
block4 <- 61:80
block5 <- 81:100
Sigma.true <- matrix(0, nrow = p, ncol = p)
Sigma.true[block1, block1] <- 0.6
Sigma.true[block2, block2] <- 0.6
Sigma.true[block3, block3] <- 0.6
Sigma.true[block4, block4] <- 0.3
Sigma.true[block5, block5] <- 0.3

Sigma.true[block1, block2] <- 0.15
Sigma.true[cbind(block1, block2)] <- 0.3
Sigma.true[block1, block3] <- 0.15
Sigma.true[block2, block3] <- 0.15

Sigma.true[block1, block4] <- 0.15
Sigma.true[cbind(block1, block4)] <- 0.3
Sigma.true[block2, block4] <- 0.15
Sigma.true[cbind(block2, block4)] <- 0.15
Sigma.true[block3, block4] <- 0.15

Sigma.true[c(block1, block2, block3, block4), block5] <- matrix(runif(80 * 20, min = 0.1, max = 0.2), nrow = 80)

Sigma.true[lower.tri(Sigma.true)] <- t(Sigma.true)[lower.tri(Sigma.true)]
diag(Sigma.true) <- 1


# Plot heatmap of the correlation matrix
require(ggplot2)

# Dummy data
x <- 0:99
y <- 99:0
corrplot <- expand.grid(X=x, Y=y)
corrplot$Correlation <- as.vector(Sigma.true)

# Heatmap 
ggplot(corrplot, aes(X, Y, fill= Correlation)) + 
  geom_tile() +
  scale_fill_gradient2(low="white", high="black", breaks = (0:5)/5) + 
  scale_x_continuous(expand = c(0,0), breaks= seq(10, 90, 20) - 0.5, limits = c(0, 99), labels = c('Block1','Block2','Block3','Block4','Block5')) +
  scale_y_continuous(expand = c(0,0), breaks= seq(10, 90, 20) - 0.5, limits = c(0, 99), labels = c('Block5','Block4','Block3','Block2','Block1')) +
  theme(panel.background = element_rect(fill = 'NA'),
        panel.ontop = TRUE,
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 


########################
## Unknown parameters ##
########################
bin.ind <- c(1:10, 21:30, 41:50, 61:70, 81:90)
con.ind <- c(11:20, 31:40, 51:60, 71:80, 91:100)
ord.ind <- NULL
con.params.true <- list()
con.params.true$mean <- rep(0, 50)
con.params.true$sd <- rep(1, 50)
bin.params.true <- as.vector(matrix(c(-1.2, -0.3, 0, 0.3, 1.2), nrow = 5, ncol = 10))
ord.params.true <- NULL
beta.true <- rep(0, 100)
# beta.true[c(1,2,21,22,41,51,61,71,81,91)] <- c(1, -1, 1, -1, 0.75, -0.75, 0.5, -0.5, 0.25, -0.25)/5

beta.true[c(1,11,22,32,43,53,64,74,85,95)] <- c(1, -1, 1, -1, 1, -1, 1, -1, 1, -1)/2

#################################
## Missing at random paramters ##
#################################
expit <- function(x) {1/(1+exp(-x))}
# MAR derived from means of the fifth-block
# c(-1, 1/20, 1/20, ..., 1/20)

####################
## IRT parameters ##
####################
a.vec <- runif(J, 0.5, 1)
d1.vec <- runif(J, -2, 0)
d2.vec <- NA * d1.vec

##############
## Booklets ##
##############
booklets <- matrix(0, nrow = 3, ncol = J)
booklets[1,1:20] <- 1
booklets[2,21:40] <- 1
booklets[3,41:60] <- 1
save( p, J, bin.ind, con.ind, ord.ind, con.params.true, bin.params.true, ord.params.true,
      beta.true, Sigma.true, block1, block2, block3, block4, block5, a.vec, d1.vec, d2.vec,
      booklets, file = './Simulation/simulation_params.RData')
