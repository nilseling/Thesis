# Load libraries
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# Read in Dictyostelium data
input <- read.table("/Users/nils/Google Drive File Stream/My Drive/BASiCS_add-on/Data/Test_Data/Dictyostelium.txt")

### Test fit on one condition first
Day0 <- input[,grepl("X0h", colnames(input))]
Day0 <- Day0[rowMeans(Day0) > 1,]

Means <- rowMeans(Day0)
CV2 <- apply(Day0, 1, function(n){(sd(n)/mean(n))^2})

plot(log(Means), log(CV2), pch = 16, xlab = "mu", ylab = "delta")

### Do this by hand
# Create Basis functions - RGB - knots
x = log(Means) # the timestamps
n = 10 # the number of basis functions over T

# Fixed m and h
range = diff(range(x))
myu = seq(min(x), max(x), length.out = n)
h = diff(myu)*1.2

B <- matrix(NA,length(x),n+2)

# Intercept
B[,1] <- 1

# Linear component
B[,2] <- x

for (j in 3:(n+2)){
  B[,j] = exp(-0.5*(x-myu[j-2])^2/(h[1]^2))
}

B.df <- data.frame(Means = log(Means), B, CV2 = log(CV2))

ggplot(data = B.df) + 
  geom_point(aes(Means, CV2)) +
  geom_line(aes(x = Means, y = X1), colour = brewer.pal(name = "Set2", n = 7)[1], size = 1) +
  geom_line(aes(x = Means, y = X2), colour = brewer.pal(name = "Set2", n = 7)[2], size = 1) +
  geom_line(aes(x = Means, y = X3), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X4), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X5), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X6), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X7), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X8), colour =brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X9), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X10), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X11), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X12), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) 


### Least square approach y(hat) = X(X'X)^-1X'
## b = (X'X)^-1X'y 
## e = y - X(X'X)^-1X'
f1 = function(data, m, B){
  Binv = solve(t(B)%*%B, diag(m)) # inverse of matrix - (X'X)^-1
  H = B%*%Binv%*%t(B) # Hat matrix X(X'X)^-1X'
  S2 = Binv%*%t(B)
  
  ### Calculate studentized residual
  eps = data - H%*%data
  lev = diag(H)
  sigm <- sqrt((1/(length(eps)-m))*sum(eps^2))
  t_eps <- eps/(sigm*(sqrt(1-lev)))
  return(list(y_hat = H%*%data, betas = S2%*%data, eps = eps, t_eps = t_eps))
}

smooth.rgb = f1(data = log(CV2),m = n+2,B = B)

df <- data.frame(Means = log(Means), CV2 = log(CV2), 
                 fit = smooth.rgb$y_hat, eps = smooth.rgb$eps, 
                 t_eps = smooth.rgb$t_eps, B)

ggplot(data = df) + 
  geom_point(aes(Means, CV2)) +
  geom_line(aes(x = Means, y = X2), colour = brewer.pal(name = "Set2", n = 7)[2], size = 1) +
  geom_line(aes(x = Means, y = X3), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X4), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X5), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X6), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X7), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X8), colour =brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X9), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X10), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X11), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X12), colour = brewer.pal(name = "Set2", n = 7)[3], size = 1) +
  geom_line(aes(x = Means, y = X1), colour = brewer.pal(name = "Set2", n = 7)[1], size = 1) +
geom_line(aes(x = Means, y = fit), size = 3, col = "red") + ylab("log(CV2)") + xlab("log(Mean Expr)")


ggplot(data = df) + 
  geom_point(aes(x = Means, y = eps), size = 1) 

ggplot(data = df) + 
  geom_point(aes(x = Means, y = t_eps), size = 1)

plot(log10(Means), smooth.rgb$betas[1]*B[,1] + smooth.rgb$betas[2]*B[,2] +smooth.rgb$betas[3]*B[,3] +
       smooth.rgb$betas[4]*B[,4] + smooth.rgb$betas[5]*B[,5] + smooth.rgb$betas[6]*B[,6])

plot(log10(Means),  smooth.rgb$eps, pch = 16, cex = 0.5)
abline(h = 0, col = "red")

