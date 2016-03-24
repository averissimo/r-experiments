#
# Taylor expansion example

# set the xdata to use in taylor prediction
xdata <- seq(-2, 2, 0.1)

# original function
f <- function(e) { 
  return(exp(e))
}

# partial taylor expansion of function for each n
t <- function(x,a,n) {
  return(exp(a) / factorial(n) * (x - a)^n)
}

# full taylor expansion
h <- function(x,a,k) { 
  sum_v <- 0
  for (i in 0:k ) {
    sum_v <- sum_v + t(x,a,i)
  }
  return(sum_v)
}

# after n = 5 it starts to 
# show plots of approximation
len <- 1:9
par(mfrow=c(2,2))

for ( i in len) {
  plot(xdata, f(xdata)    , xlim = c(-2,2), ylim = c(0, f(xdata[length(xdata)]) ), type = 'l', ylab = '', xlab ='', col = 'red')
  par(new=T)
  plot(xdata, h(xdata,0,i), xlim = c(-2,2), ylim = c(0, f(xdata[length(xdata)]) ), ylab = 'y', xlab = 'x')
  title(paste0('n = ', i, ' (residuals = ', sum(abs(f(xdata) - h(xdata,0,i))) ,')'))
}
