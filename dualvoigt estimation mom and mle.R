# packages
 require(minqa) # for bobyqa
 require(fitdistrplus) # for numerical mle and mom estimates
 require(truncnorm) # as above
 
# pdf of dual voigt
ddualvoigt<-function(x,gamma,sigma){
	erfc = 2*(1-pnorm(gamma/sigma))
	dens <- sigma/(sqrt(2*pi))*exp(-(gamma^2/(2*sigma^2)))*(1/erfc)*
	exp(-(gamma*abs(x)+x^2*sigma^2/2 ))
	return(dens)
}

# sampling function for  dual Voigt (using algorithm on appendix B)
rdualvoigt <- function(n,gamma,sigma){
	rdualsample=vector()
ngood=length(rdualsample)
while(ngood<n){
	c=rnorm(1,mean=-gamma/sigma^2,sd=1/sigma)
	if (c > 0) rdualsample = c(c,rdualsample)
	if(c < -2*gamma/sigma^2) rdualsample = c(c + 2*gamma/sigma^2 ,rdualsample)
	ngood=length(rdualsample)
}	
return(rdualsample)
}
# generate data: n-sample from dual voigt
n= 5000
gamma=3
sigma=1
x=rdualvoigt( n , gamma, sigma)
# histogram of dual voigt sample
h <- hist(x,probability =TRUE,breaks=200)
h
# compare with true density
curve(ddualvoigt(x, gamma, sigma), 
      col="darkblue", lwd=3, add=TRUE, yaxt="n")
 # add normal
curve( dnorm(x,mean=0, sd=sigma),add=TRUE,lty=2, col="black" , type="l", lwd=2)

#add legend
#legend("topright", legend = c(paste("DualVoigtsample(0,", sigma,",", gamma,")"),
#paste("N(0,",sigma^2,")"),paste("Voigt(0,", sigma,",", gamma,")") ),
# lty=c(1,2,3),lwd=2,cex=1,bty="n") 

# calculate sample second moment
sum(x^2)/length(x)
# calculate sample variance
var(x)
# true variance
gamma^2/sigma^4 + 1/sigma^2 - 1/sigma*(gamma/(sigma^2*sqrt(2*pi))* exp(-gamma^2/(2*sigma^2)))/(1-pnorm(gamma/sigma))

# calculate sample fourth moment
sum(x^4)/length(x)
# true fourth moment
gamma^4/sigma^8 + 6*gamma^2/sigma^6 + 3/sigma^4 - (5*gamma/sigma^5+gamma^3/sigma^7)* 1/sqrt(2*pi)* exp(-gamma^2/(2*sigma^2))/(1-pnorm(gamma/sigma))

# parameter estimation (for previous dataset in line 30)
# 1) non linear least square fit from empirical histogram
# find bin density at each x
bin_index <- findInterval(x, h$breaks, rightmost.closed = TRUE)
dens_at_x <- h$density[bin_index]
# find bin smoothed density at each x
spline_fit <- smooth.spline(h$mids, h$density)
dens_fun_smooth <- function(xx) predict(spline_fit, xx)$y
dens_at_x_smooth <- dens_fun_smooth(x)
# non linear fit of hist density
fit_nls <- nls( dens_at_x ~ ddualvoigt(x, gamma, sigma), data = data.frame(x,dens_at_x),lower=c(0,0), start=list(gamma=2,sigma=2),algorithm="port")
# non linear fit of smoothed hist density
fit_nls <- nls( dens_at_x_smooth ~ ddualvoigt(x, gamma, sigma), data = data.frame(x,dens_at_x_smooth), start=list(gamma=2,sigma=2),algorithm="port",lower=c(0,0))
coef(fit_nls)
# setting for R2
y_pred <- predict(fit_nls)
y_obs <- dens_at_x
#  R^2
ss_res <- sum((y_obs - y_pred)^2)
ss_tot <- sum((y_obs - mean(y_obs))^2)
R2 <- 1 - ss_res/ss_tot
R2


# transform sample
# tx is a sample from T~truncnorm(0, + infty, -gamma/sigma^2, 1/sigma^2)
# the last two parameters are mean and variance of the underlying normal
tx = x[x > 0] # right half sample
tx = c( x[x > 0] , - x[x < 0] )


# 2) mle estimates of T parameters
fit = fitdist(tx, "truncnorm", method="mle", fix.arg=list(a=0,b=Inf),
    start = list(mean = mean(tx), sd = sd(tx)))
 # mle estimates of dualvoigt parameters (via eq. 20)
egamma = - fit$estimate["mean"]*(fit$estimate["sd"])^(-2)
names(egamma)<-c("gammahat")
esigma = fit$estimate["sd"]^(-1)
 names(esigma)<-c("sigmahat")
egamma; esigma

# 3) mom estimate
# moments formulae for truncnorm (order=1,2,3 and 4)
memp <- function(x, order) mean(x^order)
mtruncnorm <- function(order,a,b,mean,sd){
	alpha= (a-mean)/sd; beta=(b-mean)/sd
	if (order==1) m<- mean + sd*(dnorm( alpha ) - dnorm(beta))/ (pnorm(beta)-pnorm(alpha))
	if (beta==Inf) fac <-0 else fac <- (mean+b)*dnorm(beta)
	if (order==2) m<- mean^2 + sd^2 - sd*( fac-(mean+a)*dnorm(alpha) )/(pnorm(beta)-pnorm(alpha))
		if (order==3) m<- mean^3 + 3*mean*sd^2 - sd*( fac-(mean^2+2*sd^2+a*mean + a^2)*dnorm(alpha) )/(pnorm(beta)-pnorm(alpha))
				if (order==4) m<- mean^4 + 6*mean^2*sd^2 + 3*sd^4 - sd*( fac-(a^3+a^2*mean +a*mean^2+sd^2*(3*a+5*mean) + mean^3)*dnorm(alpha) )/(pnorm(beta)-pnorm(alpha))
	return(m)}
	# mom estimates of T parameters
	# obtained equating T1 and T2 to theoretical moments 1 and 2 (order = c(1,2))
fit = fitdist(tx, "truncnorm", method="mme", order=c(1,2),memp=memp,fix.arg=list(a=0,b=Inf),
    start = list(mean = mean(tx), sd = sd(tx)))
fit$convergence # 0= achieved convergence
# mom estimates of dual voigt via eq 20
egamma = - fit$estimate["mean"]*(fit$estimate["sd"])^(-2)
esigma = fit$estimate["sd"]^(-1) 
names(egamma)<-c("gammahat");names(esigma)<-c("sigmahat")
egamma; esigma

# 4) mle estimate using bobyqa
# negative Log-likelihood
negloglik <- function(par) {
  gamma <- par[1]
  sigma <- par[2]
  
  # Penalize out of support par values
  if (gamma <= 0 || sigma <= 0) return(1e10)
  
  dens_vals <- ddualvoigt(x, gamma, sigma)
  
  # penalize 0 or NA par values
  if (any(dens_vals <= 0 | is.na(dens_vals))) return(1e10)
  
  return(-sum(log(dens_vals)))
}
# maximize log-lik with BOBYQA
fit <- bobyqa(
  par   = c(gamma = 1, sigma = 1),  # starting values
  fn    = negloglik,
  lower = c(1e-6, 1e-6),            # par lower bound  > 0
  upper = c(10, 10)                  # par upper bound
)
fit
egamma <- fit$par[1]
esigma <- fit$par[2]
egamma;esigma

# code for simulation study in Section 5
# set parameters
gamma = 1; sigma = 1
n = 500
# set simulation parameters, boxes..
nsim = 100
gamma.mom = gamma.mle= gamma.nls <- rep(0,nsim)
sigma.mom = sigma.mle = sigma.nls <- rep(0,nsim)
r2.mom = r2.mle = r2.nls <- rep(0,nsim)
#convergence.mom <- convergence.mle <- convergence.nls <-rep(NA,nsim)
# function to be minimized by non linear least square fit
ssq_fun <- function(par) {
  gamma <- par[1]
  sigma <- par[2] 
  # penalizza valori non validi
  if (gamma <= 0 || sigma <= 0) return(1e10)
  y_pred <- ddualvoigt(x, gamma, sigma)
  # penalizza se densità negative o NA
  if (any(is.na(y_pred)) || any(y_pred <= 0)) return(1e10)
  # somma dei quadrati degli errori tra densità teorica e densità smussata
  sum((dens_at_x_smooth - y_pred)^2)
}
# MC loop
for (j in 1:nsim){
	# sample
	x=rdualvoigt( n , gamma, sigma)
    tx = x[x > 0] 
    tx = c( x[x > 0] , - x[x < 0] )
    # mom estimates
	fit.mom = fitdist(tx, "truncnorm", method="mme", order=c(1,2),memp=memp,
	fix.arg=list(a=0,b=Inf), # known truncation values
    start = list(mean = mean(tx), sd = sd(tx)))
    gamma.mom[j] = - fit.mom$estimate["mean"]*(fit.mom$estimate["sd"])^(-2) # transform via eq 20
    sigma.mom[j] = fit.mom$estimate["sd"]^(-1)  # same as above
    # mle estimates
    fit.mle = fitdist(tx, "truncnorm", method="mle", fix.arg=list(a=0,b=Inf),
    start = list(mean = mean(tx), sd = sd(tx)))
    gamma.mle[j] = - fit.mle$estimate["mean"]*(fit.mle$estimate["sd"])^(-2)
    sigma.mle[j] = fit.mle$estimate["sd"]^(-1)	
    # nls estimates
    h <- hist(x,probability =TRUE,breaks=200)
    # find bin density at each x
    bin_index <- findInterval(x, h$breaks, rightmost.closed = TRUE)
    # dens_at_x <- h$density[bin_index]
    # find bin smoothed density at each x
    spline_fit <- smooth.spline(h$mids, h$density)
    dens_fun_smooth <- function(xx) predict(spline_fit, xx)$y
    dens_at_x_smooth <- dens_fun_smooth(x)
    fit_nls <- bobyqa(
      par   = c(1.5, 1.5),         # starting values
      fn    = ssq_fun,
     lower = c(1e-6, 1e-6),       # lower lim
      upper = c(10, 10)            # upper lim
     )
# non linear fit of smoothed hist density and R^2
# fit_nls <- nls( dens_at_x_smooth ~ ddualvoigt(x, gamma, sigma), data = data.frame(x,dens_at_x), start=list(gamma=2,sigma=2))
   gamma.nls[j] <- fit_nls$par[1]
   sigma.nls[j] <- fit_nls$par[2]
   dens_model_best <- ddualvoigt(x, fit_nls$par[1], fit_nls$par[2])
RSS <- sum((dens_at_x_smooth - dens_model_best)^2)
TSS <- sum((dens_at_x_smooth - mean(dens_at_x_smooth))^2)
r2.nls[j] = 1- RSS/TSS
}
# average point estimate and sd over nsim
mean(gamma.mom);sd(gamma.mom)
mean(sigma.mom);sd(sigma.mom)
gamma; sigma
mean(gamma.mle);sd(gamma.mle)
mean(sigma.mle);sd(sigma.mle)
gamma;sigma
mean(gamma.nls);sd(gamma.nls)
mean(sigma.nls);sd(sigma.nls); mean(r2.nls)
gamma;sigma
# Table 1 (the n^th row)
# average bias 
mean(abs(gamma.mle-gamma));mean(abs(gamma.mom-gamma));mean(abs(gamma.nls-gamma))
mean(abs(sigma.mle-sigma));mean(abs(sigma.mom-sigma)); mean(abs(sigma.nls-sigma))
# average montecarlo sd of estimates
mean(sd(sigma.mle),sd(gamma.mle))
mean(sd(sigma.mom),sd(gamma.mom))
mean(sd(sigma.nls),sd(gamma.nls))

