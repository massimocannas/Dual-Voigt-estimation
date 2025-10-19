# function for calculating density of (zero centered dual Voigt)
ddualvoigt<-function(x,gamma,sigma){
	erfc = 2*(1-pnorm(gamma/sigma))
	dens <- sigma/(sqrt(2*pi))*exp(-(gamma^2/(2*sigma^2)))*(1/erfc)*
	exp(-(gamma*abs(x)+x^2*sigma^2/2 ))
	return(dens)
}

# function for random sample from (zero centered) dual Voigt
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
# example: n size sample from dual voigt
n= 1000
gamma=3
sigma=4
x=rdualvoigt( n , gamma, sigma)
# histogram of dual voigt sample
hist(x,probability =TRUE)
# superimpose pdf of dual voigt density
curve(ddualvoigt(x, gamma, sigma), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

 # add normal?
lines( x, dnorm(x,mean=0, sd=sigma),lty=2, col="black" , type="l", lwd=1)

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
gamma^4/sigma^8 + 6*gamma^2/sigma^6 + 3/sigma^4 - (-5*gamma/sigma^5+gamma^3/sigma^7)* exp(-gamma^2/(2*sigma^2))/(1-pnorm(gamma/sigma))



