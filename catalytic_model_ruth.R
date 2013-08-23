

#LIKELIHOOD FUNCTION AND OPTIMIZATION CODE 

likeli.cons<- function(theta, data, npar) {
	
	#####    Get relevant data
	age<-1:(max(as.numeric(dimnames(data)[[1]])))
	age.dat<-as.numeric(dimnames(data)[[1]])
	nyk<-data[,which(dimnames(data)[[2]]==1)]
	N<-apply(data, 1, sum)
	nxk<-N-nyk
		
	#######Create a vector (or matrix) for your lambdas. This npar=1 for constant foi, npar>1 for time varying foi (divided equally in npar parts)	
	if (npar==1) {
	lambda<-matrix(theta, nrow=length(age), ncol=1)
	lambda[1,]<-theta*1.5 #To correct for the fact that this is the floor of the ages 
	}
	
	if (npar>1) {
		cut.seq<-floor(seq(min(age), max(age), length.out=npar+1))
		cut.ages<-as.numeric(cut(age, cut.seq, include.lowest=TRUE))	
		lambda<-matrix(NA, nrow=length(age), ncol=1)
		lambda[,1] <-theta[cut.ages]
		lambda[1,]<-theta[1]*1.5
	}
	
#####Compute cumulative FOI experienced up to age a (cum.lambda)
	xaprod2<-cumsum(lambda)
	
	x<-exp(-xaprod2)	 #Proportion that remains susceptible at age a
			
 # Get proportions susceptible for observed age groups
	x.lik<-x[which(is.element(age, age.dat))]
	y.lik<-1-x.lik
	
	#Add the log likelihood terms. This is a binomial likelihood of the form x^nxk * (y)^(N-nxk)
	
	lognxk<-nxk*log(x.lik) 
	lognyk<-(N-nxk)*log(y.lik)
	logterms<-lognxk+lognyk
	#Sum the two log likelihoods
	minuslogl<- -sum(logterms)
	return(minuslogl)
	}
	
	
#Fit the model	
fit.data <- function (data, npar, type="average", Hes=FALSE, year,...) {
	#Parscale: a technicality for keeping the guesses in reasonable bounds

initial.guess<-abs(rnorm(npar, .02, 0.01))

if(type=="serotype") {initial.guess<-c(initial.guess, 1)}

parscale <- 1e-5 + abs(initial.guess)

Lik.Fun<-switch(type, 
	"average"=likeli.cons,
	"serotype"=likeli.serotype)
	


data.ok<-ifelse((type=="serotype" & class(data)!="list"), FALSE, TRUE)
	
	
if(data.ok==TRUE) {
  opt.vals <- optim (par=initial.guess, 
  					fn=Lik.Fun,   #Which likelihood function to use likeli.cons or likeli.serotype
  					data=data,    #dataset
  					npar=npar,	#Number of parameters to fit
  					method='L-BFGS-B', 
  					lower = rep(0.000001, length(initial.guess)),  #Constrain all lambdas to be greater than 0.000001
    				control=list(trace=TRUE, parscale=parscale,  maxit=100000),  #Some tuning stuff
    				hessian=Hes, #Whether to estimate the hessian Matrix or not. Default set to FALSE
    				...)
  
opt.vals$init <- initial.guess

if(type=="average") {
opt.vals$class<-paste("type= ", type, ": Model assuming constant or varying forces of infection. Same lambda for all serotypes") }

if(type=="serotype" & npar>1) {
opt.vals$class<-paste("type= ", type, ": Model assuming lambda varies between serotypes") }

if(type=="serotype" & npar==1) {
opt.vals$class<-paste("type= ", type, ": Model assuming constant or varying forces of infection. Same lambda for all serotypes") }

opt.vals$type=type
		
		if(npar>1 &  exists("year", inherits=FALSE, mode="numeric")==TRUE & type=="average"){ 
			
		age<-1:(max(as.numeric(dimnames(data)[[1]])))
		age2<-age+.5
		seq.years<-floor(seq(2011.5-min(age2), 2011.5-max(age2+1), length.out=npar+1))
		seq.years[1]<-year
		opt.vals$year.int<-seq.years
		
		if(year != 2011.5) {opt.vals$NOTE<-"year.int needs to be interpreted with caution"}
		
		}
		}

if(data.ok==FALSE) {
	print("ERROR: For serotype specific estimates, data must be presented as a list")}		

return(opt.vals)

}


plot.fit.cons<- function(model, data, col.plot="blue") {
	
	
	if(model$type=="average") {
	theta<-model$par
	npar<-length(theta)
	
	age<-1:(max(as.numeric(dimnames(data)[[1]])))
	age.dat<-as.numeric(dimnames(data)[[1]])
	nyk<-data[,which(dimnames(data)[[2]]==1)]
	N<-apply(data, 1, sum)
	nxk<-N-nyk
		
	#Create a vector (or matrix) for your lambdas. This npar=1 for constant foi, npar>1 for time varying foi (divided equally in npar parts)	
	if (npar==1) {
	lambda<-matrix(theta, nrow=length(age), ncol=1)
	lambda[1,]<-theta*1.5 #To correct for the fact that this is the floor of the ages 
	}
	
	if (npar>1) {
		cut.seq<-floor(seq(min(age), max(age), length.out=npar+1))
		cut.ages<-as.numeric(cut(age, cut.seq, include.lowest=TRUE))	
		lambda<-matrix(NA, nrow=length(age), ncol=1)
		lambda[,1] <-theta[cut.ages]
		lambda[1,]<-theta[1]*1.5
		
			#Cumulative FOI upt to age a (cum.lambda)
	xaprod2<-cumsum(lambda)
	}
	}
	
	
if(model$type=="serotype") {
	
	age<-1:(max(as.numeric(dimnames(data[[1]])[[1]])))
	age.dat<-as.numeric(dimnames(data[[1]])[[1]])
	
	###Get relevant data from IgG
    nyk<-data[[1]][,which(dimnames(data[[1]])[[2]]==1)]
	N<-apply(data[[1]], 1, sum)
	nxk<-N-nyk
	
npar<-length(model$par)-1
	#Arrange lambda matrix. if npar=1, it will fix all of the serotypes to have the same lambda. If npar=2 it will fix serotypes 1-3 to have the same lambda, and will allow 4 to vary. If npar=4, each serotype will have it's own lambda.
	theta<-model$par[1:npar]
	lambda<-matrix(theta[1:npar], length(age), 4, byrow=TRUE)
	lambda[1,]<-theta[1:npar]*1.5
	
	if(npar==2) { 
		lambda<-matrix(NA, length(age), 4, byrow=TRUE)
		lambda[,1:3]<-theta[1]
		lambda[,4]<-theta[2]
		lambda[1,]<-lambda[1,]*1.5
}
	
	xaprod<-matrix(0, length(age),4)
	xaprod[1,]<-lambda[1,]*age[1]
	for(i in 2:length(age)) {xaprod[i,]<-lambda[i,]		}
	xaprod2<-matrix(NA, length(age), 4)
	xaprod2[1,]<-xaprod[1,]
	
	lambda.working<-lambda

	#Get cumulative lambdas up to age a	
		for (i in 2:length(age)) {
		xaprod2[i,]<-colSums(xaprod[1:i,])}	
		xaprod2<-apply(xaprod2, 1, sum)
		
		}
     

	
	x<-exp(-xaprod2)	 #Proportion susceptible
	y<- 1-x   #Proportion immune
		
	#Add the log likelihood terms. This is a binomial likelihood of the form x^nxk * (y)^(N-nxk)
	
	prop.seropos<-nyk/N
	
	plot(age.dat, prop.seropos, ylim=c(0,1), pch=19, ylab="Proporiton Seropositive", xlab="Age (years)", col="grey60", xlim=c(min(age), max(age)))
	lines(age, y, col=col.plot)	
	
	}
		
	
likeli.rat<-function(model1, model2) {
	
	val1<-model1$val
	val2<-model2$val
	
	par.diff<-abs(length(model1$par)-length(model2$par))
	
	if(par.diff==0) { cat("Models have same models of parameters. Can't compute likelihood ratio")
							likeli.rat<-NA}
	
	if(par.diff>0) {likeli.rat<-1-pchisq(abs(2*(val1-val2)), par.diff)} 
	
	return(paste("Likelihood ratio test=", likeli.rat))
	
}	


####Likelihood for serotype specific models
likeli.serotype<- function(theta, data, npar) {

#Get relevant data. for this model, data must be entered as a list. First element is binary data from IgG, second element is data from PRNT
	age<-1:(max(as.numeric(dimnames(data[[1]])[[1]])))
	age.dat<-as.numeric(dimnames(data[[1]])[[1]])
	
	###Get relevant data from IgG
    nyk<-data[[1]][,which(dimnames(data[[1]])[[2]]==1)]
	N<-apply(data[[1]], 1, sum)
	nxk<-N-nyk
	
	
	#Get relevant data from PRNT
	nik<-matrix(0, length(age.dat), ncol=4)
	m.prnt<-match(dimnames(data[[2]])[[2]], 1:4)
	m.prnt<-m.prnt[which(!is.na(m.prnt))]
	nik[, m.prnt]<-data[[2]][,which(is.element(dimnames(data[[2]])[[2]], c(1,2,3,4)))]
	nmk<-data[[2]][,which(dimnames(data[[2]])[[2]]==5)]
	
	
	#Arrange lambda matrix. if npar=1, it will fix all of the serotypes to have the same lambda. If npar=2 it will fix serotypes 1-3 to have the same lambda, and will allow 4 to vary. If npar=4, each serotype will have it's own lambda.
	
	lambda<-matrix(theta[1:npar], length(age), 4, byrow=TRUE)
	lambda[1,]<-theta[1:npar]*1.5
	
	if(npar==2) { 
		lambda<-matrix(NA, length(age), 4, byrow=TRUE)
		lambda[,1:3]<-theta[1]
		lambda[,4]<-theta[2]
		lambda[1,]<-lambda[1,]*1.5
}
	
	#This is the interaction parameter
	interpar<-theta[((npar+1):length(theta))]
	interpar<-matrix(interpar, ncol=4)
	xaprod<-matrix(0, length(age),4)
	xaprod[1,]<-lambda[1,]*age[1]
	for(i in 2:length(age)) {xaprod[i,]<-lambda[i,]		}
	xaprod2<-matrix(NA, length(age), 4)
	xaprod2[1,]<-xaprod[1,]
	
	lambda.working<-lambda

	#Get cumulative lambdas up to age a	
		for (i in 2:length(age)) {
		xaprod2[i,]<-colSums(xaprod[1:i,])}	
     
     
     #Compute proportion susceptible from the sum of the forces of infection across all serotypes
     lnx<-apply(xaprod2,1,sum)
     x<-exp(-lnx)		
	
	
	#This part of the code computes the probability of being monotypic at age a
	zi<-matrix(NA, nrow=length(age), ncol=4)
	for (primary in 1:4) {
	zimat<-matrix(NA, nrow=(length(age)), ncol=4)
	interpar2<-rep(interpar[primary],4) #Model 6a
	#interpar2<-interpar #Model 6b
	
	interpar2[primary]<-0
	interpari<-sweep(lambda.working, 2, (1-interpar2), FUN="*")
	
		zimat[1,]<-interpari[1,]
	for(w in 2:nrow(interpari)) {
		zimat[w,] <-colSums(interpari[1:w,])}

	sumzimat<-apply(zimat, 1, sum)	
			sumzimat<-c(0, sumzimat[1:(length(age)-1)])
	
	zimat13<-exp(sumzimat)		
	zimatlambda<-zimat13*lambda.working[,primary]
	
	zi[1,primary]<-zimatlambda[1]
	for (ag in 2:length(age)) {
	zi[ag, primary]<-sum(zimatlambda[1:ag])}
	
}

zi.fin<-sweep(zi, 1, x, "*")		
	

			
 # Get proportions susceptible and m onotypic for observed age groups
	x.lik<-x[which(is.element(age, age.dat))]
	y.lik<-1-x.lik
	zi.fin.lik<-zi.fin[which(is.element(age, age.dat)),]
	
	###Likelihood for exposed vs unexposed (binomial likelihood)
		lognxk<-nxk*log(x.lik)
		lognyk<-nyk*log(y.lik)
		
	###Likelihood for prnt data (conditional multimomial likelihood)
	zisum<-apply(zi.fin.lik,1,sum)
	lognik<-nik*log(zi.fin.lik/(1-x.lik))	
	lognmk<-nmk*log(((1-x.lik-zisum)/(1-x.lik)))
	logterms<-sum(lognxk) + sum(lognyk) + sum(lognik)+sum(lognmk)
	minuslogl<- -sum(logterms)
	#print(minuslogl)
	return(minuslogl)
	
	}






