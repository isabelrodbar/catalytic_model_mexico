###EXAMPLE CODE TO ESTIMATE AVERAGE FORCES OF INFECTION_P

#### SOURCE FUNCTIONS
#This line of code loads functions that will be needed to run the analysis. 
	# likeli.cons contains the likelihood function
	# fit.data fits the model
	# plot.fit.cons produces plots of fit for this average (across serotypes) model
	# likeli.dat Simple likelihood ratio test to compare two (nested)

source("~/Dropbox/Mexico/Code/catalytic_model_ruth.R") # Change this to the appropriate directory


#### Now we need to load the data for the analysis
dat<-read.csv("~/Dropbox/Mexico/Datos/Informacion_sujetos_prntcat2.csv")

dat$seropre2<-NA #Create a new variable to recode seroprevalence
dat$seropre2[which(dat$seropre==0 | dat$seropre=="Indeterminado")]<-0 # "indeterminados" are treated as negatives
dat$seropre2[which(dat$seropre==1)]<-1 # "1"s are treated as seropositives

dat2<-dat[which(!is.na(dat$seropre2)),] #Generate a new dataset excluding those observations that were not tested


###ALL DATA

# First perform the analysis for the whole dataset
loc.all<-table(dat2$edad, dat2$seropre2)

#fit.data is the function that fits the model to the data. It takes several arguments
# data= dataset (table with ages on the rows and 0, 1 on the columns)
# npar= number of parameters to fit. Use 1 for a time constant model and >1 to allow for time varying parameters. Divides all years in equal intervals

# type relates to the likelihood function to use. The options are "average" to use  likeli.cons or "serotype" to yse likeli.prnt depending of the data that is provided
#year- Year of the seroprevalence study (2011.5). Optional parameter

fit.all<-fit.data( data=loc.all, npar=1, type="average")

#Now allow for two parameters (change will occur at the mid point)

fit.all2<-fit.data( data=loc.all, npar=5, type="average")

#Explore fit of models
#plot fit function produces a plot with the data and the fitted line. 

plot.fit.cons(fit.all, loc.all)

##Likelihood ratio test
#likeli.rat does a simple likelihood ratio test to see if two nested models are significantly different. For example comparing a model with one and two parameters
likeli.rat(fit.all, fit.all2)






