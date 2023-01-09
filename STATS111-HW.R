#####################################################
# First, install the package called "epitools" and load it

# Need to only install this package once
install.packages('epitools')
install.packages('rmeta')
install.packages('pROC')


#####################################################
# Need to always load the package though
library(epitools)
library(rmeta)
library(pROC)
library(nnet)

#####################################################
#####################################################
# This is the prop.comp() function we will use
#######
prop.comp <- function( x, estimate="all", conf.level=.95, transpose=FALSE ){
	if( transpose ) x <- t(x)
	rslt <- vector( "list", length=3 )
	names( rslt ) <- c( "riskdiff", "riskratio", "oddsratio" )
	diff.rslt <- suppressWarnings(prop.test( x, conf.level=conf.level ))
	rslt[[1]] <- rslt[[2]] <- rslt[[3]] <- epitab( x, method="riskratio", pvalue="chi2", conf.level=conf.level )$tab
	colnames( rslt[[1]] )[5] <- "riskdiff"
	rslt[[1]][,5] <- c(0,diff(rev(diff.rslt$estimate)))
	rslt[[1]][2,6:7] <- diff.rslt$conf.int
	colnames( rslt[[3]] )[5] <- "oddsratio"
	rslt[[3]][,5:8] <- suppressWarnings(epitab( x, method="oddsratio", pvalue="chi2", conf.level=conf.level )$tab[,5:8])
	if(is.null(names(dimnames(x)))){
		for(i in 1:3){
			colnames(rslt[[i]])[c(1,3)] <- c("Outcome=0", "Outcome=1")
			rownames(rslt[[i]]) <- c("Group=1", "Group=2")
			}
	}
	if( is.element( estimate, c("all", "oddsratio") ) ){ 
		if(is.null(names(dimnames(x)))){
			warning( "Estimated probabilities represent Pr[ Outcome | Group ]. For estimates of 
			Pr[ Group | Outcome ], change the value of 'transpose'.")
		}
		else
			warning( paste("Estimated probabilities represent Pr[", names(dimnames(x))[2], 
			"|",names(dimnames(x))[1], "]. For estimates of 
			Pr[", names(dimnames(x))[1], "|",names(dimnames(x))[2], "], change the value of 'transpose'.") )
		}
	if( estimate == "riskdiff" ) return(rslt[[1]])
	else if( estimate == "riskratio" ) return(rslt[[2]])
	else if( estimate == "oddsratio" ) return(rslt[[3]])
	else return(rslt)
}



#####################################################
#####################################################
#####################################################

# This is the data for female and CHD example
sex.chd = matrix(c(1226,2000,823,650),2,2)
rownames(sex.chd) = c("Female no","Female yes")
colnames(sex.chd) = c("CHD no","CHD yes")

# Now use the prop.comp function on this data
prop.comp(sex.chd)


# This is the data for the blood pressure and CHD data
blood.chd = matrix(c(1663,978,395,190,534,503,255,181),4,2)
rownames(blood.chd) = c("low","above normal","high","very high")
colnames(blood.chd) = c("CHD no","CHD yes")

# Now do a chi square test of independence
chisq.test(blood.chd , correct=FALSE)

#####################################################
# The following several lines will prepare the data needed for homework 1

# Data for fiber questions will be called "fiber"
fiber = matrix(c(140,200,60,50),2,2)
rownames(fiber) = c("Low fiber","High fiber")
colnames(fiber) = c("Disease no","Disease yes")

# Data for smoking example will be callled "smoke.school"
smoke.school = matrix(c(1168,1823,1380,188,416,400),3,2)
rownames(smoke.school) = c("0 parents smoke","1 parent smokes","2 parents smoke")
colnames(smoke.school) = c("Smoke no","Smoke yes")

#####################################################
#####################################################
#####################################################

# Read in the framingham data
framingham = read.table("./framingham.txt")

# Convert sex to 0 or 1 and create a female variable
framingham$sex = framingham$sex - 1
framingham$female =framingham$sex


# Create a new variable, which categorizes the blood pressure of a subject
framingham$sbpgrp = cut( framingham$sbp,breaks=c(min(framingham$sbp),126,146,166, max(framingham$sbp)),include.lowest=TRUE )

# Create 4x2 table of chdfate and sbpgrp
chd.table = xtabs(~sbpgrp+chdfate,data=framingham)

# Trend test: manual version
epitab( chd.table, pvalue="chi2" )
r = cor(framingham$chdfate, as.numeric(framingham$sbpgrp))
n = dim(framingham)[1]
M2 = (n-1)*r^2
cbind(n,r,M2)
1- pchisq( M2, df=1 )

# Trend test: using prop.trend.test()
chd.table
n.chd = chd.table[,2]
n.strata = rowSums(chd.table) 
prop.trend.test( n.chd, n.strata )

# Create 2x2x4 table of sbp group, sex, and chd fate
chd.female.table =table(framingham$sex,framingham$chdfate, framingham$sbpgrp)
dimnames(chd.female.table)= list(chdfate=c("chdfate yes", "chdfate no"), female=c("female yes","female no"), sbpgrp=c("low","normal","high","very high"))

# Compute relevant counts to conduct odds ratio heterogeneity test
n.male = table( framingham$female, framingham$sbpgrp )[1,]
n.female = table( framingham$female, framingham$sbpgrp )[2,]
m.chd = table( framingham$chdfate, framingham$female, framingham$sbpgrp )[2,1,]
f.chd= table( framingham$chdfate, framingham$female, framingham$sbpgrp )[2,2,]

# Conduct Maentel Haenszel test
mantelhaen.test(chd.female.table)

# Test for heterogeneity with the Breslow-Day test
mh.rslt = meta.MH( n.female, n.male, f.chd, m.chd, names=levels(framingham$sbpgrp) )
summary( mh.rslt )

#####################################################
# Create death penalty data
d.penalty = array( c( 415,38,54,12,17,140,1,5),dim = c(2, 2, 2),dimnames = list( Defendant.Race = c("White","Black"),Verdict = c("No.Death", "Death"),Victim.Race = c("White","Black")))

# Combine the two tables to create a single 2x2 table
marg.data = d.penalty[,,1]+d.penalty[,,2]

# Conduct test on this combined table
prop.comp(marg.data, estimate="oddsratio")

# Conduct Maentel Haenszel test
mantelhaen.test(d.penalty)

# Conduct prop.comp function on each level of W (victim race).
prop.comp(d.penalty[,,1], estimate="oddsratio")
prop.comp(d.penalty[,,2], estimate="oddsratio")

# Test for heterogeneity with the Breslow-Day test
n.black = c( sum( d.penalty[,,1][2,] ), sum( d.penalty[,,2][2,] ) ) 
n.white = c( sum( d.penalty[,,1][1,] ), sum( d.penalty[,,2][1,] ) ) 
b.death = c( sum( d.penalty[,,1][2,2] ), sum( d.penalty[,,2][2,2] ) ) 
w.death = c( sum( d.penalty[,,1][1,2] ), sum( d.penalty[,,2][1,2] ) )

mh.rslt = meta.MH( n.black, n.white, b.death, w.death, names=c("White Victim", "Black Victim" ) )
summary( mh.rslt )

#####################################################
#####################################################
#####################################################
# Homework code

# This is the 3 way scouts table
scout = array( c(169,43,42,11,59,196,10,10),dim = c(2, 2, 2),dimnames = list( Scout = c("No","Yes"),Verdict = c("Delinquet No", "Delinquet Yes"),Socioeconomic = c("Low","High")))

# This is the two way table of Socioeconomc status and Scout stats
socio = array(c(211,69,54,206), dim=c(2,2),dimnames=list(Socioeconomic=c("Low","High"), Scout=c("No","Yes")))

# This is the two way table of Socioeconomic stats and Deliquency status
socio.deliq = array(c(212,255,53,20), dim=c(2,2),dimnames=list(Socioeconomic=c("Low","High"), Deliquent=c("No","Yes")))

# This the two way table ignoring Socioeconomic status
scout.margin = scout[,,1]+scout[,,2]

# Ordinal data from homework
cholesterol = array(c(307,246,439,245,12,8,31,41), dim=c(4,2), dimnames=list( Cholesterol=c("Normal","Above Normal","High","Very High"), CHD=c("no CHD","CHD")))

# Need to use these two as arguments in the prop.trend.test() function.
n.chol = cholesterol[,2]
n.strata = rowSums(cholesterol)

# Framingham data
framingham = read.table("./framingham.txt")

#####	Recode sex to something obvious (sex=1 -> male)
##
framingham$sex = framingham$sex - 1
names( framingham )[1] = "female"
framingham[1:10,]


##
#####	Create SBP and BMI groups
##
framingham$sbphi = cut( framingham$sbp, breaks=c(min(framingham$sbp),146, max(framingham$sbp)), include.lowest=TRUE )
framingham$bmigrp = cut( framingham$bmi, breaks=c(min(framingham$bmi, na.rm=TRUE),20, 25, 30, max(framingham$bmi, na.rm=TRUE)), include.lowest=TRUE, right=FALSE )

##
#####	Compute test of independence and test for trend for BMI and SBP
bmisbp.table = xtabs( ~ bmigrp + sbphi, data=framingham )
epitab( bmisbp.table, pvalue="chi2" )
n.hisbp = bmisbp.table[,2]
n.strata = rowSums(bmisbp.table) 
chisq.test(bmisbp.table)
prop.trend.test( n.hisbp, n.strata )

##
#####	Compute test of independence and test for trend for BMI and CHD
bmichd.table = xtabs( ~ bmigrp + chdfate, data=framingham )
epitab( bmichd.table, pvalue="chi2" )
n.chd = bmichd.table[,2]
n.strata = rowSums(bmichd.table) 
chisq.test(bmichd.table)
prop.trend.test( n.chd, n.strata )

##
#####	Compute relevant counts
##
n.sbplo = xtabs( ~ sbphi + bmigrp, data=framingham )[1,]
n.sbphi = xtabs( ~ sbphi + bmigrp, data=framingham )[2,]
lo.chd = xtabs( ~ chdfate + sbphi + bmigrp, data=framingham )[2,1,]
hi.chd = xtabs( ~ chdfate + sbphi + bmigrp, data=framingham )[2,2,]

##
#####	Compute M-H estimate of adjusted OR and test for heterogeneity
##
mh.rslt = meta.MH( n.sbphi, n.sbplo, hi.chd, lo.chd, names=levels(framingham$bmigrp) )
summary( mh.rslt )



#################################################
#################################################
# Regression review

# Import the Midwest house price data
house = read.table("./MidwestSales.txt", fill=TRUE, header=FALSE)

# This dataset does not have names, so we will add names to the variables
names(house)=c("id","price","sqft","bed","bath","ac","garage","pool","year","quality","style","lot","hwy")

# Fit a model with sqft, lot, and ac as explanatory variables and price as response
summary(lm(price~sqft+lot+ac+sqft*lot, data=house))

# Test if lot should be included in the model
full = lm(price~sqft+lot+ac+sqft*lot, data=house)
reduced = lm(price~sqft+ac, data=house)
anova(reduced,full)


#####################################################
# Outliers

crime = read.table("./crime.txt", header=TRUE)

# Fit adjusted and unadjusted models
fit.unadj = lm( crime ~ single, data=crime )
summary( fit.unadj )

fit.adj = lm ( crime ~ single + pctmetro + pcths + poverty, data=crime )
summary( fit.adj )

##### Calculation of standardized residuals
r = rstandard( fit.adj )
cbind( crime,r )[1:10,]
plot( crime$crime, r, ylim=range(r)*1.1, xlim=range(crime$crime)*1.1, 
		xlab="Crime rate (per 100,000/yr)", ylab="Standardized Residual" ) 
text( crime$crime+100, r+.2, labels=crime$state )
abline( h=c(-2.5,2.5), col="red", lwd=2 )W20cbind(crime,r)[ abs(r) > 2, ]

##### Calculation of leverage and plot of leverage by residual
lev = hatvalues( fit.adj )
cbind(crime, lev)[ lev>(2*4+2)/51, ]

plot( crime$crime, lev, ylim=range(lev)*1.1, xlim=range(crime$crime)*1.1,
		xlab="Crime rate (per 100,000/yr)", ylab="Leverage" ) 
text( crime$crime+100, lev+.02, labels=crime$state )
abline( h=(2*4+2)/51, col="red", lwd=2 )

#################################################
#################################################
#################################################
# GLM

# Defining logit and expit functions
expit = function(x) exp(x)/(1+exp(x))
logit = function(x) ifelse(x/(1-x)<0,'error',log(x/(1-x)))

# Read in the framingham data
framingham = read.table("./framingham.txt")

# Convert sex to 0 or 1 and create a female variable
framingham$sex = framingham$sex - 1
framingham$female =framingham$sex

# Create a new variable, which categorizes the blood pressure of a subject
framingham$sbpgrp = cut( framingham$sbp,breaks=c(min(framingham$sbp),126,146,166, max(framingham$sbp)),include.lowest=TRUE )

# Read in the fev data
fev = read.table("./fev.txt", header=TRUE)
fev$smoker = ifelse(fev$smoke=="nosmoker",0,1)

#####################################################
# Import the Midwest house price data
house = read.table("./MidwestSales.txt", fill=TRUE, header=FALSE)

# This dataset does not have names, so we will add names to the variables
names(house)=c("id","price","sqft","bed","bath","ac","garage","pool","year","quality","style","lot","hwy")

# Create large indicator
house$large = ifelse(house$sqft>3000,1,0)

#####################################################
# Read in the wcgs data
wcgs = read.csv("./wcgs.csv", header=TRUE)

#####################################################
# Read in MCAT data
mcat = read.table("./MedGPA.txt", fill=TRUE, header=TRUE)

#################################################

# Lecture note models

m1 = glm(Acceptance~Sex, family=binomial(link="logit"), data=mcat)
summary(m1)

m2 = glm(Acceptance~MCAT, family=binomial(link="logit"), data=mcat)
summary(m2)

m3 = glm(Acceptance~Sex+MCAT+Sex*MCAT, family=binomial(link="logit"), data=mcat)
summary(m3)

m4 = glm(chdfate~as.factor(sbpgrp), family=binomial(link="logit"),data= framingham)
summary(m4)

m5 = glm(chdfate~female+as.factor(sbpgrp)+female*as.factor(sbpgrp), family=binomial(link="logit"),data= framingham)
summary(m5)

m6 = glm(chdfate~sbp+age, family=binomial(link="logit"),data= framingham)
summary(m6)


m6.1 = glm(chdfate~sbp+age+sbp*age, family=binomial(link="logit"),data= framingham)
summary(m6.1)


# Homework models
m7 = glm(ac~pool, family=binomial(link="logit"), data=house)
summary(m7)


m8 = glm(ac~pool+sqft+pool*sqft, family=binomial(link="logit"), data=house)
summary(m8)


m8.1 = glm(ac~pool+sqft+lot, family=binomial(link="logit"), data=house)
summary(m8.1)


m9 = glm(chd~as.factor(smoke), family=binomial(link="logit"), data=wcgs)
summary(m9)

m10 = glm(chd~as.factor(smoke)+bp+bp*as.factor(smoke), family=binomial(link="logit"), data=wcgs)
summary(m10)


#####################################################
# RUN THESE FUNCTIONS ####################

# Ifelse function

ifelse1 =function(test, x, y){ if (test) x else y}

#####################################################
#  Function to exponentiate coefficients and produces CIs for GLMs

glmCI <- function( model, transform=TRUE, robust=FALSE ){
	link <- model$family$link
	coef <- summary( model )$coef[,1]
	se <- ifelse1( robust, robust.se.glm(model)[,2], summary( model )$coef[,2] )
	zvalue <- coef / se
	pvalue <- 2*(1-pnorm(abs(zvalue)))

	if( transform & is.element(link, c("logit","log")) ){
		ci95.lo <- exp( coef - qnorm(.975) * se )
		ci95.hi <- exp( coef + qnorm(.975) * se )
		est <- exp( coef )
	}
	else{
		ci95.lo <- coef - qnorm(.975) * se
		ci95.hi <- coef + qnorm(.975) * se
		est <- coef
	}
	rslt <- round( cbind( est, ci95.lo, ci95.hi, zvalue, pvalue ), 4 )
	colnames( rslt ) <- ifelse1( 	robust, 	
					c("Est", "robust ci95.lo", "robust ci95.hi", "robust z value", "robust Pr(>|z|)"),
					c("Est", "ci95.lo", "ci95.hi", "z value", "Pr(>|z|)") )			
	colnames( rslt )[1] <- ifelse( transform & is.element(link, c("logit","log")), "exp( Est )", "Est" )
	rslt
	}
	
#####################################################
#	Function to estimate linear contrasts of coefficients from a GLM fit

linContr.glm <- function( contr.names, contr.coef=rep(1,length(contr.names)), model, transform=TRUE ){
	beta.hat <- model$coef 
	cov.beta <- vcov( model )

	contr.index <- match( contr.names, dimnames( cov.beta )[[1]] )	
	beta.hat <- beta.hat[ contr.index ]
	cov.beta <- cov.beta[ contr.index,contr.index ]
	est <- contr.coef %*% beta.hat
	se.est <- sqrt( contr.coef %*% cov.beta %*% contr.coef )
	zStat <- est / se.est
	pVal <- 2*pnorm( abs(zStat), lower.tail=FALSE )
	ci95.lo <- est - qnorm(.975)*se.est
	ci95.hi <- est + qnorm(.975)*se.est
	
	link <- model$family$link
	if( transform & is.element(link, c("logit","log")) ){
		ci95.lo <- exp( ci95.lo )
		ci95.hi <- exp( ci95.hi )
		est <- exp( est )
		cat( "\nTest of H_0: exp( " )
		for( i in 1:(length( contr.names )-1) ){
			cat( contr.coef[i], "*", contr.names[i], " + ", sep="" )
			}
		cat( contr.coef[i+1], "*", contr.names[i+1], " ) = 1 :\n\n", sep="" )		
		}
	else{
		cat( "\nTest of H_0: " )
		for( i in 1:(length( contr.names )-1) ){
			cat( contr.coef[i], "*", contr.names[i], " + ", sep="" )
			}
		cat( contr.coef[i+1], "*", contr.names[i+1], " = 0 :\n\n", sep="" )
		}
	rslt <- data.frame( est, se.est, zStat, pVal, ci95.lo, ci95.hi )
	colnames( rslt )[1] <- ifelse( transform && is.element(link, c("logit","log")), "exp( Est )", "Est" )
	round( rslt, 8 )
}

#####################################################
# Function to compute deviance (LR) test p-Value

lrtest <- function( fit1, fit2 ){
	cat( "\nAssumption: Model 1 nested within Model 2\n\n" )
	rslt <- anova( fit1, fit2 )
	rslt <- cbind( rslt, c("", round( pchisq( rslt[2,4], rslt[2,3], lower.tail=FALSE ), 4 ) ) )
	rslt[,2] <- round( rslt[,2], 3 )
	rslt[,4] <- round( rslt[,4], 3 )
	rslt[1,3:4] <- c( "", "" )
	names( rslt )[5] <- "pValue"
	rslt
}

#####	H-L goodness of fit test
##
binary.gof <- function( fit, ngrp=10, print.table=TRUE ){
	y <- fit$y
	phat <- fitted( fit )
	fittedgrps <- cut( phat, quantile( phat, seq(0,1,by=1/ngrp) ), include.lowest=TRUE )
	n <- aggregate( y, list( fittedgrps ), FUN=length )[,2]
	Obs <- aggregate( y, list( fittedgrps ), FUN=sum )[,2]
	Exp <- aggregate( phat, list( fittedgrps ), FUN=sum )[,2]
	if( print.table==TRUE ){
		cat( "\nFitted Probability Table:\n\n" )
		rslt <- as.data.frame( cbind( 1:ngrp, n, Obs, Exp ) )
		names( rslt )[1] <- "group"
		print( rslt )
	}
	chisqstat <- sum( (Obs - Exp)^2 / ( Exp*(1-Exp/n) ) )
	df <- ngrp-2
	pVal <- pchisq( chisqstat, df, lower.tail=FALSE )
	cat( "\n Hosmer-Lemeshow GOF Test:\n\n" )
	cbind( chisqstat, df, pVal )
}

#####	Function to compute robust se for glms
##
robust.se.glm<-function(glm.obj){
	## 	Compute robust (sandwich) variance estimate
	if (is.matrix(glm.obj$x)) 
		xmat<-glm.obj$x
	else {
		mf<-model.frame(glm.obj)
		xmat<-model.matrix(terms(glm.obj),mf)		
	}
	umat <- residuals(glm.obj,"working")*glm.obj$weights*xmat
	modelv<-summary(glm.obj)$cov.unscaled
	robust.cov <- modelv%*%(t(umat)%*%umat)%*%modelv
	
	##	Format the model output with p-values and CIs
	s <- summary( glm.obj) 
	robust.se <- sqrt( diag( robust.cov )) 
	z <- glm.obj$coefficients/robust.se
	p <- 2*pnorm( -abs( z ) ) 
	ci95.lo <- glm.obj$coefficients - qnorm( .975 ) * robust.se
	ci95.hi <- glm.obj$coefficients + qnorm( .975 ) * robust.se
	rslt <- cbind( glm.obj$coefficients, robust.se, ci95.lo, ci95.hi, z, p ) 
	dimnames(rslt)[[2]] <- c( dimnames( s$coefficients )[[2]][1], "Robust SE", "ci95.lo", "ci95.hi", dimnames( s$coefficients )[[2]][3:4] ) 
	rslt 
	}





#########################################################################################
# Obtain exponentiated cofidence intervals for parameters
glmCI(m6.1)

glmCI(m3)


# Obtain exponentiated cofidence intervals for linear combinations of parameters
linContr.glm(c("sbp","sbp:age"), c(1,20), model=m6.1)

linContr.glm(c("sbp","sbp:age"), c(15,15*20), model=m6.1)

linContr.glm(c("MCAT","SexM:MCAT"), c(1,1), model=m3)


# Now for confidence interval for odds
linContr.glm(c("(Intercept)","sbp","age","sbp:age"), c(1,120,20,120*20),model=m6.1)

linContr.glm(c("(Intercept)","SexM","MCAT","SexM:MCAT"), c(1,1,30,30),model=m3)

linContr.glm(c("(Intercept)","SexM","MCAT","SexM:MCAT"), c(1,0,30,0),model=m3)

# Testing several coefficients at once
reduced = glm(chdfate~sbp, family=binomial(link="logit"),data= framingham)
full = glm(chdfate~sbp+age+sbp*age, family=binomial(link="logit"),data= framingham)

lrtest(reduced,full)

# Another example of the LR test
reduced1 = glm(Acceptance~GPA, family=binomial(link="logit"),data= mcat)
full1 = glm(Acceptance~GPA+Sex+MCAT, family=binomial(link="logit"),data= mcat)

lrtest(reduced1,full1)

# Conduct Hosmer Lemeshow test
binary.gof(full)

# Conduct test on a bigger model
full2 = glm(chdfate~sbp+age+female+bmi+sbp*female, family=binomial(link="logit"),data= framingham)

binary.gof(full2)

# ROC Curve for model 6.1 
roc.curve = roc(framingham$chdfate~fitted(m6.1))
plot(roc.curve)
roc.curve

roc.curve = roc(mcat$Acceptance~fitted(m3))
plot(roc.curve)
roc.curve

roc.curve = roc(house$ac~fitted(m8.1))
plot(roc.curve)
roc.curve

####################################################
# Homework code

# Read in MCAT data
mcat = read.table("./MedGPA.txt", fill=TRUE, header=TRUE)

mcat$male = mcat$Sex

mod1 =glm(Acceptance~male, data = mcat, family="binomial")

summary(mod1)


mod2 =glm(Acceptance~MCAT, data = mcat, family="binomial")

summary(mod2)


mod3 =glm(Acceptance~male+MCAT+MCAT*male, data = mcat, family="binomial")

summary(mod3)

linContr.glm(c("MCAT", "maleM:MCAT"), c(1,1), mod3)

lrtest(mod2, mod3)



# Import the Midwest house price data
midwest = read.table("./MidwestSales.txt", fill=TRUE, header=FALSE)

# This dataset does not have names, so we will add names to the variables
names(midwest)=c("id","price","sqft","bed","bath","ac","garage","pool","year","quality","style","lot","hwy")


mod = glm(ac~sqft+lot+pool, family="binomial", data=midwest)

summary(mod)

roc.curve = roc(midwest$ac~fitted(mod))
plot(roc.curve)
roc.curve

par(mfrow=c(1,2))
presids = residuals(mod, type="pearson")
muhat = fitted(mod)
plot(muhat, presids^2, xlab="Fitted expected counts", ylab="Pearson Residual Squared")
sfit = supsmu(muhat, presids^2)
lines(sfit$x[order(sfit$x)] , sfit$y[order(sfit$x)], col="red", lwd=2)

# Remove outlier
plot(muhat, presids^2, xlab="Fitted expected counts", ylab="Pearson Residual Squared", ylim=c(0,10))
sfit = supsmu(muhat, presids^2)
lines(sfit$x[order(sfit$x)] , sfit$y[order(sfit$x)], col="red", lwd=2)


summary(midwest[, c(3,8,12)])

midwest[which(hatvalues(mod) == max(hatvalues(mod))),]

linContr.glm( c("sqft" , "lot") , c(500,1500) , model=mod)


# Import nhanes data
nhanes = read.table( "./nhaneshw.txt", header=TRUE )

# Get summary of data based on male or female
lapply( split( nhanes, nhanes$male), summary )

# Create needed age group
nhanes$agegrp = cut( nhanes$age, breaks=c(0,30,40,50,60,71), right=FALSE )


fit1.full = glm( htn ~ factor(agegrp) + wt + male, family=binomial, data=nhanes )
glmCI( fit1.full )
fit1.red = glm( htn ~ wt + male, family=binomial, data=nhanes )
lrtest( fit1.red, fit1.full )


##	Test of interaction with gender
fit2 = glm( htn ~ factor(agegrp) + wt + male + factor(agegrp)*male, family=binomial, data=nhanes )
glmCI(fit2)
glmCI(fit2, transform=FALSE)
lrtest( fit1.full, fit2 )



##	Point estimates
aggregate( nhanes$wt, list(nhanes$male, nhanes$agegrp), mean )

linContr.glm( contr.names=c("(Intercept)", "factor(agegrp)[60,71)", "wt", "male", "factor(agegrp)[60,71):male"), contr.coef=c(1,1,85.543,1,1), model=fit2 )

linContr.glm( contr.names=c("(Intercept)", "factor(agegrp)[60,71)", "wt", "male", "factor(agegrp)[60,71):male"), contr.coef=c(1,1,85.543,1,1), model=fit2, transform=FALSE )

exp(-0.599)/(1+exp(-0.599))

exp(-0.802)/(1+exp(-0.802))

exp(-0.397)/(1+exp(-0.397))

linContr.glm( contr.names=c( "factor(agegrp)[40,50)", "factor(agegrp)[60,71)", "factor(agegrp)[40,50):male", "factor(agegrp)[60,71):male"), contr.coef=c(-1,1,-1,1), model=fit2 )

linContr.glm( contr.names=c( "factor(agegrp)[40,50)", "factor(agegrp)[60,71)"), contr.coef=c(-1,1), model=fit2 )


plot( hatvalues(fit2), cooks.distance(fit2), xlab="Leverage", ylab="Cooks Distance", xlim=c(0,.025) )
abline( h=qf(.5,length(fit2), fit2$df.residual), col="blue", lwd=2 )  
abline( v=2*mean(hatvalues(fit2) ), col="red", lwd=2 )
text( hatvalues(fit2)+.001, cooks.distance(fit2), nhanes$rtid )


################################################################################
# Poisson Regression section

# Read in crab data and fit some model
crab = read.csv("./crab.csv", header=TRUE)

crab.model = glm(Sa~W, family=poisson, data=crab)

crab.model2 = glm(Sa~W+Wt+W*Wt, family=poisson, data=crab)


# Read in ERtemp data and fit model
ERtemp = read.csv("./ERtemp.csv", header=TRUE)

ER.model = glm(Admissions~Temperature, family=poisson, data=ERtemp)

linContr.glm(c("(Intercept)", "Temperature"), c(1, 85), model = ER.model)

# Read in std data and fit some model
std = read.csv("./stdgrp.csv")

std.model = glm(n.reinfect~condom.always+offset(log(yrsfu)), family=poisson, data=std)

std.model2 = glm(n.reinfect~white+condom.always+edugrp+offset(log(yrsfu)), family=poisson, data=std)

# Run the next two lines for std.model3
std$edugrp = factor(std$edugrp, ordered=FALSE)
std$edugrp = relevel(std$edugrp, ref="[6,11.9]")

std.model3 = glm(n.reinfect~white+condom.always+relevel(edugrp, ref="[6,11.9]")+offset(log(yrsfu)), family=poisson, data=std)

lrtest(std.model,std.model2)

glmCI(std.model2)

# Read in Gail data, transform latitude, and fit model.
gail = read.table("./gail.txt", header=TRUE)

gail$latitude.ord = 1
gail$latitude.ord[ gail$latitude=="Middle  " ] = 2
gail$latitude.ord[ gail$latitude=="Southern" ] = 3

model.gail = glm(formula = inccases ~ factor(ageg) + latitude.ord, family = poisson(link=log), data = gail, offset = log(persyrs))

# Fit Gail model with quasipoisson
model.gail2 = glm(formula = inccases ~ factor(ageg) + latitude.ord, family = quasipoisson(link=log), data = gail, offset = log(persyrs))

# Obtain robust standard errors for model.gail
robust.se.glm(model.gail)


#################################

creditcard = read.table("./creditcard.txt")

lcases= log(creditcard[,2])
creditcard =cbind(creditcard,lcases)
colnames(creditcard)=c("income","cases","CrCards","lcases")

card.model=glm(CrCards~income+offset(lcases),family=poisson,data= creditcard)
summary(card.model)


beetles = read.table("./beetles.txt", header=TRUE)

#################################
#################################
# Multinomial

travel = read.csv("./travel.csv")
travel$Travel = rep( 1:4, 210 )
travel = travel[ travel$Mode==1, ]
travel = travel[ ,c("Travel", "Hinc", "Psize" ) ]
##
#####	
##
travel[1:10,]
dim( travel )

# Need the nnet library for multinomial regression 
library( nnet )

# Fit model for travel data
mfit = multinom( Travel ~ Hinc + Psize, data=travel )
mfit
summary( mfit )

# A function that will provide some more summary realted to multinomial
summ.mfit = function( model ){
	s = summary( model )
	for( i in 1:length(model$coef) ){
		cat( "\nLevel ", model$lev[i+1],  "vs. Level ", model$lev[1], "\n" )
		coef = s$coefficients[i,]
		rrr = exp( coef )
		se = s$standard.errors[i,]
		zStat = coef / se
		pVal = 2*pnorm( abs(zStat), lower.tail=FALSE )
		ci95.lo = exp( coef - qnorm(.975)*se )
		ci95.hi = exp( coef + qnorm(.975)*se )
		rslt = cbind( rrr, se, zStat, pVal, ci95.lo, ci95.hi )
		print( round( rslt, 3 ) )
	}
}

summ.mfit(mfit)

# Predict probabilities for Hinc=45 and Psize=3
newdata = data.frame(Hinc=45, Psize=3)
predict(mfit, type="probs" , newdata=newdata)

# Changing travel type to names (instead of numbers)
travel$Travel2[travel$Travel==1] = "Plane"
travel$Travel2[travel$Travel==2] = "Train"
travel$Travel2[travel$Travel==3] = "Bus"
travel$Travel2[travel$Travel==4] = "Car"

mfit2 = multinom( Travel2 ~ Hinc + Psize, data=travel )
mfit2
summary( mfit2 )


# Multinomial example for homework
abortion = read.table("./abortion.txt",col.names=c("year", "rel", "edu", "att", "count"))

# Fit model using weights
mfit.abort = multinom( att ~ edu+rel, data=abortion, weights=count )

summary(mfit.abort)

# Predict probabilties using multinomial distribution
newdata = data.frame(edu="Low", rel="Prot")
predict(mfit.abort, type="probs", newdata=newdata)

# Another predicition
newdata = data.frame(edu="High", rel="Prot")predict(mfit.abort, type="probs", newdata=newdata)



 