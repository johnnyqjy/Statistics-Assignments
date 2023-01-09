### GLMMs and GEEs with FLW Epilepsy Data
##

##### Load Required Libraries #####
library(lattice) # Functions: xyplot
library(lme4)    # Functions: glmer
library(geepack) # Functions: geeglm

##
#### Read in data set into R:
##
epilepsy = read.table("./epilepsy.txt", header=TRUE)

# Convert to long form:
epi.long = reshape(epilepsy, idvar="ID", varying=list(4:8), v.names="Count", timevar="Time", times = seq(0,8,by=2), direction="long")

##### Exploratory Data Analysis #####

## Number of seizures by group by time period:
boxplot(Count ~ trt+Time, at=c(1,2,4,5,7,8,10,11,13,14), data=epi.long, col=c("red","blue"), ylab="Number of Seizures", xlab="", names=rep(c("Plac","Prog"),5), main="Number of Seizures by Group and Week")
mtext( c("Baseline",paste("Week",c(2,4,6,8))), side=1, at=c(1.5,4.5,7.5,10.5,13.5), line=3 )

# or
bwplot(Count ~ trt | Time, data = epi.long)

# Calculate rates per week:
epi.long$Rate = epi.long$Count
epi.long$Rate[epi.long$Time==0] = epi.long$Count[epi.long$Time==0]/8
epi.long$Rate[epi.long$Time!=0] = epi.long$Count[epi.long$Time!=0]/2

# Plot of mean rates:
means = tapply(epi.long$Rate,list(epi.long$Time,epi.long$trt),mean)
matplot(matrix(c(0,2,4,6,8)),means,col=c(1,1),lty=c(3,1),type="o",pch=c(1,16),xlab="Time (weeks)",ylab="Mean rate of seizures (per week)",ylim=c(2.5,5.0),main="Figure 1.2: Mean Rate of Seizures by Treatment Group")
legend(3.5,3.0, c("Placebo","Progabide"), lty=c(3,1))

## Plots of individual counts:
xyplot(Count ~ Time | trt, data = epi.long)

xyplot(Count ~ Time | trt, groups=ID, type="l", data = epi.long)

matplot(matrix(c(0,2,4,6,8)),t(epilepsy[,4:8]),
	col=(as.numeric(epilepsy$trt=="Placebo")+2), type="l", ylim = c(0,180),
	xlab="Week",ylab="No. of Seizures ", main="Response Profiles by ID and Treatment")
legend(5,150,c("Placebo","Progabide"),col=c(2,3),lty=c(1,1))

 

##### Fitting GLMMs #####
# Since the number of weeks at each count differs
# (8 weeks for baseline by 2 weeks afterwards)
# we need to include an "offset" --> Model mean rate per week
epi.long$PostBase = as.numeric(epi.long$Time!=0)  # Indicator for post-baseline
epi.long$Weeks = 8*(epi.long$PostBase==0) + 2*(epi.long$PostBase==1)

# With glmer (and lmer) function, random effects specified in parentheses:
?lmer
?glmer


mod1 = glmer(Count ~ trt*PostBase + (PostBase | ID), offset=log(Weeks),
	family=poisson, data=epi.long)
summary(mod1)


mod2 = glmer(Count ~ trt+PostBase + (PostBase | ID), offset=log(Weeks),
	family=poisson, data=epi.long)
summary(mod2)

# Estimated random effects for model 1
ranef(mod1)
# Subject-specific coefficients (beta_k + b_ki)
coef(mod1)

# Treating time as quantitative:
mod3 = glmer(Count ~ trt*Time + (Time | ID), offset=log(Weeks),
	family=poisson, data=epi.long)
summary(mod3)


mod4 = glmer(Count ~ trt+Time + (Time | ID), offset=log(Weeks),
	family=poisson, data=epi.long)
summary(mod4)

 


##### Fitting Marginal Models #####
?geeglm

### IMPORTANT: geeglm needs data sorted first by ID (cluster),
### and then by time within ID
epi.sort = epi.long[order(epi.long$ID, epi.long$Time),]
# Check:
epi.sort

mod.gee = geeglm(Count ~ trt*PostBase, offset=log(Weeks),
				family=poisson(link="log"), id=ID, corstr="exchangeable",
				data=epi.sort)
# corstr="exchangeable" --> compound symmetry covariance structure.
# family=poisson --> Poisson variance function (not distribution)
summary(mod.gee)

mod.gee2 = geeglm(Count ~ trt*PostBase, offset=log(Weeks),
				family=poisson, id=ID, corstr="ar1", data=epi.sort)
summary(mod.gee2)
# Note similar coefficient estimates and Wald tests
# GEE std. errs robust to covariance structure assumption
