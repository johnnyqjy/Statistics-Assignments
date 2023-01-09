### GLMMs and GEEs with Binary Data
##

##### Load Required Libraries #####
library(lattice) 
library(lme4)    
library(geepack) 

## Data on attitudes about abortion from Agresti (2013)

## Subjects were asked whether they supported legalizing
## abortion in each of three situations: 
## 1. If the family has a very low income and cannot afford any more children.
## 2. When the woman is not married and does not want to marry the man.
## 3. When the woman wants it for any reason.
## Gender was also measured.

## Read in data:
attitudes = read.table("./abortion2.txt", header=TRUE)

# Variables:
# gender = 1 if female; 0 if male
# response = 1 if support; 0 if do not support
# question = situation (1, 2, or 3)
# case = subject id

# Code situation as factor:
attitudes$question = factor(attitudes$question)
# Make situation 3 the reference situation:
attitudes$question = relevel(attitudes$question, ref="3")
# Note that we can leave gender and response as numeric
# since they are already coded as 0 or 1.

# Plot of proportion of support in each covariate pattern:
means = tapply(attitudes$response, list(attitudes$question,attitudes$gender), mean)
plot(c(3,1,2),means[1:3,1],type="p",col="red",pch=16,
	xlab="Situation",ylab="Proportion Support",
	main="Proportion Support by Situation and Gender",
	cex=2,ylim=c(0.45,.52))
points(c(3,1,2),means[1:3,2],col="darkblue",cex=2)	
legend("bottomleft", c("Male","Female"), col=c("red","darkblue"), pch=c(16,1))

# OR barplot:
barplot(t(means), beside=TRUE, col=c("Red","Blue"),
	ylim=c(0,0.6), legend = c("Male","Female"),
	xlab="Situation",ylab="Proportion Support",
	main="Proportion Support by Situation and Gender")

# Contingency table summary of data by situation:
my.tab = xtabs(~response+gender+question, data=attitudes)
my.tab
# Visualize three-way contingency table:
library(vcd)
pairs(my.tab, diag_panel = pairs_barplot,  labeling=TRUE)

# Binary data on three occasions -->
# only 8 possible values of response vector.

# Convert data to wide form:
att.wide = reshape(attitudes, v.names="response", timevar="question", idvar="case", direction="wide")
dim(att.wide) # 1850 subjects questioned
# How many supported legalizing abortion in all three situations,
# i.e., responded (1,1,1)?
# Create sequence or responses variable:
# (this code takes a few seconds... for loops are slow)
resp.seq = NULL
for(i in 1:1850){
	x = att.wide[i,3:5]
	if(x[1]==1 & x[2]==1 & x[3]==1){resp.seq[i]="111"}
	else if(x[1]==1 & x[2]==1 & x[3]==0){resp.seq[i]="110"}
	else if(x[1]==1 & x[2]==0 & x[3]==1){resp.seq[i]="101"}
	else if(x[1]==0 & x[2]==1 & x[3]==1){resp.seq[i]="011"}
	else if(x[1]==1 & x[2]==0 & x[3]==0){resp.seq[i]="100"}
	else if(x[1]==0 & x[2]==1 & x[3]==0){resp.seq[i]="010"}
	else if(x[1]==0 & x[2]==0 & x[3]==1){resp.seq[i]="001"}
	else if(x[1]==0 & x[2]==0 & x[3]==0){resp.seq[i]="000"}
}
att.wide$resp.seq = resp.seq
xtabs(~gender+resp.seq, data=att.wide)
# Note that the majority of the sample was either 000 or 111.
# --> high association between repeated measurements.
# Sample correlations among repeated measurements:
cor(att.wide[,3:5])


### Fit GLMM by ML ###
# No interaction between gender and question -->
# gender effect assumed to be identical for each situation.
mod1 = glmer(response ~ gender+question+(1|case), family=binomial, data=attitudes, nAGQ = 5)
# Larger value of "nAGQ" =  greater accuracy in approximation
# --> longer computational time	
# Can only go up to 25.
# Default is 1 (Laplace approximation) - gives very different results!
summary(mod1)
summary(update(mod1, nAGQ = 1))

mod2 = glmer(response ~ gender*question+(1|case), family=binomial,
	data=attitudes, nAGQ = 5)
# Test for interaction effect between gender and question:	
anova(mod2)
anova(mod1,mod2)

# Estimated random effects same for all individuals with same response sequence
# and same gender, e.g.,
ranef(mod1)$case[resp.seq=="000",1]
att.wide$gender[resp.seq=="000"]


### Fit Marginal Model by GEE ##
# Remember to sort data by case and then by situation:
att.sort = attitudes[order(attitudes$case,attitudes$question),]
mod2 = geeglm(response ~ gender+question, family=binomial(link="logit"),
	corstr="exchangeable", id=case, data=attitudes)
mod2.int = geeglm(response ~ gender*question, family=binomial(link="logit"),
	corstr="exchangeable", id=case, data=attitudes)
anova(mod2,mod2.int)	
summary(mod2)	
