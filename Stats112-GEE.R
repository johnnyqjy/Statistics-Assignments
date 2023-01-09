##### Load Required Libraries #####
library(lattice) 
library(lme4)    
library(geepack) 


# Two functions needed: mymin and stacked2wide

mymin = function(x){ ifelse( sum( !is.na(x) )==0, NA, min(x, na.rm=TRUE ) ) }

stacked2wide = function( id, y, x, xvalues, half.window ){
 
# ARGUMENTS:
#            id          = id vector  (nobs x 1)
#            y           = response vector (nobs x 1)
#            x           = predictor (nobs x 1 )
#            xvalues     = vector of locations (m x 1)
#            half.window = max legal distance for |x - xvalue[j]|
#
# RETURN:
#           a matrix with y values for each subject entered 
#           as a row.  Each column, j, corresponds to the measurement
#           closest to xvalue[j], but not further than half.window
#           units from xvalue[j].  The variable x is typically time.
#
  nobs = length( id )
  nsubjects = length( table( id ) )
  m = length( xvalues )
  rmat = matrix( NA, nsubjects, m )
  nj = unlist( lapply( split( id, id ), length ) )
  for( j in 1:m ){
	legal = (x >= xvalues[j]-half.window)&(x < xvalues[j]+half.window)
	jtime = x + 0.01*rnorm(nobs)
	t0 = unlist( lapply( 
		split( abs(jtime - xvalues[j]), id ), mymin ) )
	tj = rep( t0, nj )
	keep = ( abs( jtime - xvalues[j] )==tj )&( legal )
	yj = rep( NA, nobs )
	yj[keep] = y[keep]
	yj = unlist( lapply( split( yj, id ), mymin ) )
	rmat[ , j ] = yj
  }
dimnames(rmat) = list( NULL, paste("X",xvalues) )
rmat
}



#
#### MADRAS schizophrenia symptom data: symptom = thought disorders
#### The variables are (in column order):
####     ID     = patient ID
####     Y      = symptom indicator
####     MONTH  = month since hospitalization
####     AGE    = age-at-onset (1= age<20 ; 0= age>=20)
####     GENDER = gender (1=female; 0=male)
#
#### Load in MADRAS data
#
madras = read.table( "./madras_data.txt" )[,1:5]
names( madras ) = c( "id", "y", "month", "age", "gender" )
##
#####
#####	Basic descriptive statistics
#####
#
#### Number of subjects with corresponding number of observations
table( table ( madras$id ) )
#
#### Drop individuals with less than 3 observations
drop.names = as.numeric( names( table( madras$id )[table( madras$id ) < 3 ] ) )
madras = madras[!is.element( madras$id, drop.names ),]
#
#### Number of subjects that are female/male and young/old
table( madras$age[madras$month == 0], madras$gender[madras$month == 0] )
table( madras$age[madras$month == 11], madras$gender[madras$month == 11] )
#
#### Plot of observed prevalence over time
#
prev.table = matrix( 0, 4, 12 )
for( i in 0:1 ){
	for( j in 0:1 ){
    temp.table = table( madras$y[madras$age == i & madras$gender == j], madras$month[madras$age == i & madras$gender == j] )
    prev.table[((2*i)+j+1),] = temp.table[2,] / apply( temp.table, 2, sum ) * 100
	}
}
dimnames( prev.table ) = list( c("age>=20 & male", "age>=20 & female", "age<20 & male", "age<20 & female"), paste("mon",0:11) )
round( prev.table, 0 )

plot( 0:11, prev.table[1,], xlab = "Months after Discharge", ylab = "Prevalence of Symptom", type = "n", ylim = range( prev.table ) * 1.1 )
lines( 0:11, prev.table[1,], col = 1, lwd = 3 )
lines( 0:11, prev.table[2,], col = 4, lwd = 3, lty = 5 )
lines( 0:11, prev.table[3,], col = 6, lwd = 3, lty = 3 )
lines( 0:11, prev.table[4,], col = 8, lwd = 3, lty = 4 )
legend( 6.0, 83, col = c(1,4,6,8), lwd = c(3,3,3,3), lty = c(1,5,3,4),
    	c("age>=20 & male", "age>=20 & female", "age<20 & male", "age<20 & female") )

#
#####
#####	Investigation of correlation structure
#####
#
stacked = stacked2wide( madras$id, madras$y, madras$month, 0:11, 0.5 )
cor.mat = or.mat = num.mat = matrix( 0, 12, 12 )
for( i in 1:11 ){
	for( j in i:12 ){
    temp.stacked = na.omit( stacked[,c(i,j)] )
    month.i = temp.stacked[,1]
    month.j = temp.stacked[,2]
    #
    cor.mat[i,j] = cor( month.i, month.j )
    table.ij = table( month.i, month.j )
    or.mat[i,j] = (table.ij[2,2] * table.ij[1,1]) / (table.ij[1,2] * table.ij[2,1]) 
    # or.mat[i,j] = 1 / or.mat[i,j]
    num.mat[i,j] = length( month.i )
    #
	}
}
dimnames( cor.mat ) = dimnames( or.mat ) = dimnames( num.mat ) = list( paste("mon", 0:11), paste("mon", 0:11) )
round( cor.mat, 2 )
round( or.mat, 2 )
num.mat
#
#####
#####	Model fitting with GEE
#####
#
#####	Main effects model with compound symmetry correlation structure
#
fit0 = geeglm( y ~ age + gender + month,id = id, data = madras,family = binomial(link="logit"),corstr="exchangeable" )
summary( fit0 )
#
#####	Scientific question using compound symmetry correlation structure
#	
fit1 = geeglm( y ~ age + gender + month + gender:month,id = id, data = madras, family = binomial(link="logit"), corstr="exchangeable" )
round( summary(fit1)$coef, 3 )

## To test the interaction term
anova(fit1,fit0)
