library(yarrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library (lmtest)
library(orcutt)
library(EnvStats)
library(nortest)
library(car)

testing.linearity <- function(x,y,...){
  lm.model <- lm (y~x)
  par(mfrow=c(2,2))
  
  plot(y,
       main = 'plot of scores',
       ylab = 'scores')
  text.lin1 <- paste('Do the scores follow any weird patterns?')
  mtext(text.lin1)
  abline (h=mean(y), col='red')
  text.lin.gral <- 'TESTING FOR LINEARITY'
  mtext(text.lin.gral, line = 3, side = 3, at= 175, col = 'darkblue')
  
  plot(residuals(lm.model),
       ylab = 'residuals',
       main = 'plot of residuals') 
  text.lin5 <- paste ('Are the points distributed linearly? Mean of residuals ~0? m=', round(mean(residuals(lm.model)),2))
  mtext(text.lin5)
  abline(h= mean(residuals(lm.model)), col='red')
  library (lmtest)
  lrt <- lrtest (lm.model)
  if (lrt$`Pr(>Chisq)`[2]>0.05) {col.lin <- 'red'}
  if (lrt$`Pr(>Chisq)`[2]<=0.05) {col.lin <- 'black'}
  text.lin6 <- paste('Likelihood ratio test p=', round(lrt$`Pr(>Chisq)`[2],3))
  mtext(text.lin6, line = 3, side = 3, col= col.lin)
  
  sc.fv <- lm (y~lm.model$fitted.values)
  plot(y~lm.model$fitted.values,
       ylab = 'scores',
       xlab= 'fitted values',
       main = 'scores vs fitted values')
  text.lin2 <- paste('Do the points follow a diagonal line? slope=', round(sc.fv$coefficient[2],2))
  mtext (text.lin2)
  abline(lm(y~lm.model$fitted.values), col='red')
  
  res.fv <- lm(lm.model$residuals~ lm.model$fitted.values)
  plot(lm.model$residuals~ lm.model$fitted.values,
       ylab = 'residuals',
       xlab= 'fitted values',
       main = 'residuals vs fitted values')
  text.lin3 <- paste ('Is the regression line = 0? slope=', round(res.fv$coefficients[2],5))
  mtext (text.lin3)
  abline (lm(lm.model$residuals~lm.model$fitted.values), col='red')
  text.lin7 <- paste('Your model gives signs of non-linearity! Consider logging the variables, \napplying a polynomial function or adding new variables!')
  if (lrt$`Pr(>Chisq)`[2]>0.05) {mtext(text.lin7, side = 3, line = 2.75, col = 'red')}
}

testing.homoskedasticity <- function(x,y) {
  lm.model <- lm (y~x)
  par(mfrow=c(2,2))
  
  res.fv <- lm(lm.model$residuals~ lm.model$fitted.values)
  plot(lm.model$residuals~ lm.model$fitted.values,
       ylab = 'residuals',
       xlab= 'fitted values',
       main = 'residuals vs fitted values')
  text.lin3 <- paste ('Do points grow or decrease? Is the regression line = 0? slope=', round(res.fv$coefficients[2],5))
  text.gral <- 'TESTING FOR HOMOSKEDASTICITY'
  mtext (text.lin3)
  abline (lm(lm.model$residuals~lm.model$fitted.values), col='red')
 
  gold.q.t <- gqtest (lm.model)
  breusch.p.t <- bptest (lm.model) 
  if(gold.q.t$p.value>0.05) {col.homo1 <- 'red'}
  if(gold.q.t$p.value<=0.05) {col.homo1 <- 'black'}
  if (breusch.p.t$p.value<=0.05) {col.homo2 <- 'black'}
  if (breusch.p.t$p.value>0.05) {col.homo2 <- 'red'}
  text.homo1 <- paste('Goldfeld-Quandt test p=', round(gold.q.t$p.value, 4))
  text.homo2 <- paste ('Breusch-Pagan test p=', round(breusch.p.t$p.value, 4))
  mtext(text.homo1, side = 1, line = 8, col = col.homo1)
  mtext(text.homo2, side = 1, line = 10, col = col.homo2)
  
  plot(residuals(lm.model),
       ylab = 'residuals',
       main = 'plot of residuals') 
  text.lin5 <- paste ('Are the errors spread randomly across the regression line? Mean of residuals ~0? m=', round(mean(residuals(lm.model)),2))
  mtext(text.lin5)
  mtext(text.gral, line = 3, side = 3, at=-20, col = 'darkblue')
  abline(h= mean(residuals(lm.model)), col='red')
  
  whites.errors <- coeftest(lm.model)
  text.homo2 <- paste('Your model gives signs of heteroskedasticity! \nFollowing White s standard errors (free of heteroskedasticity), \nlook at the print for your corrected regression results. \nYou could also consider running the regression on weighted least squares \nor just log all the variables')
  if (gold.q.t$p.value>0.05 | breusch.p.t$p.value >0.05) {mtext(text.homo2, side=1, line = 11, col = 'red')}
    if (gold.q.t$p.value>0.05 | breusch.p.t$p.value >0.05) {print(whites.errors)}
}

testing.autocorrelation <- function(x,y) {
  par(mfrow=c(2,2))
  lm.model <- lm(y~x)
  durbin.w.test <- dwtest(lm.model, alternative = 'two.sided')
  breusch.g.test <- bgtest(lm.model, order = 5)
  if (durbin.w.test$p.value <0.05) {col.autoc.d <- 'red'}
  if (durbin.w.test$p.value >=0.05) {col.autoc.d <- 'black'}
  if (breusch.g.test$p.value <0.05) {col.autoc.b <- 'red'}
  if (breusch.g.test$p.value >=0.05) {col.autoc.b <- 'black'}
  
  plot(lm.model$residuals,
       ylab = 'residuals')
  mtext('Does this follow any sneaky pattern?')
  text.gral <- 'TESTING FOR AUTOCORRELATION'
  mtext(text.gral, line = 3, side = 3, at=150, col = 'darkblue')
  
  text.autoc2 <- paste ('Durbin-Watson test for autocorrelation of consecutive terms p=', round (durbin.w.test$p.value,4))
  text.autoc3 <- paste ('Breusch-Godfrey test for autocorrelation of distant terms p=', round (breusch.g.test$p.value,4))
  mtext(text.autoc3, side = 1, line = 8, col = col.autoc.b)
  mtext (text.autoc2, side = 1, line = 10, col = col.autoc.d)
  
  plot(lm.model$residuals~ df$row.n,
       ylab = 'residuals',
       xlab = 'number of row')
  mtext('The residuals should be randomly and symmetrically distributed around zero')
  
  if (durbin.w.test$p.value <0.05 | breusch.g.test$p.value <0.05) {text.autoc1 <- 'Careful! There are signs of autocorrelation among \nconsecutive error terms. \nLook at the print for a generalised difference equation \nthrough Cochrane-Orcutt procedure'}
  if (durbin.w.test$p.value <0.05 | breusch.g.test$p.value <0.05) {mtext (text.autoc1, side = 1, line = 10, col = 'red')}  
  if (durbin.w.test$p.value <0.05 | breusch.g.test$p.value <0.05) {print(cochrane.orcutt(lm.model))}
}

testing.normality <- function(x,y) {
  par(mfrow=c(2,2))
  lm.model <- lm(y ~x)
  
  qqPlot (y,
          main = 'Normal Q-Q plot')
  mtext('The errors should be normally distributed along the diagonal', side = 3)
  text.gral <- 'TESTING FOR ERROR NORMALITY'
  
  
  kolmogrov <- ks.test((lm.model$residuals), rnorm(20000, mean= mean(lm.model$residuals), sd= sd(lm.model$residuals)))
  shapiro <- shapiro.test(lm.model$residuals)
  anderson <- ad.test (lm.model$residuals)
  
  # Define and add top-text
  top.text.k <- paste('Kolmogorov-Smirnov test: p-value= ', round(kolmogrov$p.value, 4))
  top.text.s <- paste ('Shapiro-Wilk test: p-value=', round(shapiro$p.value, 4))
  top.text.a <- paste ('Anderson-Darling test: p-value=', round(anderson$p.value, 4))
  if (kolmogrov$p.value <0.05) {col.k <- 'red'}
  if (kolmogrov$p.value >=0.05) {col.k <- 'black'}
  if (shapiro$p.value <0.05) {col.s <- 'red'}
  if (shapiro$p.value >=0.05) {col.s <- 'black'}
  if (anderson$p.value <0.05) {col.a <- 'red'}
  if (anderson$p.value >=0.05) {col.a <- 'black'}
  mtext(top.text.k, side = 1, line=8, col = col.k)
  mtext(top.text.s, side = 1, line = 9, col = col.s)
  mtext(top.text.a, side = 1, line = 10, col = col.a)

  hist(resid(lm.model), breaks = 10,
       main = 'Histogram of error terms',
       xlab = 'Error terms')
  mtext('Error terms should be normaly distributed', side = 3)
  mtext(text.gral, line = 3, side = 3, at=-2, col = 'darkblue')
  if (kolmogrov$p.value <0.05 | shapiro$p.value <0.05 | anderson$p.value< 0.05) {text.norm1 <- 'Careful! There are signs of residuals non-normality. \nThe distribution of the variables themselves may not be normal. \nYou may want to consider logging variables.\nYou may check or remove the outliers identified in the print!'}
  if (kolmogrov$p.value <0.05 | shapiro$p.value <0.05 | anderson$p.value< 0.05) {mtext (text.norm1, side = 1, line = 10, col = 'red')}  
  
  outlier.def <- 2
  is.outlier <- lm.model$residuals> mean(lm.model$residuals) + outlier.def * sd (lm.model$residuals) |
      lm.model$residuals < mean(lm.model$residuals) - outlier.def * sd (lm.model$residuals)
  no.outlier <- lm.model$residuals [is.outlier== FALSE]
  print (df$row.n[is.outlier==TRUE])
}
