# Value-at-Risk backtesting

# The last part of our application consists of evaluating our estimates of conditional risk 
# measures, via backtests. To carry out this procedure, two R routines were created, one 
# per risk measure to be evaluated.
# 
# The first routine performs the backtesting of VaR backtest and consists of a single 
# function called VaR_backtest_garch. This function has four parameters:
#    datos: the set of numerical observations of the time series on which the backtest
#           will be applied;
#   fechas: the set of dates associated with the datos cariable; 
#        n: the length of the time window on which we want to base the backtest; and
#      ind: dichotomous variable that indicates how volatility should  estimated. A
#           value of 1 indicates that problem cases discussed in Chapter 5 should be
#           considered, whereas a value of 0 ignores such problems and always forecasts
#           volatility for the next day based on the fitted GARCH model.
#
# The VaR_backtest_garch function uses the functions of the "Risk measure conditional 
# estimations for a time series" file to automatically estimate the VaR, from the 
# GARCH-Normal, GARCH-Student-t and GARCH-GPD models. 

# Based on the violations to these quantiles, the routine generates four vectors of 
# zeros and ones per model: the sample versions of the indicators of VaR violations.
# The routine then uses the binom.test function of R to apply a binomial test to
# these indicators
# 
# The function then generates five matrices summarizing the information from the 
# performed binomial tests: four containing the number of historical violations 
# and the p-values of the binomial tests, and one with information on the expected
# number of violations, based on the backtest lenght and the level at which the 
# quantiles were estimated.
# 
# Finally, the function generates a list containing the 5 previous matrices. This
# list is the object that the VaR_backtest_garch function returns.

VaR_backtest_garch<-function(datos,fechas,n,ind) {
  # VaR excess indicators for a GARCH model with Gaussian innovations
  ind.norm.95<-c()
  ind.norm.99<-c()
  ind.norm.995<-c()
  ind.norm.999<-c()
  # VaR excess indicators for a GARCH model with Student-t innovations
  ind.t.95<-c()
  ind.t.99<-c()
  ind.t.995<-c()
  ind.t.999<-c()
  # VaR excess indicators for a GARCH model with GPD innovations
  ind.evt.95<-c()
  ind.evt.99<-c()
  ind.evt.995<-c()
  ind.evt.999<-c()
  # Backtest length
  m2<-(length(datos)-n)
  # For loop fot backtesting through different days
  for(i in 1:(m2+1)) {
    datosfit<-datos[i:((n-1)+i)] # Selecting data to use on the i-th fit
    attr(datosfit,"times")<-fechas[i:((n-1)+i)] # Selecting dates to use on the i-th fit
    fit.qmle<-garchFit(formula=~garch(1,1), 
                       data=datosfit,cond.dist=c("QMLE"),
                       trace=F, inlude.mean=F, 
                       hessian = "ropt") # Making the i-th GARCH(1,1) fit
    # Detecting VaR violations
    # Gaussian innovations
    if((VaR_garch_norm(0.95,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.norm.95[i]<-1
    else
      ind.norm.95[i]<-0
    if((VaR_garch_norm(0.99,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.norm.99[i]<-1
    else
      ind.norm.99[i]<-0
    if((VaR_garch_norm(0.995,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.norm.995[i]<-1
    else
      ind.norm.995[i]<-0
    if((VaR_garch_norm(0.999,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.norm.999[i]<-1
    else
      ind.norm.999[i]<-0
    # Student-t innovations
    if((VaR_garch_t(0.95,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.t.95[i]<-1
    else
      ind.t.95[i]<-0
    if((VaR_garch_t(0.99,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.t.99[i]<-1
    else
      ind.t.99[i]<-0
    if((VaR_garch_t(0.995,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.t.995[i]<-1
    else
      ind.t.995[i]<-0
    if((VaR_garch_t(0.999,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.t.999[i]<-1
    else
      ind.t.999[i]<-0
    # GPD innovations
    if((VaR_garch_evt(0.95,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.evt.95[i]<-1
    else
      ind.evt.95[i]<-0
    if((VaR_garch_evt(0.99,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.evt.99[i]<-1
    else
      ind.evt.99[i]<-0
    if((VaR_garch_evt(0.995,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.evt.995[i]<-1
    else
      ind.evt.995[i]<-0
    if((VaR_garch_evt(0.999,fit.qmle,ind) < datos[(n+i)]) &
        ((n+i) <= length(datos)))
      ind.evt.999[i]<-1
    else
      ind.evt.999[i]<-0 }
  n<-m2
  # Binomial test p-values
  # Gaussian innovations
  pvalue.norm.95<-binom.test(sum(ind.norm.95),n,
                             p=0.05)$p.value
  pvalue.norm.99<-binom.test(sum(ind.norm.99),n,
                             p=0.01)$p.value
  pvalue.norm.995<-binom.test(sum(ind.norm.995),n,
                              p=0.005)$p.value
  pvalue.norm.999<-binom.test(sum(ind.norm.999),n,
                              p=0.001)$p.value
  # Student-t innovations
  pvalue.t.95<-binom.test(sum(ind.t.95),n,
                          p=0.05)$p.value
  pvalue.t.99<-binom.test(sum(ind.t.99),n,
                          p=0.01)$p.value
  pvalue.t.995<-binom.test(sum(ind.t.995),n,
                           p=0.005)$p.value
  pvalue.t.999<-binom.test(sum(ind.t.999),n,
                           p=0.001)$p.value
  # GPD innovations
  pvalue.evt.95<-binom.test(sum(ind.evt.95),n,
                            p=0.05)$p.value
  pvalue.evt.99<-binom.test(sum(ind.evt.99),n,
                            p=0.01)$p.value
  pvalue.evt.995<-binom.test(sum(ind.evt.995),n,
                            p=0.005)$p.value
  pvalue.evt.999<-binom.test(sum(ind.evt.999),n,
                             p=0.001)$p.value
  # Building the VaR 0.95 matrix
  viol.hist.95<-c(cont.norm.95,cont.t.95,cont.evt.95)
  pvalue.95<-c(pvalue.norm.95,pvalue.t.95,pvalue.evt.95)
  VAR.95<-matrix(c(viol.hist.95,pvalue.95),nrow=3,
                 ncol=2,byrow=FALSE,
                 dimnames = list(c("Modelo GARCH-Normal",
                                   "Modelo GARCH-t-Student",
                                   "Modelo GARCH-DPG"),
                                 c("Violaciones Históricas al 
                                    VaR 0.95",
                                   "P-Value Test Binomial")))
  # Building the VaR 0.99 matrix
  viol.hist.99<-c(cont.norm.99,cont.t.99,cont.evt.99)
  pvalue.99<-c(pvalue.norm.99,pvalue.t.99,pvalue.evt.99)
  VAR.99<-matrix(c(viol.hist.99,pvalue.99),nrow=3,
                 ncol=2,byrow=FALSE,
                 dimnames = list(c("Modelo GARCH-Normal",
                                   "Modelo GARCH-t-Student",
                                   "Modelo GARCH-DPG"),
                                 c("Violaciones Históricas al 
                                    VaR 0.99",
                                   "P-Value Test Binomial")))
  # Building the VaR 0.995 matrix
  viol.hist.995<-c(cont.norm.995,cont.t.995,cont.evt.995)
  pvalue.995<-c(pvalue.norm.995,pvalue.t.995,pvalue.evt.995)
  VAR.995<-matrix(c(viol.hist.995,pvalue.995),nrow=3,
                  ncol=2,byrow=FALSE,
                  dimnames = list(c("Modelo GARCH-Normal",
                                    "Modelo GARCH-t-Student",
                                    "Modelo GARCH-DPG"),
                                  c("Violaciones Históricas 
                                     al VaR 0.995",
                                   "P-Value Test Binomial")))
  # Building the VaR 0.999 matrix
  viol.hist.999<-c(cont.norm.999,cont.t.999,cont.evt.999)
  pvalue.999<-c(pvalue.norm.999,pvalue.t.999,pvalue.evt.999)
  VAR.999<-matrix(c(viol.hist.999,pvalue.999),nrow=3,
                  ncol=2,byrow=FALSE,
                  dimnames = list(c("Modelo GARCH-Normal",
                                    "Modelo GARCH-t-Student",
                                    "Modelo GARCH-DPG"),
                                  c("Violaciones Históricas 
                                     al VaR 0.999",
                                   "P-Value Test Binomial")))
  # Computing the number of expected VaR violations, based on the backtest lenght
  nombres<-c("0.95","0.99","0.995","0.999")
  viol.esp1<-c(n*(1-0.95),n*(1-0.99),n*(1-0.995),n*(1-0.999))
  viol.esp<-matrix(viol.esp1,nrow=1,ncol=4,byrow=T,
                  dimnames = list(c("Violaciones esperadas"),
                                  c("0.95","0.99","0.995",
                                    "0.999")))
  # Forming the list that our function will return
  Backtest_VaR<-list(VaR.95=VAR.95,VaR.99=VAR.99,
                    VaR.995=VAR.995,VaR.999=VAR.999,
                    Viol.Esp=viol.esp)
  return(Backtest_VaR) }