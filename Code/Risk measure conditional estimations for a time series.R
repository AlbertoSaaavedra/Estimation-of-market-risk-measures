# Risk measure conditional estimations for a financial time series

# The routine developed to estimate conditional risk measures for an analyzed 
# financial time series consists of seven functions.
# 
# The first six functions estimate the VaR and ES: there are three functions per 
# risk measure due to the different models considered for GARCH innovations.
# 
# The first six functions have the same parameters:
#   alfa: confidence level at which the risk measures are to be estimated (it
#         can be a numerical vector containing values between 0 and 1);
#   fit: R object that contains the GARCH fit made to the financial time 
#        series under study ; and
#   ind: dichotomous variable that indicates how volatility should estimated. A
#        value of 1 indicates that problem cases discussed in Chapter 5 should 
#        be considered, whereas a value of 0 ignores such problems and always
#        forecasts volatility for the next day based on the fitted GARCH model.

# Conditional VaR for Student-t innovations
VaR_garch_t<-function(alfa,fit,ind) {
  if(ind == 1) { # If GARCH is not stationary, we use an EWMA 
                 # approach to forecast volatility
    if((fit@fit$par['alpha1']+fit@fit$par['beta1'])>1 ||
        round(fit@fit$matcoef[2,4],4)>0.05) {
      vol<-sqrt((fit@fit$par['alpha1']*
                (fit@data[length(fit@data)]^2))+
                (volatility(fit)[length(volatility(fit))]^2)*
                (1-fit@fit$par['alpha1'])) }
    else
      vol<-predict(fit,1)$standardDeviation }
  if(ind == 0)
    vol<-predict(fit,1)$standardDeviation
  
  res<-residuals(fit, standardize = TRUE)
  par.t<-fitdistr(res,"t")$estimate
  nu<-par.t['df']
  qz<-qt(alfa,nu)
  VaR_alfa<-qz*vol
  return(VaR_alfa) }

# Conditional VaR for Normal innovations
VaR_garch_norm<-function(alfa,fit,ind) {
  if(ind == 1) { # If GARCH is not stationary, we use an EWMA 
                 # approach to forecast volatility
    if((fit@fit$par['alpha1']+fit@fit$par['beta1'])>1 ||
        round(fit@fit$matcoef[2,4],4)>0.05) {
      vol<-sqrt((fit@fit$par['alpha1']*
                (fit@data[length(fit@data)]^2))+
                (volatility(fit)[length(volatility(fit))]^2)*
                (1-fit@fit$par['alpha1'])) }
    else
      vol<-predict(fit,1)$standardDeviation }
  if(ind == 0)
    vol<-predict(fit,1)$standardDeviation
  
  res<-residuals(fit, standardize = TRUE)
  par.norm<-fitdistr(res,"normal")$estimate
  med<-par.norm['mean']
  sigma<-par.norm['sd']
  qz<-qnorm(alfa,med,sigma)
  VaR_alfa<-qz*vol
  return(VaR_alfa) }

# Conditional VaR for GPD innovations
VaR_garch_evt<-function(alfa,fit,ind) {
  if(ind == 1) {  # If GARCH is not stationary, we use an EWMA 
                  # approach to forecast volatility
    if((fit@fit$par['alpha1']+fit@fit$par['beta1'])>1 ||
        round(fit@fit$matcoef[2,4],4)>0.05) {
      vol<-sqrt((fit@fit$par['alpha1']*
                (fit@data[length(fit@data)]^2))+
                (volatility(fit)[length(volatility(fit))]^2)*
                (1-fit@fit$par['alpha1'])) }
    else
      vol<-predict(fit,1)$standardDeviation }
  if(ind == 0)
    vol<-predict(fit,1)$standardDeviation
  
  res<-residuals(fit, standardize = TRUE)
  gpd<-gpd(res,nextremes=100)
  qz<-riskmeasures(gpd,alfa)[2]
  VaR_alfa<-qz*vol
  return(VaR_alfa) }

# Conditional ES for Student-t innovations
ES_garch_t<-function(alfa,fit,ind) {
  if(ind == 1) {  # If GARCH is not stationary, we use an EWMA 
                  # approach to forecast volatility
    if((fit@fit$par['alpha1']+fit@fit$par['beta1'])>1 ||
        round(fit@fit$matcoef[2,4],4)>0.05) {
      vol<-sqrt((fit@fit$par['alpha1']*
                (fit@data[length(fit@data)]^2))+
                (volatility(fit)[length(volatility(fit))]^2)*
                (1-fit@fit$par['alpha1'])) }
    else
      vol<-predict(fit,1)$standardDeviation }
  if(ind == 0)
    vol<-predict(fit,1)$standardDeviation
  
  res<-residuals(fit, standardize = TRUE)
  par.t<-fitdistr(res,"t")$estimate
  nu<-par.t['df']
  ESz<-ES.t(alfa,par.t)
  ES_alfa<-ESz*vol
  return(ES_alfa) }

# Conditional ES for Normal innovations
ES_garch_norm<-function(alfa,fit,ind) {
  if(ind == 1) {  # If GARCH is not stationary, we use an EWMA 
                  # approach to forecast volatility
    if((fit@fit$par['alpha1']+fit@fit$par['beta1'])>1 ||
        round(fit@fit$matcoef[2,4],4)>0.05) {
      vol<-sqrt((fit@fit$par['alpha1']*
                (fit@data[length(fit@data)]^2))+
                (volatility(fit)[length(volatility(fit))]^2)*
                (1-fit@fit$par['alpha1'])) }
    else
      vol<-predict(fit,1)$standardDeviation }
  if(ind == 0)
    vol<-predict(fit,1)$standardDeviation
  
  res<-residuals(fit, standardize = TRUE)
  par.norm<-fitdistr(res,"normal")$estimate
  med<-par.norm['mean']
  sigma<-par.norm['sd']
  ESz<-ES.norm(alfa,med,sigma)
  ES_alfa<-vol*ESz
  return(ES_alfa) }

# Conditional ES for GPD innovations
ES_garch_evt<-function(alfa,fit,ind) {
  if(ind == 1) {  # If GARCH is not stationary, we use an EWMA 
                  # approach to forecast volatility
    if((fit@fit$par['alpha1']+fit@fit$par['beta1'])>1 ||
        round(fit@fit$matcoef[2,4],4)>0.05) {
      vol<-sqrt((fit@fit$par['alpha1']*
                (fit@data[length(fit@data)]^2))+
                (volatility(fit)[length(volatility(fit))]^2)*
                (1-fit@fit$par['alpha1'])) }
    else
      vol<-predict(fit,1)$standardDeviation }
  if(ind == 0)
    vol<-predict(fit,1)$standardDeviation
  
  res<-residuals(fit, standardize = TRUE)
  gpd<-gpd(res,nextremes=100)
  ESz<-riskmeasures(gpd,alfa)[3]
  ES_alfa<-vol*ESz
  return(ES_alfa) }

# Conditional risk measures - Synthesis

# In order to synthesize the results of the first six functions of our routine in 
# a single procedure, we created the function conditional_ risk_measures.
# 
# The conditional_measures_function has four parameters:
#     alfa: confidence level at which the risk measures are to be estimated (it 
#           can be a numerical vector containing values between 0 and 1);
#      fit: R object containing the GARCH fit made to the financial time series
#           under study; and
#      ind: dichotomous variable that indicates how volatility should  estimated. A
#           value of 1 indicates that problem cases discussed in Chapter 5 should be
#           considered, whereas a value of 0 ignores such problems and always forecasts
#           volatility for the next day based on the fitted GARCH model.
#   medida: it is a character string that indicates what risk measure we want to 
#           estimate. Valid values are: "VaR" and "ES".
# 
# The conditional_measures_function returns a matrix. The rows of this matrix correspond 
# to the models from which the VaR and ES are estimated. Its columns correspond to the
# confidence levels at which these risk measures are estimated.

medidas_riesgo_condicionales<-function(fit,medida,alfa,ind) {
  if(medida == "VaR") {
      return(matrix(c(VaR_garch_norm(alfa,fit,ind),
                      VaR_garch_t(alfa,fit,ind),
                      VaR_garch_evt(alfa,fit,ind)),
                    nrow=3, ncol=length(alfa),
                    dimnames = list(c("VaR, GARCH-Normal",
                                      "VaR, GARCH-t-Student",
                                      "VaR, GARCH-DPG"),
                                    as.character(alfa)),
                                    byrow=TRUE)) }
  if(medida == "ES") {
      return(matrix(c(ES_garch_norm(alfa,fit,ind),
                      ES_garch_t(alfa,fit,ind),
                      ES_garch_evt(alfa,fit,ind)),
                    nrow=3, ncol=length(alfa),
                    dimnames = list(c("ES, GARCH-Normal",
                                      "ES, GARCH-t-Student",
                                      "ES, GARCH-DPG"),
                                    as.character(alfa)),
                                    byrow=TRUE)) } }