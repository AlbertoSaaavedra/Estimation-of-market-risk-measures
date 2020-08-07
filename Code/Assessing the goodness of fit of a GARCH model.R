# Assessing the goodness of fit of a GARCH model

# The function ajusta_garch has two purposes: 
#   i) to fit a GARCH model of the order indicated by the user; and 
#  ii) evaluate the quality of fit of said model.
# 
# The ajusta_garch function has the following parameters:
#    datos: The set of numerical observations of the time series to analyze.
#   fechas: The set of dates associated with the data.   
#   orden: two-input vector (p, q) to indicate the order of the GARCH process 
#          that we want to fit to the data.
# 
# The options that the user has when running the ajusta_garch function are: 
#   i) graph the estimated volatility (via the adjusted GARCH) of the data; 
#  ii) graph the variability of innovations; 
# iii) graph the ACF and PACF of the innovations (raw, squared and in absolute value); 
#  iv) obtain the p-values of the Ljung-Box test applied to innovations; 
#  iv) assign the fitted model to an R object

ajusta_garch<-function(datos, fechas, orden_garch=c(p,q)) {
  # Library to fit GARCH models to a data set
  library(fGarch)
  orden_garch<-as.numeric(orden_garch)
  datos<-as.numeric(datos)
  attr(datos,"times")<-fechas
  fit<-garchFit(formula=substitute(~garch(p,q),
                                   list(p=orden_garch[1],
                                        q=orden_garch[2])),
                data=datos,cond.dist=c("QMLE"),
                trace=F, inlude.mean=F) # GARCH fitting
  # Extracting GARCH innovations
  res<-residuals(fit, standardize = TRUE)
  attr(res,"times")<-fechas 
  op<-par() # Standard graph settings
  # numbering the choices given by the function
  choices<-c("Gr醘ica de volatilidad estimada de los datos", 
            "Gr醘ica de la variabilidad de las innovaciones",
            "ACF y PACF de las innovaciones (crudas, al 
             cuadrado y en valor absoluto)",
            "Prueba de Ljung-Box (innovaciones crudas, al 
             cuadrado y en valor absoluto)",
            "Asignar el modelo a un objeto") 
  tmenu<-paste("Opci髇:", choices) 
  pick<-1 
  lastcurve<-NULL
  while(pick > 0) {
    pick<-menu(tmenu, title = "\nSelecciona una opci髇 (o 0 
                                 para salir):")
    if(pick == 1) {
      #Estimando la volatilidad de los datos
      vol<-volatility(fit) 
      D<-matrix(c(datos,vol),ncol=2,nrow=length(datos))
      vol_completo_plot<-timeSeries(D,fechas_final)
      colnames(vol_completo_plot) <- c("Datos", "Volatilidad 
                                       estimada")
      par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,3.1,2.1),
          cex.axis=0.52) # Time series plot
      plot(vol_completo_plot, 
           main="Estimaci贸n GARCH de la Volatilidad de los 
                 datos", 
           plot.type=c("multiple"),
           col="slategrey", xlab="Fecha",
           at=pretty(fechas_final,8),
           yax.flip=T,format="%h-%Y") }
    if(pick == 2) {
      res_plot<-timeSeries(res,fechas)
      C<-matrix(c(res_plot,res_plot^2,abs(res_plot)),
                ncol=3,nrow=length(res_plot))
      res_plot_comp<-timeSeries(C,fechas)
      colnames(res_plot_comp) <- c("Innovaciones",
                                   "(Innovaciones)^2",
                                   "|Innovaciones|")
      plot(res_plot_comp, main="Variabilidad de las 
                                innovaciones",
           plot.type=c("multiple"),col="slategrey", 
           xlab="Fecha",at=pretty(fechas,8),yax.flip=T,
           format="%h-%Y") } 
    if(pick == 3) {
      par(mfrow=c(3,2),mar=c(4.08,4.1,3.3,1.68))
      acf(res, main="Innovaciones",sub="Autocorrelaci贸n")
      pacf(res, main="Innovaciones",
           sub="Autocorrelaci贸n parcial")
      acf(res^2, main="(Innovaciones)^2",
          sub="Autocorrelaci贸n")
      pacf(res^2, main="(Innovaciones)^2",
           sub="Autocorrelaci贸n parcial")
      acf(abs(res), main="|Innovaciones|", 
          sub="Autocorrelaci贸n")
      pacf(abs(res), main="|Innovaciones|", 
           sub="Autocorrelaci贸n parcial") }
    if(pick == 4) {
      print("P-valores de la Prueba de Ljung-Box:")
      print(matrix(c(Box.test(res,
                              lag=round(log(length(res)),0),
                              type="Ljung")$p.value,
                     Box.test(res^2,
                              lag=round(log(length(res)),0),
                              type="Ljung")$p.value,
                     Box.test(abs(res),
                              lag=round(log(length(res)),0),
                              type="Ljung")$p.value),
                   nrow=3,ncol=1,
                   dimnames = list(c("Innovaciones crudas",
                                     "Innovaciones al 
                                      cuadrado",
                                     "Innovaciones en valor 
                                      absoluto"),
                                   c("P-valor")))) }
    if(pick == 5)
      return(fit) }
  invisible(lastcurve) }