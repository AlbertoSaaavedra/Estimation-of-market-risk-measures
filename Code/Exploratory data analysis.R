# Exploratory data analysis

# The analysis_exp function allows us to obtain the descriptive graphs
# and measurements used in our exploratory analyzes. The function 
# displays a menu for the user to choose which element of the 
# exploratory analysis they want to obtain.
# 
# The function parameters are:
#        x: the vector of (numerical) observations to which we wish to
#           carry out the exploratory analysis.
#   umbral: the threshold from which the graph of the mean excess
#           function will be drawn.
# 
# The options that the user has when running the analysis_exp function are: 
#    i) the histogram; 
#   ii) the box and arm diagram; 
#  iii) the QQ chart against a Normal reference distribution; 
#   iv) the graph of the average sample function of excess; and 
#    v) a matrix containing various measures of central tendency.
# By default, the analysis_exp function works only with the right tail 
# of the distribution of the data to be scanned, however this can be
# easily modified (line 23) to analyze the left tail or both queues of
# the observations.

analisis_explor<-function(x, umbral) {
  library(evir) # Library to  plot the mean excess function
  library(e1071) # Library with the kurtosis function 
  x<-as.numeric(x) # Transforming data to numeric values
  # Numbering the options that the function gives to the user
  choices<-c("Histograma","Boxplot","QQ-Plot (Dist. Normal)", 
               "Función muestral media de exceso",
               "Medidas de tendencia central")
  tmenu<-paste("Opción:", choices) # To deploy the menu
  pick <- 1 # Variable to store user's selection
  lastcurve <- NULL
  while (pick > 0) {
    pick<-menu(tmenu,title="\nSelecciona una opci�n (o 0 para salir):")
    if(pick == 1) {
      hist(x, ylab = "Frecuencia", xlab = "Datos", ...) }
    if(pick == 2) 
      boxplot(x, ...)
    if(pick == 3){
      qqnorm(x)
      qqline(x, col="red") }
    if(pick == 4) {
      meplot(x[x>0], umbral)
      abline(v=umbral,col="red",lty="dashed") }
    if(pick == 5){
      print(matrix(c(min(x),max(x),median(x),mean(x), 
                     var(x),kurtosis(x)),
                   nrow=1, ncol=6,
                   dimnames = list(c("Medidas"),c("Mínimo",
                   "Máximo","Mediana","Media","Varianza",
                   "Ex. Curtosis")))) } }
  invisible(lastcurve) }
