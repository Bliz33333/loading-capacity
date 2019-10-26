
#Required Parameters: Must be given values
ALIQUOT_RATIO <- 4        #metal:dendrimer ratio of each addition
my_path <- "~/R Workspace/Data Tool Project/g6nh2-15um-nickel-ratio_4-trial3"
                          #the location of the files to analyze


#Optional Parameters: Can leave them blank. Leave the brackets
MAX_SEARCH_OFFSET <- {}   #normally the program looks at the last spec to determine the peaks in the data. by setting this to some number x,
                            #you can ignore the last x specs when searching for peaks. For use in data where the last few specs are messy
                            #due to high absobances. Default is 0
HEADER_OFFSET <- {}      #The ratio of metal:dendrimer before the first measured aliquot. Is zero unless you preloaded metal
                            #default is 0
Y_MAX <- {}               #default is the max absorbance in the range of wavelength = 225 to 900
X_MAX <- {250}            #the end point of the wavelengths you want the program to consider
                            #default is 900
X_MIN <- {}               #the start point of the wavelengths you want the program to consider
                            #default is 225
lambda <- {235}           #the wavlength that will be graphed to determine loading capacity
                            #default is the absolute max of the searched spectrum in the wavelength range of X_MIN to X_MAX

  #####################################################################################################

  
  
  #read in required packages
  require(readxl)
  require(tidyverse)
  
  #set the working directory from which the files will be read from
  setwd(my_path)
  
  #create a list of the files from your target directory
  file_list <- list.files(path=my_path)
  
  file_list[1]
  
  
  if(is.null(MAX_SEARCH_OFFSET))
  {
    MAX_SEARCH_OFFSET <- 0
  }
  
  if(is.null(HEADER_OFFSET))
  {
    HEADER_OFFSET <- 0
  }
  
  TempData <- read.csv(file=file_list[1], header=TRUE, sep=",", fileEncoding = "Unicode")
  CompiledData <- TempData[,c(1,2)]
  SDData <- TempData[,c(1,3)]
  
  for (i in 2:length(file_list)){
    TempData <- read.csv(file=file_list[i], header=TRUE, sep=",", fileEncoding = "Unicode")
    CompiledData <- cbind(CompiledData, TempData[,2])
    SDData <- cbind(SDData, TempData[,3])
  }
  
  header <- (c(1:length(file_list))*ALIQUOT_RATIO) + HEADER_OFFSET
  
  names(CompiledData) <- c("Wavelength" ,header)
  names(SDData) <- c("Wavelength", header)
  
  
  #########################################################################################################


  
  if(is.null(X_MAX))
  {
    X_MAX <- 900
  }
  
  if(is.null(X_MIN))
  {
    X_MIN <- 225
  }
  
  
  if(is.null(Y_MAX))
  {
    Y_MAX <- max(CompiledData[36:911,length(file_list)])+.1
  }
  
  
  RealBound<- c(X_MIN - 189, X_MAX - 189)
  RealInterval <- c(X_MIN: X_MAX) - 189
  
  plot(CompiledData[,1], CompiledData[,2],type="l", ylim = c(-.1,Y_MAX), xlim = c(X_MIN, X_MAX))
  for (i in 3:length(file_list)){
    lines(CompiledData[,1], CompiledData[,i],type="l")
  }
  
  
  CompiledData[100,length((file_list))-2]
  
  maxes<- {}
  abs_max <- RealInterval[1]
  j <- 6
  
  
  for(j in RealInterval){
    if(CompiledData[j,length(file_list)-MAX_SEARCH_OFFSET] > CompiledData[j+1,length(file_list)-MAX_SEARCH_OFFSET] & CompiledData[j,length(file_list)-MAX_SEARCH_OFFSET] > CompiledData[j-1,length(file_list)-MAX_SEARCH_OFFSET])
    {
      maxes <- c(maxes,j)
    }
    if(CompiledData[j,length(file_list)-MAX_SEARCH_OFFSET] > CompiledData[abs_max,length(file_list)-MAX_SEARCH_OFFSET])
    {
      abs_max <- j
    }
    
    
  }
  
  points(CompiledData[maxes,1],CompiledData[maxes,length(file_list)-MAX_SEARCH_OFFSET])
##########################################################################################################

  if(!is.null(lambda))
  {
    wl <- (lambda-189)
  }
  if(is.null(lambda))
  {
    wl <- abs_max
  }
  
  
  
wl_data <- as.numeric(CompiledData[wl,2:(length(file_list)+1)])
wl_sd <- as.numeric(SDData[wl,2:(length(file_list)+1)])
wl_data

plot(header, wl_data, xlim = c(0,max(header)), ylim = c(0,max(wl_data)))

arrows(header, wl_data-wl_sd, header, wl_data+wl_sd, length=0.05, angle=90, code=3)

#################################################################################################################
best <- {}


avg_rsq <- {}
rms_rsq <- {}
intersections_x <- {}
intersections_y <- {}

for(i in 1:(length(file_list)))
{
  avg_rsq <- c(avg_rsq, (summary(lm(wl_data[1:i]~header[1:i]))$r.squared + summary(lm(wl_data[(i+1):length(wl_data)]~header[(i+1):length(wl_data)]))$r.squared)/2)
  rms_rsq <- c(rms_rsq, sqrt( (summary(lm(wl_data[1:i]~header[1:i]))$r.squared^2 + summary(lm(wl_data[(i+1):length(wl_data)]~header[(i+1):length(wl_data)]))$r.squared^2)/2 ))

  temp_first_line <- lm(wl_data[1:i]~header[1:i])
  temp_second_line <- lm(wl_data[(i+1):length(wl_data)]~header[(i+1):length(wl_data)])
  
  cm <- rbind(coef(temp_first_line),coef(temp_second_line)) # Coefficient matrix
  intersections_x <- c(intersections_x, c(-solve(cbind(cm[,2],-1)) %*% cm[,1])[1])
  intersections_y <- c(intersections_y, c(-solve(cbind(cm[,2],-1)) %*% cm[,1])[2])
}

plot(header,avg_rsq, ylim = c(.97,1))
#plot(header, -log10(1-avg_rsq))

#plot(header,wl_data)
#plot(diff(wl_data))
#plot(diff(wl_data,differences = 2))

#plot(header,rms_rsq, ylim = c(.97,1))

if(is.null(best))
{
  best <- which(avg_rsq==max(avg_rsq))
}

plot(header, wl_data, xlim = c(0,max(header)), ylim = c(0,max(wl_data)))
arrows(header, wl_data-wl_sd, header, wl_data+wl_sd, length=0.05, angle=90, code=3)

first_line <- lm(wl_data[1:best]~header[1:best])
lines(header[1:best],predict(first_line,data.frame(header[1:best])))

second_line <- lm(wl_data[(best+1):length(wl_data)]~header[(best+1):length(wl_data)])
lines(header[(best+1):length(wl_data)],predict(second_line,data.frame(header[(best+1):length(wl_data)])))

#cm <- rbind(coef(first_line),coef(second_line)) # Coefficient matrix
#intersection <- c(-solve(cbind(cm[,2],-1)) %*% cm[,1])
#points(intersection[1],intersection[2])


plot(header,intersections_x)
points(header[best],intersections_x[best],bg="red",col="red")

#loading capacity vs breakpoint
intersections_x[best]
intersections_y[best]
