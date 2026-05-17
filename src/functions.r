# Create categories of trajectories
# decreasing, increasing, convex, concave
##Splines for an observation over 3 months: 3 internal knots and 2 external at the 2% and 98% percentiles of measurment times
trend<-function(y,time){
    #' To idendity trend in the time series
    #'@param y: measurments variable
    #'@param time: time variable
    #'@return Caterory of trajectory
    tryCatch({
        ## Define splines
        #fit spline
        fit1 <- lm(y ~ time+time^2)
        fit2 <- lm(y ~ time+time^2+time^3)
        # ##fit sinusoidal function
        # period=365
        # fit2 <- lm(y ~ sin(2*pi*time/period) + cos(2*pi*time/period))

        ##Predictions
        grid_time <- seq(min(time), max(time), length.out = 500)
        y_pred1 <- predict(fit1, newdata = data.frame(time = grid_time))
        y_pred2 <- predict(fit2, newdata = data.frame(time = grid_time))

        ##Model performances
        if(AIC(fit2)<= AIC(fit1)){
                dy <- diff(y_pred2) / diff(grid_time) #first derivative
                d2y <- diff(dy) / diff(grid_time[-1]) #second derivative
                if(all(dy <= 0) & all(abs(d2y) < 1e-2)) {
                        return ("decreasing")
                }else if(all(dy >= 0) & all(abs(d2y) < 1e-2)) {
                        return ("increasing")
                } else if (all(d2y<0)){
                        return ("concave")
                }else if (all(d2y>0)) {
                        return ("convex")
                }
        }else{
                dy <- diff(y_pred1) / diff(grid_time) #first derivative
                d2y <- diff(dy) / diff(grid_time[-1]) #second derivative
                if(all(dy <= 0) & all(abs(d2y) < 1e-3)) {
                        return ("decreasing")
                }else if(all(dy >= 0) & all(abs(d2y) < 1e-3)) {
                        return ("increasing")
                } else if (all(d2y<0) & all(d2y>= -1e-3) ){
                        return ("concave")
                }else if (all(d2y>=1e-3)) {
                        return ("convex")
                }
        }
        return("unidentified")
}, error = function(e) {
    # If anything fails, return unidentified
    return("unidentified")
  })
}


auc_manual <- function(response, predictor) {
  # remove NA
  df <- data.frame(response, predictor)
  df <- df[complete.cases(df), ]
  
  # rank predictions
  r <- rank(df$predictor)
  
  n1 <- sum(df$response == 1)
  n0 <- sum(df$response == 0)
  
  auc <- (sum(r[df$response == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  return(auc)
}
