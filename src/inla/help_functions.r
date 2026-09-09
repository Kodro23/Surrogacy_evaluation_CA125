# Create categories of trajectories
trend<-function(y,time){
    #' To idendity trend in the time series (decreasing, increasing, convex, concave)
    #'@param y: longitudinal variable
    #'@param time: time variable
    #'@return Category of trajectory
    tryCatch({
        ## Define splines
        #fit spline
        fit1 <- lm(y ~ time+time^2) #square
        fit2 <- lm(y ~ time+time^2+time^3) #cubic

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


#extract coefficients from samples: 
extract_samples <- function(tag) {
  #' Extract coefficients from samples
  #' @param tag: coefficient name
  #' @return Numeric vector of coefficients for the given tag
  idx <- idx_list[[tag]]

  vapply(
    SMP,
    function(s)
      as.numeric(s$latent[idx, 1]),
    numeric(length(idx))
  )
}

#predict trajectory per sample and trial
predict_one_draw <- function(s, j) {
  #' Reconstitute the trajectory for a given sample and trial using estimated fixed and random effects
  #'@param s: sample index
  #'@param j: trial index
  #'@return "Predicited" trajectory for control and treatment arms

  y_ctrl <-theta0[s] +theta3[s] * B1 +theta2[s] * B2 +u0[j, s]
  y_trt <-y_ctrl +(theta4[s] + u1[j, s]) * B1 +(theta5[s] + u2[j, s]) * B2
  list(ctrl = as.numeric(y_ctrl),trt  = as.numeric(y_trt))
}

#compute summaries
compute_trial_summary <- function(s, j) {
  #' Compute trial-specific summaries for a given sample and trial
  #' @param s: sample index
  #' @param j: trial index
  pred <- predict_one_draw(s, j) # compute predicted trajectories
  logHR <- beta_trt_surv[s] + b_trt_surv[j, s]
  y_ctrl <- pred$ctrl
  y_trt  <- pred$trt
  #Retrieve monthly values
  month_days <- c(30, 60, 90, 120, 150, 180)
  nearest_idx <- sapply(month_days,function(t)which.min(abs(time_grid - t))) #closest index to each month

  # -------------------------
  # Nadir
  # -------------------------

  nadir_ctrl <- min(y_ctrl, na.rm = TRUE)
  nadir_trt  <- min(y_trt, na.rm = TRUE)

  delta_nadir <- nadir_trt - nadir_ctrl
  ratio_nadir <- log(nadir_trt / nadir_ctrl)


  # -------------------------
  # Time to nadir
  # -------------------------

  t_nadir_ctrl <- time_grid[which.min(y_ctrl)]
  t_nadir_trt  <- time_grid[which.min(y_trt)]

  delta_time_nadir <- t_nadir_trt - t_nadir_ctrl


  # -------------------------
  # AUC
  # -------------------------
 
  auc_trap <- function(time, y) {
    #' Calculate area under the curve
    #'@param y: longitudinal variable
    #'@param time: time variable
    #'@return Area under the curve

    sum(diff(time) *(head(y, -1) + tail(y, -1)) / 2,na.rm = TRUE)
  }
  auc_ctrl <- auc_trap(time_grid, y_ctrl)
  auc_trt  <- auc_trap(time_grid, y_trt)
  delta_auc <- auc_trt - auc_ctrl


  # -------------------------
  # Percent reductions
  # -------------------------

  pct_ctrl <- (y_ctrl[1] - y_ctrl[nearest_idx]) / y_ctrl[1]
  pct_trt <- (y_trt[1] - y_trt[nearest_idx]) / y_trt[1]

  delta_pct <- pct_trt - pct_ctrl
  ratio_pct <- log(pct_trt / pct_ctrl)


  # -------------------------
  # Monthly slopes
  # -------------------------

  idx0 <- which.min(abs(time_grid - 0)) #Find the grid index closest to day 0
  interval_idx <-c(idx0, nearest_idx) #corresponding number of days

  # Compute slopes for each month
  slope_ctrl <- numeric(6)
  slope_trt  <- numeric(6)

  for (m in 1:6) {

    i0 <- interval_idx[m]
    i1 <- interval_idx[m + 1]

    dt <- time_grid[i1] - time_grid[i0]

    slope_ctrl[m] <-(y_ctrl[i1] - y_ctrl[i0]) / dt
    slope_trt[m] <-(y_trt[i1] - y_trt[i0]) / dt
  }

  delta_slope <- slope_trt - slope_ctrl
  ratio_slope <- log(slope_trt / slope_ctrl)


  # -------------------------
  # Absolute differences
  # -------------------------

  abs_diff <-y_trt[nearest_idx] -y_ctrl[nearest_idx]
  relative_diff <-log(y_trt[nearest_idx]/y_ctrl[nearest_idx])


  #Combine all summaries into a data frame
  data.frame(
    sample = s,
    trialid_num = j,

    logHR=logHR,

    delta_nadir = delta_nadir,
    ratio_nadir = ratio_nadir,
    delta_time_nadir = delta_time_nadir,
    delta_auc = delta_auc,

    delta_pct_m1 = delta_pct[1],
    delta_pct_m2 = delta_pct[2],
    delta_pct_m3 = delta_pct[3],
    delta_pct_m4 = delta_pct[4],
    delta_pct_m5 = delta_pct[5],
    delta_pct_m6 = delta_pct[6],

    ratio_pct_m1 = ratio_pct[1],
    ratio_pct_m2 = ratio_pct[2], 
    ratio_pct_m3 = ratio_pct[3],
    ratio_pct_m4 = ratio_pct[4],
    ratio_pct_m5 = ratio_pct[5],
    ratio_pct_m6 = ratio_pct[6],

    delta_slope_m1 = delta_slope[1],
    delta_slope_m2 = delta_slope[2],
    delta_slope_m3 = delta_slope[3],
    delta_slope_m4 = delta_slope[4],
    delta_slope_m5 = delta_slope[5],
    delta_slope_m6 = delta_slope[6],

    ratio_slope_m1 = ratio_slope[1],
    ratio_slope_m2 = ratio_slope[2],
    ratio_slope_m3 = ratio_slope[3],
    ratio_slope_m4 = ratio_slope[4],
    ratio_slope_m5 = ratio_slope[5],
    ratio_slope_m6 = ratio_slope[6],

    absolutediff_m1 = abs_diff[1],
    absolutediff_m2 = abs_diff[2],
    absolutediff_m3 = abs_diff[3],
    absolutediff_m4 = abs_diff[4],
    absolutediff_m5 = abs_diff[5],
    absolutediff_m6 = abs_diff[6],

    relativediff_m1 = relative_diff[1],
    relativediff_m2 = relative_diff[2],
    relativediff_m3 = relative_diff[3],
    relativediff_m4 = relative_diff[4],
    relativediff_m5 = relative_diff[5],
    relativediff_m6 = relative_diff[6]
  )
}