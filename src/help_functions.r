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




#Calculate AUC
auc_trap <- function(time, y) {

  sum(
    diff(time) *
      (head(y, -1) + tail(y, -1)) / 2,
    na.rm = TRUE
  )
}

#extract coefficients from samples: 
extract_samples <- function(tag) {

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

  y_ctrl <-
    beta0[s] +
    beta1[s] * B1 +
    beta2[s] * B2 +
    beta3[s] * B3 +
    u0[j, s]

  y_trt <-
    y_ctrl +
    (gamma1[s] + u1[j, s]) * B1 +
    (gamma2[s] + u2[j, s]) * B2 +
    (gamma3[s] + u3[j, s]) * B3

  list(
    ctrl = as.numeric(y_ctrl),
    trt  = as.numeric(y_trt)
  )
}

#compute summaries
compute_trial_summary <- function(s, j) {

  pred <- predict_one_draw(s, j)

  y_ctrl <- pred$ctrl
  y_trt  <- pred$trt

  # -------------------------
  # Nadir
  # -------------------------

  nadir_ctrl <- min(y_ctrl, na.rm = TRUE)
  nadir_trt  <- min(y_trt, na.rm = TRUE)

  delta_nadir <- nadir_trt - nadir_ctrl


  # -------------------------
  # Time to nadir
  # -------------------------

  t_nadir_ctrl <- time_grid[which.min(y_ctrl)]
  t_nadir_trt  <- time_grid[which.min(y_trt)]

  delta_time_nadir <-
    t_nadir_trt - t_nadir_ctrl


  # -------------------------
  # AUC
  # -------------------------

  auc_ctrl <- auc_trap(time_grid, y_ctrl)
  auc_trt  <- auc_trap(time_grid, y_trt)

  delta_auc <- auc_trt - auc_ctrl


  # -------------------------
  # Percent reductions
  # -------------------------

  # CA125 is log-transformed in the model.
  # For a genuine percentage reduction,
  # go back to the original CA125 scale.

  ca_ctrl <- exp(y_ctrl)
  ca_trt  <- exp(y_trt)

  pct_ctrl <- (
    ca_ctrl[1] - ca_ctrl[nearest_idx]
  ) / ca_ctrl[1]

  pct_trt <- (
    ca_trt[1] - ca_trt[nearest_idx]
  ) / ca_trt[1]

  delta_pct <- pct_trt - pct_ctrl


  # -------------------------
  # Monthly slopes
  # -------------------------

  idx0 <- which.min(abs(time_grid - 0))

  interval_idx <-
    c(idx0, nearest_idx)

  slope_ctrl <- numeric(6)
  slope_trt  <- numeric(6)

  for (m in 1:6) {

    i0 <- interval_idx[m]
    i1 <- interval_idx[m + 1]

    dt <- time_grid[i1] - time_grid[i0]

    slope_ctrl[m] <-
      (y_ctrl[i1] - y_ctrl[i0]) / dt

    slope_trt[m] <-
      (y_trt[i1] - y_trt[i0]) / dt
  }

  delta_slope <- slope_trt - slope_ctrl


  # -------------------------
  # Absolute differences
  # -------------------------

  abs_diff <-
    y_trt[nearest_idx] -
    y_ctrl[nearest_idx]


  # -------------------------
  # Trial-specific treatment
  # spline coefficients
  # Henderson-type variables
  # -------------------------

  effect1 <- gamma1[s] + u1[j, s]
  effect2 <- gamma2[s] + u2[j, s]
  effect3 <- gamma3[s] + u3[j, s]


  data.frame(
    sample = s,
    trialid_num = j,

    delta_nadir = delta_nadir,
    delta_time_nadir = delta_time_nadir,
    delta_auc = delta_auc,

    pct_m1 = delta_pct[1],
    pct_m2 = delta_pct[2],
    pct_m3 = delta_pct[3],
    pct_m4 = delta_pct[4],
    pct_m5 = delta_pct[5],
    pct_m6 = delta_pct[6],

    slope_m1 = delta_slope[1],
    slope_m2 = delta_slope[2],
    slope_m3 = delta_slope[3],
    slope_m4 = delta_slope[4],
    slope_m5 = delta_slope[5],
    slope_m6 = delta_slope[6],

    absolutediff_m1 = abs_diff[1],
    absolutediff_m2 = abs_diff[2],
    absolutediff_m3 = abs_diff[3],
    absolutediff_m4 = abs_diff[4],
    absolutediff_m5 = abs_diff[5],
    absolutediff_m6 = abs_diff[6],

    effect_trt_ca1251 = effect1,
    effect_trt_ca1252 = effect2,
    effect_trt_ca1253 = effect3
  )
}