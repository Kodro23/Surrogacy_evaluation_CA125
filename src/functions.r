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

# -----------------------------
# Helper functions to compute summaries
# -----------------------------

get_monthly_slopes <- function(time_grid, y_ctrl, y_trt, months = 1:6) {

  slopes <- lapply(months, function(m) {

    target_time <- m * 30
    idx1 <- which.min(abs(time_grid - target_time)) # current month
    idx0 <- which.min(abs(time_grid - (target_time - 30))) # previous month

    if (idx1 == idx0) {return(c(trt = NA, ctrl = NA))}
    dy_c <- (y_ctrl$pred[idx1] - y_ctrl$pred[idx0]) /(time_grid[idx1] - time_grid[idx0])
    dy_t <- (y_trt$pred[idx1] - y_trt$pred[idx0]) /(time_grid[idx1] - time_grid[idx0])
    list(trt = dy_t, ctrl = dy_c)
  })
  slopes
}


compute_auc <- function(y, time_grid) {
  sum(diff(time_grid) * (head(y$pred[,1], -1) + tail(y$pred[,1], -1)) / 2)
}
get_weekly_nadir <- function(time_grid, pred, week_size = 7) {

  week_id <- floor(time_grid / week_size)
  aggregate(pred[,1],by = list(week = week_id),FUN = min,na.rm = TRUE)
  }

get_monthly_pct_reduction <- function(time_grid, y_ctrl, y_trt, months = 1:6) {

  baseline_ctrl <- y_ctrl$pred[1, 1]
  baseline_trt  <- y_trt$pred[1, 1]
  pct_reduction <- lapply(months, function(m) {
    idx <- which.min(abs(time_grid - m * 30))
    pct_ctrl <- (baseline_ctrl - y_ctrl$pred[idx, 1]) / baseline_ctrl
    pct_trt  <- (baseline_trt - y_trt$pred[idx, 1]) / baseline_trt
    list(trt  = pct_trt,ctrl = pct_ctrl)
  })

  names(pct_reduction) <- paste0("month", months)

  pct_reduction
}
# -----------------------------
# Helper functions for bootstrap
# -----------------------------
get_r2_models <- function(data, formulas = model_formulas) {

  sapply(formulas, function(f) {

    tryCatch({

      model <- suppressWarnings(lm(f,data = data,weights = n_events))
      summary(model)$r.squared
    }, error = function(e) {
      NA_real_
    })
  })
}
# -----------------------------
# Helper functions to plot LOSO
# -----------------------------

plot_loso <- function(summary_name, loso_df) {

  df <- loso_df %>%
    dplyr::filter(summary == summary_name) %>%
    dplyr::arrange(R2)

  ggplot(
    df,
    aes(
      x = reorder(left_out_trial, R2),
      y = R2
    )
  ) +
    geom_point(size = 3) +
    geom_segment(
      aes(
        x = left_out_trial,
        xend = left_out_trial,
        y = 0,
        yend = R2
      )
    ) +
    coord_flip() +
    labs(
      title = paste0(
        "Leave-one-study-out analysis: ",
        summary_name
      ),
      x = "Trial omitted",
      y = expression(R^2)
    ) +
    theme_bw()
}