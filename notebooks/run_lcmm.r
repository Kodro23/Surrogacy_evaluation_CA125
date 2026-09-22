#Set directory
install.packages("here")
library(here)
setwd(here())
#Install and load libraires
invisible(source("./src/install_librairies.r"))
install_packages(directory)
invisible(source("./src/load_librairies.r"))
# Load functions
invisible(source("./src/functions.r"))

#Load dataset
ca125 <- read.csv("./data/MAOV_ca125_all.csv", sep=";")
ca125 <- ca125%>%mutate(id=paste0(trialid,"_",patid)) 
ca125 %>% summarize(
    n_patients=n_distinct(id),
    n_trials=n_distinct(trialid),
    n_measures=n()
)
#Clean dataset
#check time between first surgery and first measure
Sys.setlocale("LC_TIME", "C")
ca125 <- ca125 %>%
  mutate(
    dosurg1_date = as.Date(toupper(dosurg1), format = "%d%b%Y"),
    do_CA125_date = as.Date(toupper(do_CA125), format = "%d%b%Y")
  ) %>%
  group_by(id) %>%
  mutate(
    delay_surg_ca125 = if (all(is.na(do_CA125_date)) ||
                           all(is.na(dosurg1_date))) {
      NA_real_
    } else {
      as.numeric(
        min(do_CA125_date, na.rm = TRUE) -
          min(dosurg1_date, na.rm = TRUE)
      )
    }
  ) %>%
  ungroup() %>%
  filter(is.na(delay_surg_ca125) | delay_surg_ca125 >= 7)
# Clean first progression column
ca125 <- ca125 %>%
  mutate(
    dor = as.Date(toupper(dor), format = "%d%b%Y"),
    dofirst_prog = as.Date(toupper(dofirst_prog), format = "%d%b%Y"),
    dofirst_prog = if_else(
      is.na(dofirst_prog) & !is.na(dor),
      dor + pfs_time_m * 30.4375,
      dofirst_prog
    )
  )
ca125 <- ca125 %>%
  mutate(
    dotrt_begin = as.Date(toupper(dotrt_begin), format = "%d%b%Y"),
    donext_trt = as.Date(toupper(donext_trt), format = "%d%b%Y"),
    do_CA125 = as.Date(toupper(do_CA125), format = "%d%b%Y")
  )
#delete measures after progression 
ca125<-ca125 %>%
  filter( (trialid != "MRC-ICON7" & do_CA125 <= dofirst_prog) |
        (trialid == "MRC-ICON7" & time_ca125<=pfs_time_m * 30.4375))
#Create variables
ca125$trt_bin<- if_else(ca125$trt=="Standard regimen", 0,1)
ca125["os_time_d"]<-ca125$os_time_m*30.4375
ca125["pfs_time_d"]<-ca125$pfs_time_m*30.4375
ca125["os_statut"]<-factor(ifelse(ca125$os_status=='Yes', 1, 0), levels=c(0,1))
ca125["pfs_statut"]<-factor(ifelse(ca125$PFS_status=='Yes', 1, 0), levels=c(0,1))
#restrict data                             
colonnes<-c("trialid","patid","id",'trt','dor','CA125_mes','do_CA125','id_temp','dotrt_begin', "donext_trt",'dotrt_end',"time_ca125",'os_status','OS_time','PFS_status','PFS_time',"nb_mesures",
"ca125_base",'Visit','pfs_time_m',"pfs_time_d",'os_time_m',"os_time_d",'group','country','maov_group','dob','age','figo','figo4c',"trt_bin","os_statut","pfs_statut")
ca125_all_6 <- ca125[colonnes]%>%filter(time_ca125 <= 188)
ca125_all_6["id_num"]<- as.numeric(as.factor(ca125_all_6$id))
ca125_all_6["trialid_num"]<- as.numeric(as.factor(ca125_all_6$trialid))
ca125_all_6["log_ca125"]<-log(ca125_all_6$CA125_mes)


#Set up model
internal_knots <- quantile(ca125_all_6$time_ca125, probs = c(0.22,0.66))
boundary_knots <- quantile(ca125_all_6$time_ca125, probs = c(0.18,0.98))
ns_basis <- ns(ca125_all_6$time_ca125,knots = internal_knots,Boundary.knots = boundary_knots)
ca125_all_6["id_num"]<- as.integer(as.numeric(as.factor(ca125_all_6$id)))
ca125_all_6["trialid_num"]<- as.integer(as.numeric(as.factor(ca125_all_6$trialid)))
ca125_all_6["log_ca125"]<-log(ca125_all_6$CA125_mes)
#Load kelim data
kelim <- read.csv("./data/data_GCIG_KELIM_survie.csv", sep=";")
kelim <- kelim%>%rename(id=ID_unique)%>%
                 mutate(trt_bin=as.factor(if_else(ttt_arm=="ctrl", 0,1))) %>%
                 select(id,trial_name,trt_bin, KELIM_continuous,KELIM_score,OS_status,OS_time,PFS_status,PFS_time)
print(paste0(c("Number of patients:", nrow(kelim))))
print(paste0(c("Number of trials:", length(unique(kelim$trial_name)))))

ca125_all_6<-ca125_all_6 %>%arrange(id_num)
trial_list <- split(ca125_all_6, ca125_all_6$trialid)
time_grid <- seq(min(ca125_all_6$time_ca125),max(ca125_all_6$time_ca125),length.out = 200)
# -----------------------------
# storage
# -----------------------------

results_list <- list()
weekly_predictions<-list()
average_values<-list()
model_store <- list()
nadirs<-list()
nadirs_week<-list()
wald_list <-list()

# -----------------------------
# loop
# -----------------------------

for (trial_name in names(trial_list)) {

  cat("Fitting:", trial_name, "\n")
  
  # Filter dataset for the current trial
  df <- trial_list[[trial_name]]

  # -------------------------
  # spline basis
  # -------------------------
  internal_knots <- quantile(df$time_ca125, c(0.22, 0.66))
  boundary_knots <- quantile(df$time_ca125, c(0.18, 0.98))

  ns_basis <- ns(df$time_ca125,
                 knots = internal_knots,
                 Boundary.knots = boundary_knots)

  colnames(ns_basis) <- paste0("ns", 1:ncol(ns_basis))
  df <- cbind(df, ns_basis)
  df$week <- round(df$time_ca125/7)

  # -------------------------
  # model
  # -------------------------
  model <- tryCatch({

    Jointlcmm(
      fixed = log_ca125 ~ ns1 + ns2 + ns3 +
        trt_bin:(ns1 + ns2 + ns3),
      random = ~ ns1 + ns2 + ns3,
      subject = "id_num",
      ng = 1,
      survival = Surv(pfs_time_d, pfs_statut) ~ trt_bin,
      hazard = "splines",
      data = df,
      maxiter = 1500
    )
  
  }, error = function(e) NULL)

  if (is.null(model)) next
  model_coefs<-model$best
  log_ca125_pred <- predictY(model,newdata = df,var.time = "time_ca125")$pred
  weekly_prediction <- tryCatch({
    ca125_patient <- data.frame(trialid = trial_name,
                                trt=df$trt_bin,
                                log_ca125=as.vector(df$log_ca125),
                                week=df$week, 
                                log_ca125_pred =as.vector(log_ca125_pred)) %>%
      group_by(trialid, trt, week) %>%
      summarise(
        log_ca125_pred_week =
          mean(
            log_ca125_pred,
            na.rm = TRUE
          ),
        log_ca125_week =
          mean(
            log_ca125,
            na.rm = TRUE
          ),
        .groups = "drop"
      )
    
  }, error = function(e) {
    
    cat(
      "❌ Failed:",
      trial_name,
      "\n"
    )
    
    print(e$message)
    
    NULL
  })

  average_value <- tryCatch({
    ca125_patient <- data.frame(id_num = df$id_num,log_ca125_pred =as.vector(log_ca125_pred)) %>%
      group_by(id_num) %>%
      summarise(
        log_ca125_pred =
          mean(
            log_ca125_pred,
            na.rm = TRUE
          ),
        .groups = "drop"
      )
    
    # -------------------------
    # Number of events
    # -------------------------
    n_events_trial <- df %>%
      distinct(
        id_num,
        .keep_all = TRUE
      ) %>%
      summarise(
        n_events =
          sum(
            pfs_statut == 1,
            na.rm = TRUE
          )
      ) %>%
      pull(n_events)
    
    # -------------------------
    # Trial summary
    # -------------------------
    trial_summary <- ca125_patient %>%
      summarise(
        mean_log_ca125 =
          mean(
            log_ca125_pred,
            na.rm = TRUE
          ),
        .groups = "drop"
      ) %>%
      mutate(
        trialid = trial_name,
        log_hr = model_coefs["trt_bin"][1],
        hr = exp(model_coefs["trt_bin"][1]),
        n_events = n_events_trial
      )
    
    trial_summary
    
  }, error = function(e) {
    
    cat(
      "❌ Failed:",
      trial_name,
      "\n"
    )
    
    print(e$message)
    
    NULL
  })
  
  
  weekly_predictions[[trial_name]] <- weekly_prediction
  model_store[[trial_name]] <- model
  average_values[[trial_name]] <- average_value
  wald_list[[trial_name]] <-WaldMult(model, pos=c(13,14,15))
  beta_longitudinal_trt <- model_coefs[grep("ns[1-3]:trt_bin", names(model_coefs))]

 

  # -------------------------
  # predictions
  # -------------------------
  y_ctrl <- predictY(model, newdata = pred_ctrl, var.time = "time_ca125")
  y_trt  <- predictY(model, newdata = pred_trt,  var.time = "time_ca125")

  # -------------------------
  # NADIR
  # -------------------------
  nadir_ctrl <- min(y_ctrl$pred[,1], na.rm=TRUE)
  nadir_trt  <- min(y_trt$pred[,1], na.rm=TRUE)
  delta_nadir <- nadir_trt - nadir_ctrl
  ratio_delta_nadir <- nadir_trt/nadir_ctrl
  nadirs[[trial_name]]<-data.frame(nadir_trt=nadir_trt,nadir_ctrl=nadir_ctrl)

  # -------------------------
  # weekly NADIR
  # -------------------------
  weekly_ctrl <- get_weekly_nadir(time_grid, y_ctrl$pred, week_size = 7)
  weekly_trt  <- get_weekly_nadir(time_grid, y_trt$pred,  week_size = 7)

  # global nadir (from weekly values)

  nadirs_week[[trial_name]] <- data.frame(
    nadir_trt = min(weekly_trt$x, na.rm = TRUE),
    nadir_ctrl = min(weekly_ctrl$x, na.rm = TRUE)
  )
  # -------------------------
  # time-to-nadir
  # -------------------------
  t_nadir_ctrl <- time_grid[which.min(y_ctrl$pred[,1])]
  t_nadir_trt  <- time_grid[which.min(y_trt$pred[,1])]
  delta_t_nadir <- t_nadir_trt - t_nadir_ctrl

  # -------------------------
  # AUC
  # -------------------------
  auc_ctrl <- compute_auc(y_ctrl, time_grid)
  auc_trt  <- compute_auc(y_trt, time_grid)
  delta_auc <- auc_trt - auc_ctrl

  # -------------------------
  # percent reduction
  # -------------------------
  pcts<- get_monthly_pct_reduction(
    time_grid,
    y_ctrl,
    y_trt,
    months = 1:6)

  # -------------------------
  # 6-month slopes
  # -------------------------
  slopes <- get_monthly_slopes(time_grid, y_ctrl, y_trt)

  # -------------------------
  # Differences
  # -------------------------

  #1 month
  idx_30  <- which.min(abs(time_grid - 30))
  D_30 <- list(trt=y_trt$pred[idx_30, 1], ctrl=y_ctrl$pred[idx_30, 1])
  #2 months
  idx_60  <- which.min(abs(time_grid - 60))
  D_60 <- list(trt=y_trt$pred[idx_60, 1], ctrl=y_ctrl$pred[idx_60, 1])
  #3 months
  idx_90  <- which.min(abs(time_grid - 90))
  D_90 <- list(trt=y_trt$pred[idx_90, 1], ctrl=y_ctrl$pred[idx_90, 1])
  #4 months
  idx_120  <- which.min(abs(time_grid - 120))
  D_120 <- list(trt=y_trt$pred[idx_120, 1], ctrl=y_ctrl$pred[idx_120, 1])
  #5 months
  idx_150  <- which.min(abs(time_grid - 150))
  D_150 <- list(trt=y_trt$pred[idx_150, 1], ctrl=y_ctrl$pred[idx_150, 1])
  #6 months
  idx_180 <- which.min(abs(time_grid - 180))
  D_180 <- list(trt=y_trt$pred[idx_180, 1], ctrl=y_ctrl$pred[idx_180, 1])

  #Effect of treatment on Kelim:
  kelim_data<- kelim[kelim$trial_name==trial_name,]
  if (nrow(kelim_data)>0){
    linear_reg<-lm(KELIM_continuous ~ trt_bin, data = kelim_data)
    beta_kelim<-linear_reg$coefficients["trt_bin"]
    delta_kelim<-mean(kelim_data$KELIM_continuous[kelim_data$trt_bin==1], na.rm=TRUE)-mean(kelim_data$KELIM_continuous[kelim_data$trt_bin==0], na.rm=TRUE)
    ratio_kelim<-mean(kelim_data$KELIM_continuous[kelim_data$trt_bin==1], na.rm=TRUE)/mean(kelim_data$KELIM_continuous[kelim_data$trt_bin==0], na.rm=TRUE)
  }else{
    beta_kelim<-NA
    delta_kelim<-NA
    ratio_kelim<-NA
  }
  p_logHR <- 2 * pnorm(-abs(model_coefs["trt_bin"] /sqrt(diag(vcov(model)))["trt_bin"])) #pvalue for HR
  se_trt_ca125 <- sqrt(diag(vcov(model)))[names(beta_longitudinal_trt)] #standard error
  p_trt_ca125 <- 2 * pnorm(-abs(beta_longitudinal_trt / se_trt_ca125)) #pvalue for effect on CA125

  # -------------------------
  # store
  # -------------------------
  results_list[[trial_name]] <- data.frame(
    trialid = trial_name,
    logHR = model_coefs["trt_bin"] ,
    SE_logHR = sqrt(diag(vcov(model)))["trt_bin"],
    p_logHR = p_logHR,
    effect_trt_ca1251 = beta_longitudinal_trt[1],
    effect_trt_ca1252 = beta_longitudinal_trt[2],
    effect_trt_ca1253 = beta_longitudinal_trt[3],
    p_effect_trt_ca1251 = p_trt_ca125[1],
    p_effect_trt_ca1252 = p_trt_ca125[2],
    p_effect_trt_ca1253 = p_trt_ca125[3],
    beta_kelim = beta_kelim,
    delta_kelim = delta_kelim,
    ratio_kelim = ratio_kelim,

    delta_nadir = delta_nadir,
    delta_time_nadir = delta_t_nadir,
    delta_auc = delta_auc,
    
    pct_m1 = I(list(pcts[1])),
    pct_m2 = I(list(pcts[2])),
    pct_m3 = I(list(pcts[3])),
    pct_m4 = I(list(pcts[4])),
    pct_m5 = I(list(pcts[5])),
    pct_m6 = I(list(pcts[6])),

    diff_m1=I(list(D_30)),
    diff_m2=I(list(D_60)),
    diff_m3=I(list(D_90)),
    diff_m4=I(list(D_120)),
    diff_m6=I(list(D_150)),
    diff_m5=I(list(D_180)),

    slope_m1 = I(list(slopes[1])),
    slope_m2 = I(list(slopes[2])),
    slope_m3 = I(list(slopes[3])),
    slope_m4 = I(list(slopes[4])),
    slope_m5 = I(list(slopes[5])),
    slope_m6 = I(list(slopes[6])),

    n=nrow(df),
    n_events = nrow(df%>%filter(pfs_statut==1)%>% distinct(id))
  )
}

#Save results
save(model_store, file = "results/summaries/model_store.RData")
save(wald_list, file = "results/summaries/wald_list.RData")
save(results_list, file = "results/summaries/results_list.RData")
save(weekly_predictions, file = "results/summaries/weekly_predictions.RData")
save(average_values, file = "results/summaries/average_values.RData")
trial_results <- bind_rows(results_list)
trial_weekly_predictions <- bind_rows(weekly_predictions)
trial_summary <- bind_rows(average_values)


#Surrogacy
# Format
trial_results_ <- trial_results %>%
  mutate(

    # % reduction
    delta_pct1 = map_dbl(pct_m1, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_pct2 = map_dbl(pct_m2, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_pct3 = map_dbl(pct_m3, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_pct4 = map_dbl(pct_m4, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_pct5 = map_dbl(pct_m5, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_pct6 = map_dbl(pct_m6, ~ .x[[1]]$trt - .x[[1]]$ctrl),

    ratio_pct1 = map_dbl(pct_m1, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_pct2 = map_dbl(pct_m2, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_pct3 = map_dbl(pct_m3, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_pct4 = map_dbl(pct_m4, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_pct5 = map_dbl(pct_m5, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_pct6 = map_dbl(pct_m6, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),

    # Slopes
    delta_slope1 = map_dbl(slope_m1, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_slope2 = map_dbl(slope_m2, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_slope3 = map_dbl(slope_m3, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_slope4 = map_dbl(slope_m4, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_slope5 = map_dbl(slope_m5, ~ .x[[1]]$trt - .x[[1]]$ctrl),
    delta_slope6 = map_dbl(slope_m6, ~ .x[[1]]$trt - .x[[1]]$ctrl),

    ratio_slope1 = map_dbl(slope_m1, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_slope2 = map_dbl(slope_m2, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_slope3 = map_dbl(slope_m3, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_slope4 = map_dbl(slope_m4, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_slope5 = map_dbl(slope_m5, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),
    ratio_slope6 = map_dbl(slope_m6, ~ log(.x[[1]]$trt / .x[[1]]$ctrl)),

    # Absolute CA125 values
    abs1 = map_dbl(diff_m1, ~ .x$trt - .x$ctrl),
    abs2 = map_dbl(diff_m2, ~ .x$trt - .x$ctrl),
    abs3 = map_dbl(diff_m3, ~ .x$trt - .x$ctrl),
    abs4 = map_dbl(diff_m4, ~ .x$trt - .x$ctrl),
    abs5 = map_dbl(diff_m5, ~ .x$trt - .x$ctrl),
    abs6 = map_dbl(diff_m6, ~ .x$trt - .x$ctrl),

    relativ1 = map_dbl(diff_m1, ~ log(.x$trt / .x$ctrl)),
    relativ2 = map_dbl(diff_m2, ~ log(.x$trt / .x$ctrl)),
    relativ3 = map_dbl(diff_m3, ~ log(.x$trt / .x$ctrl)),
    relativ4 = map_dbl(diff_m4, ~ log(.x$trt / .x$ctrl)),
    relativ5 = map_dbl(diff_m5, ~ log(.x$trt / .x$ctrl)),
    relativ6 = map_dbl(diff_m6, ~ log(.x$trt / .x$ctrl))
    )
# Compute r2

model_formulas <- list(
  nadir = logHR ~ delta_nadir,
  time_nadir = logHR ~ delta_time_nadir,
  auc = logHR ~ delta_auc,
  delta_kelim = logHR ~ delta_kelim,
  ratio_kelim = logHR ~ ratio_kelim,

  delta_pct1 = logHR ~ delta_pct1,
  ratio_pct1 = logHR ~ ratio_pct1,
  delta_pct2 = logHR ~ delta_pct2,
  ratio_pct2 = logHR ~ ratio_pct2,
  delta_pct3 = logHR ~ delta_pct3,
  ratio_pct3 = logHR ~ ratio_pct3,
  delta_pct4 = logHR ~ delta_pct4,
  ratio_pct4 = logHR ~ ratio_pct4,
  delta_pct5 = logHR ~ delta_pct5,
  ratio_pct5 = logHR ~ ratio_pct5,
  delta_pct6 = logHR ~ delta_pct6,
  ratio_pct6 = logHR ~ ratio_pct6,

  delta_slope1 = logHR ~ delta_slope1,
  ratio_slope1 = logHR ~ ratio_slope1,
  delta_slope2 = logHR ~ delta_slope2,
  ratio_slope2 = logHR ~ ratio_slope2,
  delta_slope3 = logHR ~ delta_slope3,
  ratio_slope3 = logHR ~ ratio_slope3,
  delta_slope4 = logHR ~ delta_slope4,
  ratio_slope4 = logHR ~ ratio_slope4,
  ratio_slope6 = logHR ~ ratio_slope6,
  delta_slope5 = logHR ~ delta_slope5,
  delta_slope6 = logHR ~ delta_slope6,

  abs1 = logHR ~ abs1,
  relativ1 = logHR ~ relativ1,
  abs2 = logHR ~ abs2,
  relativ2 = logHR ~ relativ2,
  abs3 = logHR ~ abs3,
  relativ3 = logHR ~ relativ3,
  abs4 = logHR ~ abs4,
  relativ4 = logHR ~ relativ4,
  abs5 = logHR ~ abs5,
  relativ5 = logHR ~ relativ5,
  abs6 = logHR ~ abs6,
  relativ6 = logHR ~ relativ6
)

#Boostrap
set.seed(123)
B <- 2000
n <- nrow(trial_results_)
boot_r2 <- matrix(NA_real_,nrow = B,ncol = length(model_formulas))
colnames(boot_r2) <- names(model_formulas)
for (b in seq_len(B)) {
  idx <- sample.int(n,size = n,replace = TRUE)
  data_b <- trial_results_[idx, ]
  boot_r2[b, ] <- get_r2_models(data_b)
}
boot_summary <- data.frame(
  model = colnames(boot_r2),
  R2 = colMeans(boot_r2,na.rm = TRUE),
  R2_low = apply(boot_r2,2,quantile,probs = 0.025,na.rm = TRUE),
  R2_high = apply(boot_r2,2,quantile,probs = 0.975,na.rm = TRUE)
)

#KELIM
id_map <- ca125_all_6 %>%distinct(id_num, id)
trial_map <- ca125_all_6 %>%distinct(trialid_num, trialid)
id_map <- ca125_all_6 %>%distinct(id_num, id)
kelim_data <- kelim %>%
  filter(!is.na(KELIM_continuous),!is.na(PFS_time))%>%
  mutate( 
         pfs_time_d=as.numeric(PFS_time*365.2425 ),
         pfs_statut=factor(ifelse(PFS_status=='Yes', 1, 0), levels=c(0,1))
         ) %>%
  mutate(
    trt_bin = unname(trt_bin),
    pfs_time_d = unname(pfs_time_d),
    pfs_statut = unname(pfs_statut),
    KELIM_continuous = unname(KELIM_continuous)
  ) %>%left_join(trial_map,by = c("trial_name"= "trialid")) %>%left_join(id_map,by = "id")  %>%
  mutate(
    trialid_num = unname(trialid_num),
    id_num = unname(id_num))
kelim_data$trt_bin <- as.numeric(as.character(kelim_data$trt_bin))
kelim_data$pfs_statut <- as.numeric(as.character(kelim_data$pfs_statut))
print(paste0(c("Number of measurements with KELIM data:", nrow(kelim_data), "\nNumber of patients with KELIM data:", length(unique(kelim_data$id)))))
set.seed(123)

B <- 200

trials <- unique(kelim_data$trialid_num)
n_trials <- length(trials)
boot_r2 <- matrix(NA_real_,nrow = B,ncol = 1)
colnames(boot_r2) <- "plackett_weighted"

for (b in seq_len(B)) {

  sampled_trials <- sample(trials,size = n_trials,replace = TRUE)
  data_b <- do.call(
    rbind,
    lapply(seq_along(sampled_trials), function(j) {
      trial_j <- kelim_data[kelim_data$trialid_num == sampled_trials[j],]
      # Give each sampled copy a new trial ID
      trial_j$trialid_num <- j
      # Give patients unique IDs across bootstrap trial copies
      trial_j$id_num <- paste0(j, "_", trial_j$id_num)
      trial_j
    })
  )

  # --------------------------------
  # Fit model
  # --------------------------------
  model <- tryCatch(
    MetaAnalyticSurvCont(
      data = data_b,
      true = pfs_time_d,
      trueind = pfs_statut,
      surrog = KELIM_continuous,
      trt = trt_bin,
      center = trialid_num,
      trial = trialid_num,
      patientid = id_num,
      copula = "Plackett",
      adjustment = "weighted"
    ),
    error = function(e) NULL
  )

  # --------------------------------
  # Extract R2
  # --------------------------------
  if (!is.null(model)) {
    boot_r2[b, 1] <- model$Trial.R2["R2 Trial (weighted)"][[1]]
}
}
boot_summary_kelim <- data.frame(
  model = colnames(boot_r2),
  R2 = colMeans(boot_r2,na.rm = TRUE),
  R2_low = apply(boot_r2,2,quantile,probs = 0.025,na.rm = TRUE),
  R2_high = apply(boot_r2,2,quantile,probs = 0.975,na.rm = TRUE)
)
save(boot_summary_kelim, file = "results/summaries/boot_summary_kelim.RData")

