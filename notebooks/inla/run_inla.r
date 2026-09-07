# Set directory
install.packages("here")
library(here)
directory=here()
setwd(directory)

#Install packages
invisible(source("./src/inla/nstall_librairies.r"))
install_packages(directory)
# Load libraries
invisible(source("./src/load_librairies.r"))

# Load functions
invisible(source("./src/help_functions.r"))

#Load dataset
ca125 <- read.csv("./data/MAOV_ca125_all.csv", sep=";")
ca125 <- ca125%>%mutate(id=paste(trialid,"_",patid)) 

#Clean dataset
#check time between first surgery and first measure
Sys.setlocale("LC_TIME", "C")
ca125$delay_surg_trt<-(as.Date(toupper(ca125$dotrt_begin), format = "%d%b%Y")-as.Date(toupper(ca125$dosurg1), format = "%d%b%Y"))
#delete patients with delay <7
ca125<- ca125[!ca125$id %in% unique(ca125$id[ca125$delay_surg_trt < 7 & !is.na(ca125$delay_surg_trt)]),]
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
ca125_all_6 <- ca125[colonnes]%>%
  filter(time_ca125 <= 188) %>%
  group_by(id) %>%
  filter(n() > 1) %>%
  ungroup()
ca125_all_6["id_num"]<- as.numeric(as.factor(ca125_all_6$id))
ca125_all_6["trialid_num"]<- as.numeric(as.factor(ca125_all_6$trialid))
ca125_all_6["log_ca125"]<-log(ca125_all_6$CA125_mes)

#COMPUATION
pak::pak("DenisRustand/INLAjoint")
library(INLAjoint)
library(INLA)

#Model input
ca125_all_6<-ca125_all_6%>% arrange(id_num, time_ca125)
ca125_all_6_surv<-ca125_all_6[!duplicated(ca125_all_6$id),]
internal_knots <- quantile(ca125_all_6$time_ca125, probs = c(0.5))
boundary_knots <- quantile(ca125_all_6$time_ca125, probs = c(0.18,0.98))
ns_basis <- ns(ca125_all_6$time_ca125,knots = internal_knots,Boundary.knots = boundary_knots)
f1 <- function(x) predict(ns_basis, x)[,1]
f2 <- function(x) predict(ns_basis, x)[,2]
n_trials<- length(unique(ca125_all_6$trialid_num))
library(parallel)
INLA::inla.setOption(num.threads = detectCores()-5)
#model setup
formLong <- log_ca125 ~ f1(time_ca125) + f2(time_ca125) + trt_bin:f1(time_ca125)+ trt_bin:f2(time_ca125)+ (1 + f1(time_ca125) + f2(time_ca125)|id_num)
formSurv <- INLA::inla.surv(pfs_time_d, pfs_statut) ~ trt_bin + (1|trialid_num)
model_setup <- joint(
  formSurv  = list(formSurv),
  formLong  = list(formLong),
  dataSurv  = as.data.frame(ca125_all_6_surv),
  dataLong  = as.data.frame(ca125_all_6),
  id        = "id_num",
  timeVar   = "time_ca125",
  family    = "gaussian",
  basRisk   = "rw2", 
  NbasRisk = 50,
  assoc     = "SRE_ind",
  control   = list(int.strategy = "eb", verbose=TRUE),
  run=FALSE, 
)
formLong_trial <- log_ca125 ~ trt_bin + f1(time_ca125) + f2(time_ca125) + trt_bin:f1(time_ca125)+ trt_bin:f2(time_ca125)+(1 + f1(time_ca125)+ f2(time_ca125)|trialid_num)
model_setup_trial <- joint(
  formSurv  = list(formSurv),
  formLong  = list(formLong_trial),
  dataSurv  = as.data.frame(ca125_all_6_surv),
  dataLong  = as.data.frame(ca125_all_6),
  id        = "trialid_num",
  timeVar   = "time_ca125",
  family    = "gaussian",
  basRisk   = "rw2", 
  NbasRisk = 50,
  assoc     = "SRE_ind",
  corRE= FALSE, 
  control   = list(int.strategy = "eb", verbose=TRUE), 
  run=FALSE,
)
#Treatment * splines random effects id
model_setup$.args$data$IDIntercept_S1 <- model_setup_trial$.args$data$IDIntercept_S1
model_setup$.args$data$WIntercept_S1 <- model_setup_trial$.args$data$WIntercept_S1*model_setup_trial$.args$data$trt_bin_S1
environment(model_setup$.args$formula)$IDIntercept_S1 <-model_setup$.args$data$IDIntercept_S1
environment(model_setup$.args$formula)$WIntercept_S1 <-model_setup$.args$data$WIntercept_S1
model_setup$.args$data$IDTrial0 <- model_setup_trial$.args$data$IDIntercept_L1
model_setup$.args$data$IDTrial1 <- model_setup_trial$.args$data$IDf1time_ca125_L1
model_setup$.args$data$IDTrial2 <- model_setup_trial$.args$data$IDf2time_ca125_L1

# Treatment * splines random effects weights
model_setup$.args$data$WTrial0 <- model_setup_trial$.args$data$WIntercept_L1
model_setup$.args$data$WTrial1 <- model_setup_trial$.args$data$Wf1time_ca125_L1 * model_setup_trial$.args$data$trt_bin_L1
model_setup$.args$data$WTrial2 <- model_setup_trial$.args$data$Wf2time_ca125_L1 * model_setup_trial$.args$data$trt_bin_L1

# Share and scale the trial level random effects in survival
model_setup$.args$data$IDTrial0_SRE <- model_setup_trial$.args$data$SRE_Intercept_L1_S1
model_setup$.args$data$IDTrial1_SRE <- model_setup_trial$.args$data$SRE_f1time_ca125_L1_S1
model_setup$.args$data$IDTrial2_SRE <- model_setup_trial$.args$data$SRE_f2time_ca125_L1_S1

# Priors
prior_prec <- list(prior = "loggamma",param = c(1.5, 0.5))

# Set environnement
form_env <- environment(model_setup$.args$formula)

form_env$IDTrial0 <- model_setup$.args$data$IDTrial0
form_env$IDTrial1 <- model_setup$.args$data$IDTrial1
form_env$IDTrial2 <- model_setup$.args$data$IDTrial2

form_env$WTrial0 <- model_setup$.args$data$WTrial0
form_env$WTrial1 <- model_setup$.args$data$WTrial1
form_env$WTrial2 <- model_setup$.args$data$WTrial2

form_env$IDTrial0_SRE <- model_setup$.args$data$IDTrial0_SRE
form_env$IDTrial1_SRE <- model_setup$.args$data$IDTrial1_SRE
form_env$IDTrial2_SRE <- model_setup$.args$data$IDTrial2_SRE

form_env$n_trials <- n_trials
form_env$prior_prec <- prior_prec

# Update formula
f_txt <- deparse1(model_setup$.args$formula)

f_txt <- sub(
  "n = 7933",
  "n = n_trials",
  f_txt,
  fixed = TRUE
)


model_setup$.args$formula <-
  as.formula(f_txt, env = form_env)


model_setup$.args$formula <- update(
  model_setup$.args$formula,
  ~ . +
    f(IDTrial0,WTrial0,model = "iid",n = n_trials,constr = FALSE,hyper = list(prec = prior_prec)) +
    f(IDTrial1,WTrial1,model = "iid",n = n_trials,constr = FALSE,hyper = list(prec = prior_prec)) +
    f(IDTrial2,WTrial2,model = "iid",n = n_trials,constr = FALSE,hyper = list(prec = prior_prec)) +
    f(IDTrial0_SRE,copy = "IDTrial0",hyper = list(beta = list(fixed = FALSE,param = c(0, 1),initial = 0.1))) +
    f(IDTrial1_SRE,copy = "IDTrial1",hyper = list(beta = list(fixed = FALSE,param = c(0, 1),initial = 0.1))) +
    f(IDTrial2_SRE,copy = "IDTrial2",hyper = list(beta = list(fixed = FALSE,param = c(0, 1),initial = 0.1)))
)

# Safety checks before running model
# Should be 13
stopifnot(length(unique(na.omit(model_setup$.args$data$IDIntercept_S1))) == n_trials)
# Formula environment should also contain 13 levels
stopifnot(length(unique(na.omit(environment(model_setup$.args$formula)$IDIntercept_S1))) == n_trials)
# Formula should no longer contain n = 7933
stopifnot(!grepl("n = 7933",deparse1(model_setup$.args$formula),fixed = TRUE))
# run the model
model_fit <- joint.run(model_setup, slientMode=TRUE)
#save model
saveRDS(model_fit, file = "jointed_inla_pfs.rds")
