#Intall packages
install_packages<-function(path_to_requirement){
    packages <- readLines(paste(path_to_requirement, "requirements.txt"))
    install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
    }
    }
    invisible(lapply(packages, install_if_missing))
    install.packages(c("survminer","markdown","performance"), type = "binary")
}
# Required packages
#Handling data
library(dplyr)
library(tidyverse)
library(tidyr)
library(writexl)
library(openxlsx)
library(reshape2)
library(lubridate)
library(Dict)
library(broom.mixed)
#Plots
library(ggplot2)
library(gridExtra)
library(patchwork)
#CA125 modeling
library(splines)
library(statmod)
library(msm)
library(Matrix)
library(pracma)
library(lme4)
library(performance)
library(purrr)
library(lmerTest)
library(Metrics)
#Survival modeling
library(survival)
library(survminer)
library(flexsurv)
library(rstpm2)
#Joint modeling
library(meta)
library(gtools)
library(lcmm)
library(INLA)
library(INLAjoint)