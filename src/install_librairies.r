#Install packages
install.packages(
  "INLA",
  repos = c(
    INLA = "https://inla.r-inla-download.org/R/stable",
    CRAN = "https://cloud.r-project.org"
  ),
  dependencies = TRUE
)
pak::pak("DenisRustand/INLAjoint")
install_packages<-function(path_to_requirement){
    packages <- readLines(file.path(path_to_requirement, "requirements.txt"))
    packages <- packages[!grepl("^#", packages)]
    packages <- trimws(packages)
    packages <- packages[packages != ""]
    install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
    }
    }
    invisible(lapply(packages, install_if_missing))
    install.packages(c("survminer","markdown","performance"),  repos = "https://cloud.r-project.org", type = "binary")
}
