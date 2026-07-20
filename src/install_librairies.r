#Intall packages
install_packages<-function(path_to_requirement){
    packages <- readLines(file.path(path_to_requirement, "requirements.txt"))
    packages <- packages[!grepl("^#", packages)]
    packages <- trimws(packages)
    packages <- packages[packages != ""]
    install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
    }
    }
    invisible(lapply(packages, install_if_missing))
    install.packages(c("survminer","markdown","performance"), type = "binary")
}
