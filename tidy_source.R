
if (!require("pacman")) install.packages("pacman")

# Declare CRAN dependencies :
r_dependencies <- c("here", "formatR")

# load dependencies via pacman
pacman::p_load(r_dependencies, character.only = TRUE)

source_dir <-  here("profile_binr/_R")
tidy_dir(source_dir)

