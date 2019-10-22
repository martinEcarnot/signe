library(tourr)

if(Sys.getenv("RSTUDIO") == "1" & # check if running in RStudio
   .Platform$OS.type == "unix") quartz() else X11()
library(tourr)
animate_xy(flea[, 1:6])
