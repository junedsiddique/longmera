## code to prepare gruder dataset

setwd("C:/Users/jsi894/OneDrive - Northwestern University/_FromNUBox/longmera/data-raw")
gruder <- read.table(file="SmkStudy_Covariates.dat", quote="\"", comment.char="")
names(gruder) <- c("id", "quit", "time", "group", "grouptime", "race", "tv", "manual")

# quit = 0 if smoking, 1 if abstinent
# time = 0, 1, 2, 4 for post-intervention, 6-, 12-, and 24-month follow-ups
# grp  = 0 for control 1 for intervention
# grptime = grp x time
# racew = 0 for non-white, 1 for white
# tv = 0 for not watching tv, 1 for watching
# manual = 0 for not reading manual, 1 for reading manual

usethis::use_data(gruder, overwrite = TRUE)

