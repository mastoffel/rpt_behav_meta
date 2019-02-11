# analysis of data in papers to obtain repeatabilities

library(tidyverse)
library(readxl)
library(XLConnect)
library(httr)
library(rptR)

# variables from the meta table
meta_table_template <- data.frame("Key" = NA,                # identifier in meta-table "data/meta_table_filled.xlsx"
                         "species_common" = NA, 
                         "sample_size" = NA, 
                         "sex" = NA,                # 0 = both, 1 = females, 2 = males
                         "behaviour" = NA,          # measured behaviour as stated by authors
                         "context" = NA,            # 1 = lab exp. / lab-reared, 2 = lab exp. / wild-caught, 3 field exp
                         "type_of_treatment"= NA,  # 0 = no treatment, 1 = between-subject treatment, 2 = within-subject
                         "treatment"= NA,          # Verbal description
                         "life_stage"= NA,         # "juvenile", "adult", "both"
                         "event"= NA,              # Major life-event, such as metamorphosis between measurements
                         "R"= NA,
                         "R_se"= NA,
                         "CI_lower"= NA,
                         "CI_upper"= NA,
                         "p_val"= NA,
                         "t1"= NA,                  # timepoint of first measurement in days old (or mean of measurements)
                         "t2"= NA,                  # timepoint of second measurement in days old# (or mean of measurements)
                         "delta_t"= NA,             # difference between timepoints in days
                         "remarks"= NA)


##### Paper 1 ##########
# Baker, Matthew R.; Goodman, Alex; C., er; Santo, Jonathan B.; Wong, Ryan Y. (2018)
# Repeatability and reliability of exploratory behavior in proactive and reactive zebrafish, Danio rerio
# 
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-018-30630-3/MediaObjects/41598_2018_30630_MOESM1_ESM.xlsx"

# download xlsx file
tmp <- tempfile(fileext = ".xlsx")
httr::GET(url = url,
    write_disk( tmp) )
# read data table
dat <- read_excel(tmp) # skip = n

# coding used in the table:
# Strain (1 = LSB 2 = HSB 3 = WC)
# Sex (1 = Male 2 = Female)
# Standard Length (cm)
# speed in (cm/s) / time in (s)

# reformat data into longformat for rptR using data.table
require(data.table)
dat_tidy <- as_tibble(melt(setDT(dat), id=c(1,2,3,19), measure=patterns("^Swim", "^Stationary", "^Time"), 
    value.name=c("swim_speed", "stationary_time", "time_in_center"), variable.name="week"))
names(dat_tidy)[c(2,3,4)] <- c("strain", "sex", "length")

# select only wild-caught strain as treatment groups are selected for or against stationary behaviour
dat_wc <- dat_tidy %>% filter(strain == 3) %>% 
            mutate(week = as.numeric(as.character(week))) %>% 
            mutate(day = 1 + (week - 1) * 7)

# calculate all possible pairs of repeatabilities for the 5 measurements
combinations <- as.data.frame(t(combn(names(table(dat_wc$day)), m = 2)))

get_rpts <- function(timepoints, dep_var) {
    # filter only two timepoints to calculate rpt over
    dat_filtered <- dat_wc %>% 
                        filter(day %in% as.numeric(!!timepoints))
    rpt_formula <- formula(paste0(dep_var, '~ (1|ID)'))
    # calculate repeatability
    rpt_out <- rptGaussian( rpt_formula, grname = c("ID"), data = dat_filtered, nboot = 1000)
    # create df ready to be in the meta-table
    out_df <- data.frame(R = unlist(rpt_out$R), R_se = unlist(rpt_out$se), 
                         t1 = timepoints[1], t2 = timepoints[2], row.names = NULL,
                         sample_size = length(table(dat_filtered$ID)), p_val = unlist(rpt_out$P$LRT_P),
                         behaviour = dep_var)
}

# combine repeatabilities for different variables into one
all_rpts_var1 <- bind_rows(apply(combinations, 1, get_rpts, "swim_speed"))
all_rpts_var2 <- bind_rows(apply(combinations, 1, get_rpts, "stationary_time"))
all_rpts_var3 <- bind_rows(apply(combinations, 1, get_rpts, "time_in_center"))

all_rpts <- rbind(all_rpts_var1, all_rpts_var2, all_rpts_var3)

meta_table <- data.frame(Key = "9Z5A7EQB", species_common = "zebrafish",
     all_rpts,  
    "sex" = 0,                # 0 = both, 1 = females, 2 = males
    "context" = 2,            # 1 = lab exp. / lab-reared, 2 = lab exp. / wild-caught, 3 field exp
    "type_of_treatment"= 0,  # 0 = no treatment, 1 = between-subject treatment, 2 = within-subject
    "treatment"= NA,          # Verbal description
    "life_stage"= "adult",         # "juvenile", "adult", "both"
    "event"= NA,              # Major life-event, such as metamorphosis between measurements
    "CI_lower"= NA,
    "CI_upper"= NA,
               # timepoint of second measurement in days old# (or mean of measurements)
    "delta_t"= NA,             # difference between timepoints in days
    "remarks"= "Authors declare everything as exploratory behaviour") %>% 
    mutate(t1 = as.numeric(t1)) %>% 
    mutate(t2 = as.numeric(t2)) %>% 
    mutate(delta_t = sqrt((t2 - t1)^2))

# decided to always write first two authors and year
meta_table <- select(meta_table, names(meta_table_template))
write_delim(meta_table, path = "output/Baker_Goodman_2018.txt", delim = " ", col_names = TRUE)
# meta_table <- read_delim(file = "output/Baker_Goodman_2018.txt", delim = " ", col_names = TRUE)



##### Paper 2, incomplete #######
# Fisher, David N.; James, Adele; Rodriguez-Munoz, Rol; , o; Tregenza, Tom (2015)
# behaviour in captivity predicts some aspects of natural behaviour, but not others, in a wild cricket population
# FVQNEVJS

# data set: bushcricket personality in the wild
url <- "https://ore.exeter.ac.uk/repository/bitstream/handle/10871/16930/bbivar.txt?sequence=6&isAllowed=y"
dat <- read_delim(url, delim = "\t")

url2 <- "https://ore.exeter.ac.uk/repository/bitstream/handle/10871/16930/ae.l-w.txt?sequence=2&isAllowed=y"
dat <- read_delim(url2, delim = "\t")
# how many individuals?
length(table(dat$tag))


###### Paper 3, incomplete #######
### Bierbach, David; Laskowski, Kate L.; Wolf, Max
### behavioural individuality in clonal fish arises despite near-identical rearing conditions
### MV5GA8PC

url <- "https://datadryad.org/bitstream/handle/10255/dryad.140031/Bierbach%20et%20al%20clonal%20molly%20behav%20development_data%20for%20deposit.xlsx?sequence=1"

# download xlsx file
tmp <- tempfile(fileext = ".xlsx")
httr::GET(url = url,
    write_disk( tmp) )
# read data table
dat <- read_excel(tmp) # skip = n
dat


####### Paper 4 ##########
# additional calculations to fill meta table, sorted by study

library(lubridate)

# Aplin, L. M.; Firth, J. A.; Farine, D. R.; Voelkl, B.; Crates, R. A.; Culina, A.; Garroway, C. J.; Hinde, C. A.; Kidd, L. R.; Psorakis, I.; Milligan, N. D.; Radersma, R.; Verhelst, B. L.; Sheldon, B. C.
# consistent individual differences in the social phenotypes of wild great tits, parus major
# Key: JQJKK3Y6
# DOI: 10.1016/j.anbehav.2015.07.016

# Data collection
# Winter 1 (2011): 27 February
winter1_start <- dmy(03122011)
winter1_end <- dmy(27022012)
time_diff_winter1 <- as.numeric(winter1_end - winter1_start)
average_time_diff_winter1 <- time_diff_winter1/ 11 # see paper, median of 11 measurements
# 7.8

winter2_start <- dmy(01122012)
winter2_end <- dmy(03032013)
time_diff_winter2 <- as.numeric(winter2_end - winter2_start) 
average_time_diff_winter2 <- time_diff_winter2/ 10


winter3_start <- dmy(30112013)
winter3_end <- dmy(02032014)
time_diff_winter3 <- as.numeric(winter3_end - winter3_start)
average_time_diff_winter3 <- time_diff_winter3 / 12

# repeatability across all 3 winters (~ 365*2 = 730 days in between)
# 565 tits observed in 2 winters, 210 in 3 winters (--> maybe take mean??)

# extract errorbars data from figure 2 to get se's for within year repeatability
library(metaDigitise)
?metaDigitise
data <- metaDigitise(dir = "to_digitise/")

# question: sample size is number of individuals or number of measurements?

# figure out what t1, t2 delta_t are when age of birds not known (we know its some juveniles, some adults)
# se extracted from figure (see below)
meta_table <- tribble(
    ~behaviour,              ~R,  ~R_se, ~sample_size, ~t1, ~t2, ~delta_t, ~remarks,
    "degree",                0.46, 0.01,      1053,     NA,  NA, time_diff_winter1, "From start to end of season",  
    "group_size",            0.43, 0.01,      1053,     NA,  NA, time_diff_winter1, "From start to end of season",  
    "association_strength",  0.41, 0.01,      1053,     NA,  NA, time_diff_winter1, "From start to end of season",  
    "betweenness",           0.19, 0.01,      1053,     NA,  NA, time_diff_winter1, "From start to end of season",  
    
    "degree",                0.61, 0.02,      729,      NA,  NA, time_diff_winter2, "From start to end of season",  
    "group_size",            0.64, 0.01,      729,      NA,  NA, time_diff_winter2, "From start to end of season",  
    "association_strength",  0.64, 0.01,      729,      NA,  NA, time_diff_winter2, "From start to end of season",  
    "betweenness",           0.18, 0.01,      729,      NA,  NA, time_diff_winter2, "From start to end of season",  
    
    "degree",                0.58, 0.02,      816,      NA,  NA, time_diff_winter3, "From start to end of season",  
    "group_size",            0.60, 0.01,      816,      NA,  NA, time_diff_winter3, "From start to end of season",  
    "association_strength",  0.63, 0.003,     816,      NA,  NA, time_diff_winter3, "From start to end of season",  
    "betweenness",           0.38, 0.02,      816,      NA,  NA, time_diff_winter3, "From start to end of season",  
    
    "degree",                0.55, 0.005,     210,      NA,  NA, 730,               "Over three years, only minor differences in age classes and sexes", 
    "group_size",            0.51, 0.005,     210,      NA,  NA, 730,               "Over three years, only minor differences in age classes and sexes", 
    "association_strength",  0.57, 0.00175,   210,      NA,  NA, 730,               "Over three years, only minor differences in age classes and sexes", 
    "betweenness",           0.33, 0.02,      210,      NA,  NA, 730,               "Over three years, only minor differences in age classes and sexes"
)

meta_table <- meta_table %>% 
                mutate(Key = "JQJKK3Y6",
                       species_common = "Great tit",
                       sex = 0,
                       context = NA,
                       type_of_treatment = NA,
                       treatment = NA,
                       life_stage = "both",
                       event = NA,
                       CI_lower = NA,
                       CI_upper = NA,
                       p_val = NA)

# order columns by meta_table_template
meta_table <- select(meta_table, names(meta_table_template))
meta_table
write_delim(meta_table, path = "output/Aplin_Firth_2015.txt", delim = " ", col_names = TRUE)
