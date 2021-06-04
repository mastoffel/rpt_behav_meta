# analysis of data in papers to obtain repeatabilities

library(tidyverse)
library(readxl)
library(XLConnect)
library(httr)
library(rptR)
library(tabulizer)
library(magrittr)
library(broom)
library(tidyverse)
library(metaDigitise)
library(AnAgeScrapeR)
library(lubridate)
library(digitize)
library(datapasta)

# variables from the meta table
meta_table_template <- tibble("Key" = NA,               # identifier in meta-table "data/meta_table_filled.xlsx"
                         "species_common" = NA, 
                         "species_latin" = NA,
                         "sample_size" = NA, 
                         "measurements_per_ind" = NA,   # new, check papers 1-6 again
                         "sex" = NA,                    # 0 = both, 1 = females, 2 = males
                         "behaviour" = NA,              # measured behaviour as stated by authors
                         "context" = NA,                # 1 = lab exp. / lab-reared, 2 = lab exp. / wild-caught, 3 field exp / maybe another category: 4 field behaviour?
                         "type_of_treatment"= NA,       # 0 = no treatment, 1 = between-subject treatment, 2 = within-subject
                         "treatment"= NA,               # Verbal description
                         "life_stage"= NA,              # "juvenile", "adult", "both"
                         "event"= NA,                   # Major life-event, such as metamorphosis between measurements
                         "R"= NA,
                         "R_se"= NA,
                         "CI_lower"= NA,
                         "CI_upper"= NA,
                         "p_val"= NA,
                         "t1"= NA,                      # timepoint of first measurement in days old (or mean of measurements)
                         "t2"= NA,                      # timepoint of second measurement in days old# (or mean of measurements)
                         "delta_t"= NA,                 # difference between timepoints in days
                         "remarks"= NA,
                         "max_lifespan_days" = NA)


###### Paper 1: Baker_Goodman 2018 ##########
# Baker, Matthew R.; Goodman, Alex; C., er; Santo, Jonathan B.; Wong, Ryan Y. (2018)
# Repeatability and reliability of exploratory behavior in proactive and reactive zebrafish, Danio rerio
# 9Z5A7EQB
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
# write_delim(meta_table, path = "output/Baker_Goodman_2018.txt", delim = " ", col_names = TRUE)
# meta_table <- read_delim(file = "output/Baker_Goodman_2018.txt", delim = " ", col_names = TRUE)

max_life <- scrape_AnAge(latin_name = c("Danio rerio"), vars = "maximum_longevity_yrs", download_data = FALSE)
max_lifespan <- as.numeric(max_life$maximum_longevity_yrs) * 365

# add on
meta_table_new <- read_delim("data/meta_tables_to_be_updated/Baker_Goodman_2018_old.txt", delim = " ") %>% 
        mutate(species_latin = "Danio_rerio",
               max_lifespan_days = max_lifespan,
               t1 = t1 + 21*30,
               t2 = t2 + 21*30,
               measurements_per_ind = 2)

write_delim(meta_table_new, path = "output/Baker_Goodman_2018.txt", delim = " ", col_names = TRUE)






###### Paper 2: Fisher_Adele 2015 incomplete #######
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


###### Paper 4: Aplin_Firth 2015 ##########
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
data <- metaDigitise(dir = "data/to_digitise/")

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
                       context = 3,
                       type_of_treatment = 0,
                       treatment = NA,
                       life_stage = "both",
                       event = NA,
                       CI_lower = NA,
                       CI_upper = NA,
                       p_val = NA)

# order columns by meta_table_template
meta_table <- select(meta_table, names(meta_table_template))
meta_table
# write_delim(meta_table, path = "output/Aplin_Firth_2015.txt", delim = " ", col_names = TRUE)

# update meta table
dat_tit <- scrape_AnAge(latin_name = "Parus major", vars = "maximum_longevity_yrs", download_data = FALSE)
tit_lifespan <- as.numeric(dat_tit$maximum_longevity_yrs) * 365
avg_age_adult <- tit_lifespan*0.25

# measurements_per_ind 
mes_per_ind <- (9835 + 6853 + 7940) / (1053 + 729 + 816)

meta_table_new <- read_delim("data/meta_tables_to_be_updated/Aplin_Firth_2015_old.txt", delim = " ") %>% 
                        mutate(t1 = avg_age_adult,
                               t2 = avg_age_adult + delta_t,
                               species_latin = "patus_major",
                               measurements_per_ind = mes_per_ind,
                               max_lifespan_days = tit_lifespan)

write_delim(meta_table_new, path = "output/Aplin_Firth_2015.txt", delim = " ", col_names = TRUE)

###### Paper 5, Arroyo_Mougeot 2017 #######
# Arroyo, Beatriz; Mougeot, Francois; Bretagnolle, Vincent
# individual variation in behavioural responsiveness to humans leads to differences in breeding success and long-term population phenotypic changes
# Key: V5R8XCKE
# DOI: 10.1111/ele.12729

# a bit tricky: Sample sizes are a bit unclear (number of measurements are reported only, see table)
# between years consists of 30 individuals with 2-7 years observation
# as only range given, simulate potential poisson
num_obs <- (rpois(n = 10000, lambda = 3 ))
delta_years_in_days <- (mean(num_obs[num_obs >= 2]) - 1) * 365

# within year delta (breeding season) ~ 75 days
# https://onlinelibrary.wiley.com/doi/full/10.1046/j.0019-1019.2001.00009.x

# total: 833 nests, 2402 visity
# number of visits per nest: 3.1, n = 98 


harrier_data <- scrape_AnAge("Circus pygargus", vars = "maximum_longevity_yrs", download_data = FALSE)
harrier_ls_max <- as.numeric(harrier_data$maximum_longevity_yrs) * 365
avg_adult_age <- harrier_ls_max * 0.25

tribble(
    ~behaviour,                   ~R,  ~R_se,        ~Nobs, ~t1,               ~t2,                             ~delta_t,                           ~remarks, ~sex,
    "nest_departure_distance", 0.001,  0.023,         1662,  avg_adult_age,  avg_adult_age + 75,                     75, 'within-year-repeatability', 1,
    "fleeing probability",     0.541,  0.049,         2399,  avg_adult_age,  avg_adult_age + 75,                     75, 'within-year-repeatability', 1,    
    "passive_if_present",      0.520,  0.059,         1844,  avg_adult_age,  avg_adult_age + 75,                     75, 'within-year-repeatability', 1, 
    "defence_intensity_PCA",   0.518,  0.033,         1014,  avg_adult_age,  avg_adult_age + 75,                     75, 'within-year-repeatability', 1,        
    "nest_departure_distance", 0.001,  0.017,          197,  avg_adult_age,  avg_adult_age +delta_years_in_days,      delta_years_in_days, 'between-year-repeatability (range 2-7 years)', 1,
    "fleeing probability",     0.341,  0.135,          299,  avg_adult_age,  avg_adult_age +delta_years_in_days,      delta_years_in_days, 'between-year-repeatability (range 2-7 years)', 1,    
    "passive_if_present",      0.522,  0.194,          252,  avg_adult_age,  avg_adult_age +delta_years_in_days,      delta_years_in_days, 'between-year-repeatability (range 2-7 years)', 1, 
    "defence_intensity_PCA",   0.131,  0.073,          189,  avg_adult_age,  avg_adult_age +delta_years_in_days,      delta_years_in_days, 'between-year-repeatability (range 2-7 years)', 1
) %>% 
    mutate(Key = "V5R8XCKE",
           species_common = "montagus_harrier",
           species_latin = "Circus_pygargus",
           sample_size = Nobs / 3.1,
           measurements_per_ind = 3.1,
           context = 3,
           type_of_treatment = NA,
           treatment = NA,
           life_stage = "adult",
           event = remarks,
           CI_lower = NA,
           CI_upper = NA,
          p_val = NA,
          max_lifespan_days = harrier_ls_max) %>% 
    select(-Nobs) -> meta_table
    
write_delim(meta_table, path = "output/Arroyo_Mougeot 2017 .txt", delim = " ", col_names = TRUE)
    


###### Paper 6: Ballew_Mittelbach 2017 ########
# Ballew, Nicholas G.; Mittelbach, Gary G.; Scribner, Kim T. 2017	
# fitness consequences of boldness in juvenile and adult largemouth bass
# RR968ED6		

# note: juvenile: 1-3 years, afterwards adult
meta_table <- tribble(
    ~behaviour,    ~R,  ~R_se,  ~CI_lower, ~CI_upper, ~sample_size,     ~t1,     ~t2, ~delta_t, ~life_stage, ~sex, ~remarks, 
    "boldness",  0.54,     NA,       0.38,      0.67,           93, 36*30.5, 48*30.5,       NA,   "both",       0, "Boldness was a PC from 4 different behaviours",    
    "boldness",  0.75,     NA,       0.63,      0.84,           71, 48*30.5, 60*30.5,       NA,   "adult",      0, "Boldness was a PC from 4 different behaviours", 
    "boldness",  0.50,     NA,       0.30,      0.65,           71, 36*30.5, 60*30.5,       NA,   "both",       0, "Boldness was a PC from 4 different behaviours"   
)
    
bass_data <- scrape_AnAge(latin_name = "Micropterus salmoides", vars = "maximum_longevity_yrs", download_data = FALSE)
max_bass <- as.numeric(bass_data$maximum_longevity_yrs) * 365

meta_table <- meta_table %>% 
    mutate(Key = "RR968ED6",
           species_common = "largemouth_bass",
           species_latin = "Micropterus_salmoides",
           context = 2,
           type_of_treatment = 0,
           treatment = NA,
           p_val = NA,
           event = NA, 
           measurements_per_ind = 2,
           max_lifespan_days = max_bass)

# order columns by meta_table_template
meta_table <- select(meta_table, names(meta_table_template))
meta_table
write_delim(meta_table, path = "output/Ballew_Mittelbach_2017.txt", delim = " ", col_names = TRUE)




    
###### Paper 7: Bosco et al 2017, nearly done: waiting for authors to write about ages #######
# Bosco, Jennifer M.; Riechert, Susan E.; O'Meara, Brian C.	2017
# the ontogeny of personality traits in the desert funnel-web spider, agelenopsis lisa (araneae: agelenidae)

library(tabulizer)
library(magrittr)
library(broom)
library(tidyverse)
location <- "data/papers/Bosco et al. - 2017 - the ontogeny of personality traits in the desert f.pdf"
# Extract the table
out <- extract_tables(location, pages = 6)
out
locate_areas(location, pages = 6)    

R_table <- extract_tables(location,
    output = "data.frame",
    pages = c(6, 6), # include pages twice to extract two tables per page
    area = list(
        c(82.72286, 43.16004, 643.55584, 295.53488 ),
        c(58.88746, 305.34945, 256.58109, 547.90972)
    
    ),
    guess = FALSE
)

str(R_table[1])
R_table[[1]] <- mutate(R_table[[1]], p = as.numeric(p))
meta_table <- R_table %>% 
                reduce(bind_rows) %>% 
                as_tibble() %>% 
                rename(behaviour = Trait,
                       life_stage = Life.stage,
                       R = Radj, 
                       CI = `CI..95..`, 
                       p_val = p) %>% 
                separate(CI, into = c("CI_lower", "CI_upper"), sep = ",\\s*") %>% 
                mutate(life_stage = str_replace(life_stage, " \\+ ", "_")) %>% 
                separate(life_stage, into = c("life_stage", "sex"), sep = " ") %>% 
                dplyr::filter(life_stage != "All") %>%  # filter all, otherwise pseudoreplication
                filter(life_stage != "") %>% 
                filter(life_stage != "female") %>% 
                mutate(sex = replace_na(sex, "female")) %>% 
                mutate(sex = case_when(sex == "male" ~ 2,
                                       sex == "female" ~ 1)) %>% 
                #seperate(life_stage, into = c("life_stage", ""))
                print(n = 60)

# info about timings:
# We subjected each of 98 of 120 surviving spiders to a battery of 11 behavioral trait tests at three stages in the life cycle: juvenile (third–
# fifth instar), penultimate (one molt removed from sexual maturity), and sexually mature. Two weeks separated replicate within stage tests.
# We completed the first trial within a particular life stage seven days following an individual’s molt to a new stage and three days after an ad libitum feeding.

# life time of Agelenopsis lisa: probably a year with 10 instars
# rough calculation assuming birth in spring, mature late summer:
time_one_molt <- (365 * 3/4) / 10
juvenile_age <- time_one_molt * 4
penultimate_age <- time_one_molt * 9
adult_age <- (365 * 3/4) + (365 - (365 * 3/4))/2

meta_table <- meta_table %>% 
    mutate(Key = "V5SYT8AL",
           species_common = "desert_funnel_web_spider",
           sample_size = 98,
            # see paper
           measurements_per_ind = rep(c(19, 54, 22, 67, 22, 54, 19, 54, 22, 54), 4),
           behaviour = c(rep("aggressiveness_or_foraging_behaviour", 10), rep("exploration_or_activity", 10),
                         rep("latency_exploring_new_envir", 10), rep("latency_to_return_after_predatory_cue", 10)),
           context = 1,
           type_of_treatment = 0,
           treatment = NA,
           R_se = NA,
           event = rep(c(rep(NA, 6), "molts", "molts", "molts_sexual_maturity", "molts_sexual_maturity"), 4),
           t1 = rep(c(rep(juvenile_age, 2), rep(penultimate_age, 2), rep(adult_age, 2), rep(juvenile_age, 2), rep(penultimate_age, 2)), 4),
           # 14 days between measurements within stages)
           t2 = rep(c(rep(juvenile_age, 2) + 14, rep(penultimate_age, 2) + 14, rep(adult_age, 2) + 14, 
                    rep(penultimate_age, 2), rep(adult_age, 2)), 4)) %>% 
    mutate(delta_t = t2 - t1,
           remarks = "timings in days assumed from standard life-cycle of Agelenopsis")
            
#write_delim(meta_table, path = "output/Bosco_Riechert_2016.txt", delim = " ", col_names = TRUE)


scrape_AnAge("agelenopsis lisa", download_data = FALSE, vars = "maximum_longevity_yrs")

# lifespan roughly a year (see paper discussion)
meta_table_new <- read_delim("data/meta_tables_to_be_updated/Bosco_Riechert_2016_old.txt", delim = " ") %>% 
                        mutate(species_latin = "Micropterus_salmoides",
                            max_lifespan_days = 365)

write_delim(meta_table_new, path = "output/Bosco_Riechert_2016.txt", delim = " ", col_names = TRUE)



###### Paper 9: Bubac_Coltman et 2018 #######

# Bubac, Christine M.; Coltman, David W.; Bowen, W. Don; Lidgard, Damian C.; Lang, Shelley L. C.; den Heyer, Cornelia E.	2018
# AKG6LAL3
# repeatability and reproductive consequences of boldness in female gray seals	BEHAVIORAL ECOLOGY AND SOCIOBIOLOGY

# Boldness was measured for a total of 469 branded females
# during nine consecutive breeding seasons (2008 to 2016).
# Over a 9-year period, 2504 observations were made for between year repeatability
# measures. Boldness scores were obtained within a single year on 35 females between 2009 and 2016 (105 total
# observations). On average, females were tested 5.5 ± 0.09
# SE and 3.0 ± 0.30 SE times per female for between and within
# year sampling efforts, respectively.

# Only females with at least two observations were included in between (n = 458 females) and within year (n = 35 females)
# repeatability analyses.

# adult female grey seals were 6 to 43 years old, so mean age is 18.5 years or 6753 days
# data was collected over a 9 year period with an average of 5.5 measurments per female, so over a minimum of 5.5 years or 2008 days (t2 = 6753 + 2008)

max_age <- scrape_AnAge("Halichoerus grypus", vars = "maximum_longevity_yrs", download_data = FALSE)
max_age_greyseal <- as.numeric(max_age$maximum_longevity_yrs) * 365

meta_table <- tribble(
    ~behaviour,              ~R,  ~R_se, ~CI_lower, ~CI_upper, ~sample_size,  ~t1,   ~t2, ~delta_t, ~measurements_per_ind, ~remarks,
    "boldness",            0.61,     NA,      0.57,      0.66,          458, 6753,   8760,    2008,                    5.5, "one measurement per year and average of 5.5 equals 5.5 years span of rpt",
    "boldness",            0.82,     NA,      0.54,      0.92,           35, 6753,  7118,      365,                      3, "measurement times start based on mean adult life time" 
)

meta_table <- meta_table %>% 
    mutate(Key = "AKG6LAL3",
           species_common = "grey_seal",
           species_latin = "Halichoerus grypus",
           sex = 1,
           context = 3,
           type_of_treatment = NA,
           treatment = NA,
           life_stage = "adult",
           event = NA,
           p_val = NA,
        max_lifespan_days = max_age_greyseal)

write_delim(meta_table, path = "output/Bubac_Coltman_2018.txt", delim = " ", col_names = TRUE)







###### Paper 10: Burtka_Grindstaff 2013 ######
# Burtka, Jennifer L.; Grindstaff, Jennifer L.	2013
# LFNHZNJS		
# repeatable nest defense behavior in a wild population of eastern bluebirds (sialia sialis) as evidence of personality

location <- "data/papers/Burtka and Grindstaff - 2013 - repeatable nest defense behavior in a wild populat.pdf"
# Extract the table
loc_area <- locate_areas(location, pages = 7)    

R_table <- extract_tables(location,
    output = "data.frame",
    pages = c(7), # include pages twice to extract two tables per page
    area = list(
        unlist(loc_area)
    ),
    guess = FALSE
)

# Eastern Bluebirds live around 6-10 years, average age assumed: 4 years = 1460 days
# within year season is 230 - 90 = 140 days. (max)


meta_table <- R_table[[1]] %>% 
                select(2,5, 6, 8) %>% 
                rename(R = `τ`,
                       sample_size = N,
                       R_se = SE,
                       p_val = p) %>% 
                mutate(t1 = c(rep(1460, 6)),
                       t2 = c(rep(1460+140, 4), rep(1460 + 365, 2)),
                       sex = rep(c(2,1), 3),
                       behaviour = "nest_defence_aggressiveness",
                    # see table 3
                       measurements_per_ind = c(2.2, 2.3, 2.1, 2.3, 2, 2),
                       context = 3,
                       type_of_treatment = 0,
                       treatment = NA,
                       life_stage = "adult",
                       event = c(rep("within_years", 4), rep("between_years", 2)),
                       CI_lower = NA,
                       CI_upper = NA,
                       delta_t = c(rep(140, 4), rep(365, 2)),
                       remarks = NA,
                       Key = "LFNHZNJS",
                       species_common = "eastern_bluebird")


# write_delim(meta_table, path = "output/Burtka_Grindstaff_2013.txt", delim = " ", col_names = TRUE) 


# lifespan
max_lifespan_bluebird <- as.numeric(scrape_AnAge("sialia sialis", vars = "maximum_longevity_yrs", download_data = FALSE)$maximum_longevity_yrs) * 365

meta_table_new <- read_delim("data/meta_tables_to_be_updated/Burtka_Grindstaff_2013_old.txt", delim = " ") %>% 
    mutate(species_latin = "sialia_sialis",
        max_lifespan_days = max_lifespan_bluebird)

write_delim(meta_table_new, path = "output/Burtka_Grindstaff_2013.txt", delim = " ", col_names = TRUE)


# check again from here 

###### Paper 11: Cabrera_Andres_2017 #######
 #Cabrera, Doreen; Andres, Daniel; McLoughlin, Philip D.; Debeffe, Lucie; Medill, Sarah A.; Wilson, Alastair J.; Poissant, Jocelyn	2017	
# ESF6PD88	
# island tameness and the repeatability of flight initiation distance in a large herbivore

location <- "data/papers/Cabrera et al. - 2017 - island tameness and the repeatability of flight in-rotated.pdf"

# average time of measurements among days: 
delta_t_among <- 10.2
# all data from foles (maximum of one year old, but rather a few days to a few months)
# few days to few month old, so take 3 month
age_fole <- 90
# measurements_per_ind
378/105 # 3.6

max_ls_data <- scrape_AnAge("Equus caballus", vars = "maximum_longevity_yrs", download_data = FALSE)
max_ls_horse <- as.numeric(max_ls_data$maximum_longevity_yrs) * 365

tribble(
    ~timespan,         ~sample_size, ~measurements_per_ind,    ~R, ~R_se,     ~t1,                           ~t2,      ~delta_t,
  #  "within_among_days",          103,               376/103,   0.42,  0.06,  age_fole, age_fole + delta_t_among, delta_t_among,
    "within_days",                156,               375/156,   0.55,  0.05,  age_fole,           age_fole + 0.5,           0.5,
    "among_days",                  45,                 99/45,   0.39,  0.12,  age_fole, age_fole + delta_t_among, delta_t_among
) %>% 
    mutate(Key = "ESF6PD88",
           species_common = "horse",
           species_latin = "Equus caballus",
           sex = 0,
           behaviour = "flight_initiation_distance",
           context = 3,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "juvenile",
           event = NA,
           CI_lower = NA,
           CI_upper = NA,
           p_val = NA,
           remarks = "foles under 1 year",
           max_lifespan_days = max_ls_horse) %>% 
    select(-timespan) -> meta_table

write_delim(meta_table, path = "output/Cabrera_Andres_2017.txt", delim = " ", col_names = TRUE)


###### Paper 12: D'Eath_2004 #######
# D'Eath, RB 2004	
# 2WGN2QQC		
# consistency of aggressive temperament in domestic pigs: the effects of social experience and social disruption

pig_ls <- scrape_AnAge("Sus scrofa", vars = "maximum_longevity_yrs", download_data = FALSE)
pig_max_ls <- as.numeric(pig_ls$maximum_longevity_yrs) * 365

# pigs were tested with two different intruders on the same day
tribble(
     ~sample_size, ~measurements_per_ind,    ~R, ~R_se, ~p_val,    ~t1,   ~t2,      ~delta_t,
          112,                         2,   0.303,  NA,  0.002,     48,  48.1,       0.1,
          112,                         2,   0.453,  NA,  0.001,     80,  80.1,       0.1,
          112,                         2,   0.629,  NA,  0.001,    113, 113.1,       0.1,
          112,                         6,   0.413,  NA,  NA,        48,     113,      65
) %>% 
    mutate(Key = "2WGN2QQC",
           species_common = "domestic_pig",
           species_latin = "Sus_scrofa",
           sex = 0,
           behaviour = "aggressiveness",
           context = 1,
           type_of_treatment = c(0,0,0,2),
           treatment = c(NA, NA, NA, "litter_change"),
           life_stage = c("juvenile", "adult", "adult", "both"),
           event = NA,
           CI_lower = NA,
           CI_upper = NA,
           remarks = c("correlation", "correlation", "correlation", "repeatability"),
           max_lifespan_days = pig_max_ls
           ) -> meta_table

write_delim(meta_table, file = "output/D'eath_2004.txt", delim = " ", col_names = TRUE)

###### Paper 13: David_Auclair_2012 ######
# David, Morgan; Auclair, Yannick; Cezilly, Frank	2012	
# assessing short- and long-term repeatability and stability of personality in captive zebra finches using longitudinal data
# Z9FL9L9P	

# long term  exploration delta_t = 263 d
# long term struggling rate detla_t = 209 d
# short_term = 7 day first session, 3 day second session

# zebra finch life span : 2-3 years (Zann 1996)
# e.g. long term intervall ~23 % of the life of a wild zebrafinch
# short term

# only females, two month of age
# 20 for explor, 17 for struggling

zf <- scrape_AnAge("Taeniopygia guttata", "maximum_longevity_yrs", download_data = FALSE)
zf_lv <- as.numeric(zf$maximum_longevity_yrs) * 365

tribble(
    ~behaviour,   ~sample_size, ~measurements_per_ind,    ~R, ~CI_lower, ~CI_upper, ~p_val,    ~t1,   ~t2,      ~delta_t,
    "exploration",          17,                     2,    0.81,    0.59,     0.95,  0.001,     61,    68,            7, 
    "exploration",          17,                     2,    0.35,       0,     0.73,  0.14,   322.5, 325.5,            3,   
    "exploration",          17,                     4,    0.76,    0.46,     0.92,  0.001,     61,   324,           263,
    "struggling_rate",      20,                     2,    0.54,    0.21,     0.88,  0.005,     61,    68,            7, 
    "struggling_rate",      20,                     2,    0.88,    0.78,     0.99,  0.001,  268.5, 271.5,            3,
    "struggling_rate",      20,                     4,    0.15,   -0.31,     0.61,  0.25,      61,   270,           209
) %>% 
    mutate(Key = "Z9FL9L9P",
           species_common = "zebra_finch",
           species_latin = "Taeniopygia_guttata",
           sex = 1,
           context = 1,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = NA,
           R_se = NA,
           remarks = NA,
        max_lifespan_days = zf_lv
         ) -> meta_table
    
write_delim(meta_table, path = "output/David_Auclair_2012.txt", delim = " ", col_names = TRUE)

###### Paper 14: DeWitt_Sih_1999 ####
# DeWitt, TJ; Sih, A; Hucko, JA	1999
# 7QAACTY6	
# trait compensation and cospecialization in a freshwater snail: size, shape and antipredator behaviour

# snail were wild caught

# longevity: 14.5 months
# citation: DeWitt, R. M. (1954). Reproductive capacity in a pulmonate snail (Physa gyrina Say). The American Naturalist, 88(840), 159-164.

max_lifespan_snail <- 14.5*30
avg_adult_age <- max_lifespan_snail * 0.25

# To quantify the repeatability of individual variation in antipredator behaviour, we repeated our
# behavioural assays on 6 days spread over a 13-day period.
tribble(
    ~behaviour,   ~sample_size, ~measurements_per_ind,    ~R, ~CI_lower, ~CI_upper,  ~p_val,    ~t1,   ~t2, ~delta_t,
    "antipredator_behaviour",         96,           3,   0.33,      NA,        NA,   0.0001,     avg_adult_age,    NA,         3,  
    "antipredator_behaviour",         96,           6,   0.27,      NA,        NA,   0.0001,     avg_adult_age,    NA,        13     
) %>% 
    mutate(Key = "7QAACTY6",
           species_common = "physid_snail",
           species_latin = "Physa_gyrina",
           t2 = t1 + delta_t,
           context = 2,
           remarks = "No SE/CI, also data was obtained by pooling across video observation, i.e. measurement per ind not reliable",
           sex = 0,
           context = 2, 
           type_of_treatment = 0,
           treatment = NA, 
           life_stage = NA,
           event = NA,
           R_se = NA,
           max_lifespan_days = max_lifespan_snail) -> meta_table
write_delim(meta_table, path = "output/DeWitt_Sih_1999.txt", delim = " ", col_names = TRUE)

###### Paper 15: Debeffe_Lemaitre_2015 #######
# Debeffe, L.; Lemaitre, J. F.; Bergvall, U. A.; Hewison, A. J. M.; Gaillard, J. M.; Morellet, N.; Goulard, M.; Monestier, C.; David, M.; Verheyden-Tixier, H.; Jaederberg, L.; Vanpe, C.; , Kjell; er, P.	
# 4ABCSFR6	2015	
# short- and long-term repeatability of docility in the roe deer: sex and age matter

# winter 2007e2008 and winter 2012e2013
# N = 130 individuals for a total of 532 capture events, mean ± SE number of captures per individual = 4.09 ± 2.81, range 2e22

# Between-winter adjusted repeatability of standardized handling scores was assessed as a measure of long-term repeatability while controlling for sex differences in repeatability
# and using a subsample of 65 individuals recaptured during successive winters (mean ± SE ¼ 2.68 ± 0.97 winters, range 2e5).
# winter capuring season delta_t = 185 days, between 12 November and 25 March
# between years: 365 + 185

# within season : N = 117
# between seasons: N = 65

# roe deer life expectancy at 1 year of age: 10-12 years
# sample: 66 juveniles (<1 year), 87 adults, i.e. roughly average age: (66 * 1/2 + 87 * 5) / (66 + 87) = 3.1 years = 1132 days
# 2.7: average number of winters an individual is caught: so delta_t (between 2.7 winters) = 1.7 * 365 + 185

#### sex difference is reported but not split into long and short term

deer_data <- scrape_AnAge("Capreolus capreolus", vars = "maximum_longevity_yrs", download_data = FALSE)
deer_max <- as.numeric(deer_data$maximum_longevity_yrs) * 365

delta_t_between <- (1.7*365)+185

tribble(
    ~behaviour,   ~sample_size, ~measurements_per_ind,    ~R, ~CI_lower, ~CI_upper,  ~p_val,    ~t1,                    ~t2,        ~delta_t, ~remarks,
    "docility",         116,           4,               0.41,      0.30,      0.51,      NA,   1132,              1132+185,             185, "delta_t is full season",
    "docility",          65,           2.7,             0.31,      0.13,      0.49,      NA,   1132,  1132+delta_t_between, delta_t_between, "delta_t is full season plus an 2.7*365 days"   
) %>% 
    mutate(Key = "4ABCSFR6",
           species_common = "roe_deer",
           species_latin = "Capreolus_capreolus",
           sex = 0,
           context = 3,
           type_of_treatment = 0,
           life_stage = "both",
           event = NA,
           treatment = NA,
           R_se = NA,
        max_lifespan_days = deer_max) -> meta_table

write_delim(meta_table, path = "output/Debeffe_Lemaitre_2015.txt", delim = " ", col_names = TRUE)








###### Paper 17: Erhard_Mendl_1997 #######

# Erhard, HW; Mendl, M		1997	
# SX3CYX5C
# measuring aggressiveness in growing pigs in a resident-intruder situation

# 1 and 2: 11 weeks of age and 1 day later
pig_ls <- scrape_AnAge("Sus scrofa", vars = "maximum_longevity_yrs", download_data = FALSE)
pig_max_ls <- as.numeric(pig_ls$maximum_longevity_yrs) * 365

tribble(
    ~behaviour,         ~R, ~R_se, ~p_val, ~sample_size, ~measurements_per_ind, ~t1, ~t2, ~delta_t,
    "aggressiveness",  0.56, NA,     0.01,   85,                              2,  77,  78,        1,
    "aggressiveness",  0.73, NA,     0.01,   78,                              2,  77,  78,        1,
    "aggressiveness",  0.57, NA,     0.01,   53,                              2,  49,  77,        28
) %>% 
    mutate(Key = "SX3CYX5C",
           species_common = "pig",
           species_latin = "sus_scrofa",
           sex = 0,
           context = 1,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "juvenile",
           event = NA,
           CI_lower = NA,
           CI_upper = NA,
           remarks = "All spearman rank correlations, no SE, not the same individuals for short and long term rpt",
           max_lifespan_days = pig_max_ls) -> meta_table

write_delim(meta_table, path = "output/Erhard_Mendl_1997.txt", delim = " ", col_names = TRUE)






###### Paper 18: Fisher_James_2015 incomplete, data_to_be_analysed #####

# Fisher, David N.; James, Adele; Rodriguez-Munoz, Rol; , o; Tregenza, Tom	2015	
# behaviour in captivity predicts some aspects of natural behaviour, but not others, in a wild cricket population
# FVQNEVJS
# url with data
# "https://ore.exeter.ac.uk/repository/handle/10871/16930"
url <- "https://ore.exeter.ac.uk/repository/bitstream/handle/10871/16930/bae.w2.txt?sequence=4&isAllowed=y"
wild_shy <- read_delim(url, delim = "\t")
names(wild_shy)
table(wild_shy$tag)
table(wild_shy$date)

###### Paper 19: Gabriel_Black_2010 ######
# Gabriel, Pia O.; Black, Jeffrey M.	2010	
# NWYGHZWI
# behavioural syndromes in steller's jays: the role of time frames in the assessment of behavioural traits


# long term: measurments_per_ind = 3.4, N = 109
# age: 1-11 years, so take 5 years, 1825 days
# sample size long term is guessed
# could be double checked
# I think measures have been averaged across years for long term, so only one datapoint for 2006 and one for 2008

jay_data <- scrape_AnAge(latin_name = "Cyanocitta stelleri", vars = "maximum_longevity_yrs", download_data = FALSE)
max_ls_jay <- as.numeric(jay_data$maximum_longevity_yrs)*365

tribble(
      ~R, ~CI_lower, ~CI_upper, ~sample_size, ~measurements_per_ind,        ~t1,       ~t2,   ~delta_t,  ~p_val, 
    0.49,      0.33,      0.64,          44,                   3.4,        1825,    1825+365,      365,  0.0001,
    0.41,      0.28,      0.54,          65,                   3.4,  1825+2*365,  1825+3*365,      365,  0.0001,
    0.74,      0.49,      0.88,          20,                     2,        1825,  1825+3*365,    3*365,  0.0001   
) %>% 
    mutate(Key = "NWYGHZWI",
           species_common = "stellers_jay",
           species_latin = "Cyanocitta_stelleri",
           behaviour = "risk_taking",
           sex = 0,
           context = 3,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = "between_season",
           R_se = NA,
           remarks = "sample sizes and measurement were a bit of guesswork",
           max_lifespan_days = max_ls_jay ) -> meta_table

write_delim(meta_table, path = "output/Gabriel_Black_2010.txt", delim = " ", col_names = TRUE)



###### Paper 20: Garamszegi_Marko_2015 ######

# Garamszegi, L.Z.; Mark, G.; Szsz, E.; Zsebk, S.; Azcrate, M.; Herczeg, G.; Trk, J. 2015
# among-year variation in the repeatability, within- and between-individual, and phenotypic correlations of behaviors in a natural population
# WFBJ35EC	

# juveniles and adults
# collared flycatcher
# not clear from article how far apart between year repeatabilities are
# basically its within year but different individuals in the 5 years
# and then between years it can be any 2 or more years within these 8 (so minimum 2, maximum 8 years)

# average measurements within year: 2.74

location <- "data/papers/Garamszegi et al. - 2015 - among-year variation in the repeatability, within-.pdf"
# Extract the table
loc_area <- locate_areas(location, pages = 10)    

R_table <- extract_tables(location,
    output = "data.frame",
    pages = c(10), # include pages twice to extract two tables per page
    area = list(
        unlist(loc_area)
    ),
    guess = FALSE
)

# breeding season measurements within one month

# not mentioned how far apart between year measurements were, so 1.5 years between measurements on average assumed
# 8 year study period, most of the 19 between year individuals will have been observed in consecutive years,
# but some further apart, so 1.5*365 days between measurements seems reasonable

flycatch_lf <- scrape_AnAge("Ficedula albicollis", "maximum_longevity_yrs", download_data = FALSE)
max_lifespan <- as.numeric(flycatch_lf$maximum_longevity_yrs) * 365
avg_adult_age <- max_lifespan*0.25

# over course of 8 years, but unclear whats the average distance, so minimum 2 years taken
extract_table <- function(trait) {
    trait <- enquo(trait)
    R_table[[1]] %>% 
        filter(!(Year == "")) %>% 
        select(!!trait, Year) %>% 
        separate(!!trait, into = c("sample_size", "R", "CI"), sep = " ") %>% 
        separate(CI, into = c("CI_lower", "CI_upper"), sep = "/") %>% 
        mutate(CI_lower = str_replace(CI_lower, "\\(", ""), 
            CI_upper = str_replace(CI_upper, "\\)", ""),
            delta_t = c(rep(365, 5), 365 * 2),
            remarks = paste0("year_", Year),
            measurements_per_ind = c(rep(2.7, 5), 2), 
            behaviour = str_replace(!!trait, "\\.", "_"),
            sex = c(rep(0, 5), 2)) %>% 
        select(-Year)
} 

# list to table
meta_table <- bind_rows(lapply(c("Novelty.avoidance", "Aggression", "Risk.taking"), extract_table)) %>% 
                mutate(
                    Key = "WFBJ35EC",
                    species_common = "collared_flycatcher",
                    context = 3,
                    type_of_treatment = 0,
                    treatment = NA,
                    life_stage = "both",
                    event = NA,
                    R_se = NA,
                    p_val = NA,
                    t1 = "average_age",
                    t2 = rep(c(rep("average_age_plus_365", 5), "average_age_plus_at_least_2times365"), 3)
                )

# write_delim(meta_table, path = "output/Garamszegi_Mark_2015.txt", delim = " ", col_names = TRUE)

# redone 
dat <- read_delim("data/Garamszegi_Mark_2015_old.txt", delim = " ")
dat$t2
dat %>% 
    mutate(t1 = avg_adult_age,
           t2 = ifelse(delta_t == 365, avg_adult_age + 30, avg_adult_age + 547.5)) %>% 
    mutate(delta_t = t2 - t1,
           sex = 2,
           species_latin = "Ficedula_albicollis",
           max_lifespan_days = max_lifespan) -> meta_table

write_delim(meta_table, path = "output/Garamszegi_Mark_2015.txt", delim = " ", col_names = TRUE)




###### Paper 21: Gifford_Clay_2014 ###### ######

# Gifford, Matthew E.; Clay, Timothy A.; Careau, Vincent 2014	
# D4KRF54L	
# individual (co) variation in standard metabolic rate, feeding rate, and exploratory behavior in wild-caught semiaquatic salamanders

# adult individuals caught from the wild
# lifetime: up to 10 years?

# three years here
avg_adult <- 5*365

tribble(
    ~behaviour,                 ~R,  ~R_se,  ~sample_size,~measurements_per_ind,  ~delta_t,            ~t1,           ~t2,
    "feeding_rate",          0.778, 0.101,             19,                   2,         21,      avg_adult, avg_adult + 3,
    "feeding_rate",          0.543, 0.192,             19,                   2,         42,  avg_adult + 3, avg_adult + 9,
    "feeding_rate",          0.451, 0.217,             19,                   2,         63,      avg_adult, avg_adult + 9, 
    "exploratory_behaviour", 0.389, 0.237,             19,                   2,         21,      avg_adult, avg_adult + 3,
    "exploratory_behaviour", 0.302, 0.251,             19,                   2,         42,  avg_adult + 3, avg_adult + 9,
    "exploratory_behaviour", 0.147, 0.258,             19,                   2,         63,      avg_adult, avg_adult + 9
) %>% 
    mutate(Key = "D4KRF54L",
           species_common = "ouchita_dusky_salamander",
           species_latin  = "Desmognathus_brimleyorum",
           sex = 0,
           context = 2,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = NA,
           CI_lower = NA,
           CI_upper = NA,
           p_val = NA,
           remarks = "lifetime without academic source",
           max_lifespan_days = 3650) -> meta_table

write_delim(meta_table, path = "output/Gifford_Clay_2014.txt", delim = " ", col_names = TRUE)











###### Paper 23: Grace_Anderson_2014 #####
# Grace, Jacquelyn K.; Anderson, David J. 2014
# X75NVMBP
# personality correlates with contextual plasticity in a free-living, long-lived seabird

# nazga booby
# lifespan up to at least 26 years
# ~43 days egg incubation 15 days nestling brooding, all testing during incubation
# first round of testin:
# nest intruder 444 birds on eggs, 35 on chicks
# novel object: 418 on eggs, 61 on chicks
# second novel object and social stimulus: 409 on eggs, 70 on chicks

# behaviours:
# Gardening / Moving nest material / Mateadvertisement (males), Territy display (both), anxiety or agitation 
# Shaking / Body shakes and shivers / Aggressive signal, anxiety or agitation related
# Aggression / Biting and jabbing intruder, novel object, simulated conspecific / Aggressivenes or Boldness

# four tests :nest intruder, two novel object, one social stimulus
# first novel obj: red bull can, second: plastic crate

# test sequence (all during incubation): 
# session 1: 2008-2009, sample size 51
# session 2: 2009-2010, sample size 71
# session 3: 2010       sample size 15
# session 4: 2011-2012 (repeated birds from session 1-3) sample size 20 / but for long term: 86
# overall sample size short term repeatability: 51+71+15+20 = 157
# overall samlpe size long term repeatability : 86
# intervals: test-retest: average: 14 days (mean sd = 13.8 +- 3.8), N = 193, two measurements per session
#            test-retest long term: min: 1 year, max: 2 years, so lets take mean: 548 days
todigit <- metaDigitise("data/to_digitise/study3_grace2014/")
 # write_delim(todigit, path = "data/to_digitise/study3_grace2014/digitized_figure.txt")
todigit <- read_delim("data/to_digitise/study3_grace2014/digitized_figure.txt", delim = " ")


boobie_ls <- scrape_AnAge(latin_name = "sula granti", vars = "maximum_longevity_yrs", download_data = FALSE)
# to be figured out
booby_max <- as.numeric(boobie_ls$maximum_longevity_yrs)*365
avg_adult_age <- booby_max * 0.25

todigit %>% 
    select(group_id, mean, se, n) %>% 
    rename(R = mean,
           behaviour = group_id,
           R_se = se,
           sample_size = n) %>% 
    mutate(delta_t = c(rep(548, 12), rep(14, 12)),
           measurements_per_ind = c(rep(2, 12), rep(2, 12)), 
           sex = 0,
           context = 3,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = c(rep("between_season", 12), rep("within_season", 12)),
           CI_lower = NA,
           CI_upper = NA,
           p_val = NA,
           t1 = avg_adult_age,
           t2 = t1 + delta_t,
           remarks = "age_unknown",
           Key = "X75NVMBP",
           species_common = "nazca_booby",
           species_latin = "sula_granti",
           max_lifespan_days = booby_max) -> meta_table

write_delim(meta_table, path = "output/Grace_Anderson_2014.txt", delim = " ", col_names = TRUE)



###### Paper 24: Greggor_Jolles_2016 ######

# Greggor, Alison L.; Jolles, Jolle W.; Thornton, Alex; Clayton, Nicola S.	2016
# 6GYLQSZ7		
# seasonal changes in neophobia and its consistency in rooks: the effect of novelty type and dominance position

# rooks
# risk-taking 
# captive rooks
# long term: 4 year period / neophobia and dominance
# 19 initially, 16 in the end
# birds were 7 years old during first testing 2010 (and then 2014 )
# breeding season: end of feburary / early march
# 3 measurements in 2014: breeding season 2014, summer and winter // dominance hierarchies

7*365
11*365

# breeding to non-breeding season ~ 365/2
# assumed length of breeding season ~1 month
# assumed length of non-breeding season measurements: ~1 month

rook_data <- scrape_AnAge(latin_name = "Corvus frugilegus", vars = "maximum_longevity_yrs", download_data = FALSE)
rook_max <- as.numeric(rook_data$maximum_longevity_yrs) * 365

tribble(
    ~behaviour,          ~R, ~CI_lower, ~CI_upper,  ~t1,  ~t2, ~delta_t,                                              ~event,   ~remarks,
    "dominance_ranks", 0.77,      0.45,      0.91, 2555, 4015,     1460,             "across_breeding_seasons_in_both_years",   "spearman_corr",
    "dominance_ranks", 0.42,     -0.10,      0.76, 4015, 4198,      183, "within_year_change_breeding_to_nonbreeding_season",   "spearman_corr",
    "novel_object",    0.28,      0.16,      0.49, 4015, 4045,       30,                            "within_breeding_season",   "neophobia",
    "novel_object",    0.15,      0.06,      0.34, 4198, 4228,       30,                         "within_nonbreeding_season",   "neophobia",
    "novel_object",    0.49,     -0.01,      0.79, 4015, 4198,      183,             "change_breeding_to_nonbreeding_season",   "neophobia",
    "novel_object",    0.48,      0.00,      0.78, 2555, 4015,     1460,             "across_breeding_seasons_in_both_years",   "spearman_corr",
    "novel_people",    0.16,      0.06,      0.36, 4015, 4045,       30,                            "within_breeding_season",   "neophobia",
    "novel_people",    0.00,      0.00,      0.07, 4198, 4228,       30,                         "within_nonbreeding_season",   "neophobia",
    "novel_people",    0.55,      0.07,      0.82, 4015, 4198,      183, "within_year_change_breeding_to_nonbreeding_season",   "spearman_corr"
) %>% 
    mutate(Key = "6GYLQSZ7",
           species_common = "rook",
           species_latin = "Corvus_frugilegus",
           sample_size = 16,
           measurements_per_ind = c(2, 2, 2, 2, 4, 2, 2, 2, 4),
           sex = 0,
           context = 1,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           R_se = NA,
           p_val = NA,
           max_lifespan_days = rook_max) -> meta_table

write_delim(meta_table, path = "output/Greggor_Jolles_2016.txt", delim = " ", col_names = TRUE)



###### Paper 25: Guenther_Groothuis_2018,  data_to_be_analysed ############
# Guenther, A.; Groothuis, A. G. G.; Krueger, O.; Goerlich-Jansson, V. C.	2018	
# 	4FC9GYVQ

# data to be analysed and models are described in paper

# cortisol during adolescence organises personality traits and behavioural syndromes
dat <- read_xlsx("data/downloaded_from_papers/raw_data_Guenther_et_al_Hormones&Behaviour.xlsx") %>% 
    rename(testround = age)
table(dat$ID)

# measures:
# novel_object / boldness: No_touches
# exploration: LF_trips
# struggle docility: Struggle

dat$testround
# models: (testround is testround here, going from 1-3)
# fixed: Treatment + sex + testround + treatment*sex + testround*sex + (1|mother/ID)
# for sociapositive and aggressivenes: + (1|Stimulus)

# gaussian: hand-escape latency, struggle docility
# poisson: boldness, exploration, sociopositive
# binary: aggressiveness

# tests: first test: 61 days old, (early adolescence)
# second: 92 days (late adolescence)
# third: 123 days (adult)

# measures detail:
# Struggle_docility: Gaussian
dat$Struggle
# Hand_escape_latency: Gaussian
dat$Hand
# Boldness, Novel object: Poisson
dat$NO_touches
# Free exploration (number of trips to new area): Poisson
dat$LF_trips
# sociopositive behaviour: Poisson
dat$Socpos
# aggressive: Binary
dat$Aggression_binary

dat_final <- dat %>% 
    select(ID, sex, Treat, mother, Stimulus, testround, Struggle, Hand, NO_touches, 
           LF_trips, Socpos, Aggression_binary) %>% 
    mutate_at(vars(6:11), funs(as.numeric))


# not_roundx is which measurement round from those three not to take
# creates meta_table
calc_rpt <- function(resp, treat = "T", not_roundx, dat){
    
    nboot <- 100
   # dat <- dplyr::filter(dat, ((Treat == treat) & (testround != not_roundx)))
    dat <- dplyr::filter(dat, testround != not_roundx)
    
    if (resp %in% c("Struggle", "Hand")) {
        mod_form <- as.formula(paste0(resp, " ~ sex + testround + testround * sex + (1|ID) + (1|mother)"))
        res <- rpt(mod_form, data = dat, grname = "ID",
            datatype = "Gaussian", nboot = nboot)
    } 
    if (resp %in% c("NO_touches", "LF_trips")) {
        mod_form <- as.formula(paste0(resp, " ~ sex + testround + testround * sex + (1|ID) + (1|mother)"))
        res <- rpt(mod_form, data = dat, grname = "ID",
            datatype = "Poisson", nboot = nboot)
    }
    
    if (resp == "Socpos") {
        mod_form <- as.formula(paste0(resp, " ~ sex + testround + (1|ID) + (1|mother)")) #(1|Stimulus) + testround * sex 
        res <- rpt(mod_form, data = dat, grname = "ID",
            datatype = "Poisson", nboot = nboot)
    }
    
    if (resp == "Aggression_binary") {
        mod_form <- as.formula(paste0(resp, " ~ sex + testround + (1|ID) + (1|mother)")) # (1|Stimulus) + testround * sex
        res <- rpt(mod_form, data = dat, grname = "ID",
            datatype = "Binary", nboot = nboot)
    }

    treat_name <- ifelse(treat == "T", "cortisol", "control")
    
    ages_at_roundx <- data.frame("round_1" = 61, "round_2" = 92, "round_3" = 123)
    ages <- ages_at_roundx %>% 
               # unselect the round not used 
               select(-grep(pattern = not_roundx, x = names(ages_at_roundx))) 
    
    if (resp %in% c("Struggle", "Hand")) {
        R <- unname(res$R)
        R_se <- unname(res$se)
    } else {
        R <- unname(res$R[2, ])
        R_se <- unname(res$se[2, ])
    }
    out <- data.frame(R = R, 
                      R_se = R_se, 
                      sample_size = res$nobs,  
                      behaviour = resp, 
                      treatment = treat_name,
                      t1 = unname(ages[[1]]),
                      t2 = unname(ages[[2]])) %>% 
                 mutate(delta_t = t2 - t1)
    out <- as_tibble(out)
}


# set up variables for all models

# takes a long time so split up:
# first four responses
resps <- c("Struggle", "Hand", "NO_touches", "LF_trips" ) #"Socpos", "Aggression_binary"
treats <- c("T", "C")
not_rounds <- c(3,1,2) # first measurements close in time and at last the two further away

df_vars <- data.frame(all_resps = rep(resps, 6), 
                 all_treats = rep(rep(treats, each = length(resps)), 3), 
                 all_not_rounds = rep(not_rounds, each = length(resps)*2), stringsAsFactors = FALSE)


all_rpts <- apply(df_vars, 1, function(x) calc_rpt(unname(x[[1]]), unname(x[[2]]), unname(x[[3]]), dat_final))

# reshape some stuff which wasn't right
all_rpts_df <- rbind_list(all_rpts) %>% 
                mutate(R_se = ifelse(is.na(se), R_se, se)) %>% 
                select(-se)
# write_delim(all_rpts_df, path = "data/anja_paper_rpts1.txt")

# next response:
resps <- c("Socpos") #"Socpos", "Aggression_binary"
# treats <- c("T", "C")
not_rounds <- c(3,1,2) # first measurements close in time and at last the two further away

df_vars <- data.frame(all_resps = rep(resps, 3), 
    # all_treats = rep(rep(treats, each = length(resps)), 3), 
    all_not_rounds = rep(not_rounds, each = length(resps)), stringsAsFactors = FALSE)

all_rpts2 <- apply(df_vars, 1, function(x) calc_rpt(unname(x[[1]]), not_roundx = unname(x[[2]]), dat = dat_final))
all_rpts3 <- apply(df_vars, 1, function(x) calc_rpt(unname(x[[1]]), not_roundx = unname(x[[2]]), dat = dat_final))



#  "Socpos", "Aggression_binary"

rpt(Struggle ~ sex + testround + testround * sex + (1 | ID) + (1 | mother), data = dat, grname = "ID",
    datatype = "Gaussian", nboot = 10)



###### Paper 26: Haage_Bergvall_2013 ######
# Haage, Marianne; Bergvall, Ulrika A.; Maran, Tiit; Kiik, Kairi; Angerbjorn, Anders	2013	
# situation and context impacts the expression of personality: the influence of breeding season and test context	
# GZFPB2J7
# european_mink
# captive bred

# within season: November/December 2009, n = 80 (males and females) 
# between season also: march april 2010 (breeding season): N = 68
# from juvenile to adult (0 to 6 years old) / average 3 years old?

avg_adult <- 1095

# mink_data <- scrape_AnAge("european mink", vars = "maximum_longevity_yrs", download_data = FALSE)
max_age_mink <- 3650 # citation https://en.wikipedia.org/wiki/Mink

tribble(
    ~R,    ~behaviour,     ~sample_size, ~measurements_per_ind, ~delta_t,   ~t1, ~t2,  ~event,
    0.59, "novel_object",            80,                     2,       30,  1095, 1125, "non_breeding_season",
    0.69, "novel_object",            68,                     4,      120,  1095, 1215, "non_breeding_to_breeding_season_rpt" 
) %>%
    mutate(
        Key = "GZFPB2J7",
        species_common = "european_mink",
        species_latin = "Mustela_lutreola",
        sex = 0,
        context  = 1,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "both",
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        remarks = "no se or ci or pval, maybe simulation?",
        max_lifespan_days = max_age_mink) -> meta_table
    
write_delim(meta_table, path = "output/Haage_Bergvall_2013.txt", delim = " ", col_names = TRUE)


###### Paper 27: Hammond-Tooke_Cally_2012, metadigitise here #######   
### metadigitise
# Hammond-Tooke, Cally A.; Nakagawa, Shinichi; Poulin, Robert 2012	
# 6TGIGAET	
# parasitism and behavioural syndromes in the fish gobiomorphus cotidianus

# max age giant bully: 10 years #Jellyman, D. J., Sagar, P. M., Glova, G. J., & Sykes, J. R. E. (2000). Age, growth, and movements of giant bullies (Gobiomorphus gobioides) in the Kakanui River estuary, South Island, New Zealand. New Zealand Journal of Marine and Freshwater Research, 34(3), 523-530.
max_age_bully <- 3650
# common_bullies
# Gobiomorphus cotidianus
# "Young bullies" between 25 and 55 mm in length

# fish from wild to lab
 # three sessions, three weeks between each
# one session : fish tested twice, 7 days apart
# second session with predator cue added

# all young fish, so guess: half a year old (they mature at 1 year)
avg_juvenile <- 183

scrape_AnAge(latin_name = "Gobiomorphus cotidianus", vars = "maximum_longevity_yrs", download_data = FALSE)

# digitalise
# todigit <- metaDigitise("data/to_digitise/study4_hammond-tooke_cally2012/")
# write_delim(x = todigit, path = "data/to_digitise/study4_hammond-tooke_cally2012/rpt_table.txt")
todigit <- read_delim("data/to_digitise/study4_hammond-tooke_cally2012/rpt_table.txt", delim = " ")

todigit %>% 
    select(group_id, mean, n, se) %>% 
    rename(R = mean,
           R_se = se,
           sample_size = n) %>% 
    separate(group_id, into = c("session", "behaviour"), sep = "_") %>% 
    mutate(t1 = c(rep(avg_juvenile, 6), rep(avg_juvenile + 3*7, 3), rep(avg_juvenile + 7*7, 3)),
           t2 = c(rep(avg_juvenile + 8*7, 3), rep(avg_juvenile + 7, 3), rep(avg_juvenile + 4*7, 3), rep(avg_juvenile + 8*7, 3)),
           delta_t = t2-t1,
           Key = "6TGIGAET",
           species_common = "common_bully",
           species_latin = "Gobiomorphus_cotidianus",
           measurements_per_ind = c(rep(6, 3), rep(2, 9)),
           sex = 0,
           context = 2,
           type_of_treatment = 2,
           treatment = c(rep("predator_cue", 3), rep("no_cue", 3), rep("predator_cue", 3), rep("no_cue", 3)),
           life_stage = "juvenile",
           event = NA, 
           p_val = NA,
           CI_lower = NA,
           CI_upper = NA,
           remarks = "age_guessed, double check",
           max_lifespan_days = max_age_bully) -> meta_table
    
write_delim(meta_table, path = "output/Hammond-Tooke_Cally_2012.txt", delim = " ", col_names = TRUE)



 



###### from now on: avg_adult_age = 0.25 * max_longevity (from http://genomics.senescence.info/species/) ######
#####              avg_juvenile_age = 0.5 * age_maturity
##### variable    max_lifespan_days added to meta_table



###### Paper 28: Herde_Eccard_2013, tabulizer here ##############     
## tabulized extractio here
# Herde, Antje; Eccard, Jana A.	2013
# consistency in boldness, activity and exploration at different stages of life
# CVPFUYZS	

?tabulizer

location <- "data/papers/Herde and Eccard - 2013 - consistency in boldness, activity and exploration .pdf"
# Extract the table
# out <- extract_tables(location, pages = 3)
# out
# locate_areas(location, pages = 3)    

R_table <- extract_tables(location,
    output = "data.frame",
    pages = c(3), # include pages twice to extract two tables per page
    area = list(c(101.44020 , 55.68943, 434.54786, 546.44064)),
    guess = FALSE
) 

# get number of measurements
only_num_meas <- R_table[[1]] %>% 
    as_tibble() %>% 
    rename(mean_sd = 4) %>% 
    filter(grepl("\\(", mean_sd)) %>% 
    select(mean_sd) %>% 
    mutate(mean_sd = str_replace(mean_sd, "\\(", "")) %>% 
    mutate(mean_sd = str_replace(mean_sd, "\\)", "")) %>% 
    mutate(mean_sd = as.numeric(mean_sd)) %>% 
    rename(num_measurements = mean_sd) 

# transform 
R_table[[1]] %>% 
    as_tibble() %>% 
    rename(mean_sd = 4,
           over_mat = 5,
           adult_short = 6,
           adult_long = 7) %>% 
    filter(!grepl("\\(", mean_sd)) %>% 
    slice(2:n()) %>% 
    mutate(num_measurements = only_num_meas$num_measurements) %>% 
    mutate(Test = na_if(Test, ""),
           over_mat = na_if(over_mat, "")) %>% 
    fill(Test) %>% 
    # behaviour names according to cluster analysis
    mutate(behaviour = c("boldness1", "activity1", "activity2", "boldness2", "activity3", 
         "boldness3",  "boldness4", "activity4", "exploration1", "exploration2", "activity5")) %>% 
    # mutate(Definition = str_replace_all(Definition, " ", "_"),
    #        Variable = str_replace_all(Variable, " ", "_")) %>% 
    # unite(Test, Variable, Definition, col = "behaviour", sep = "_") %>% 
    gather(over_mat, adult_short, adult_long, key = "timespan", value = "rpt") %>% 
    select(behaviour, timespan, num_measurements, rpt) %>% 
    filter(!is.na(rpt)) %>% 
    separate(rpt, into = c("sample_size", "R", "p_val"), sep = " ") -> meta_table_raw

# the minus in -0.031 isn't recognized properly, replace here with actual minus string
meta_table_raw[meta_table_raw$behaviour == "activity5" & meta_table_raw$timespan == "adult_long", "R"] <- "-0.031"

# common_voles 
# Microtus_arvalis
# 1) before and after maturation over three month 62 +- 20 days old and 90 days later after mat. lab-born , 17 voles 9 males 8 females
# 2) adult during one week / different voles, adult voles (avg life-span 4.5 month), 88 males, 80 females, wild trapped and lab tested
# 3) adult over three month, 48 adults males and females 2.5 month apart wild trapped and lab tested
# 
# assumed average adult age: 0.25 * 4.8 year = 1.2 years = 438 days

avg_adult_age <- 438

meta_table <- meta_table_raw %>% 
    mutate(Key = "CVPFUYZS",
           sample_size = as.numeric(sample_size),
           R = as.numeric(R),
           measurements_per_ind = num_measurements/sample_size,
           t1 = c(62, 62, rep(avg_adult_age , 22)),
           t2 = case_when(
               timespan == "over_mat" ~ 152,
               timespan == "adult_short" ~ avg_adult_age  + 7,
               timespan == "adult_long" ~ avg_adult_age + 76
           ),
           delta_t = t2 - t1,
           event = c(rep("maturation", 2), rep(NA, 22)),
           context = c(1, 1, rep(2, 22)),
           species_common = "common_vole",
           species_latin = "Microtus_arvalis",
           sex = 0,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = c(rep("both", 2), rep("adult", 22)),
           R_se = NA,
           CI_lower = NA,
           CI_upper = NA,
           remarks = "no uncertainty yet",
           max_lifespan_days = 438) %>% 
        select(-c(num_measurements, timespan))

meta_table

write_delim(meta_table, path = "output/Herde_Eccard_2013.txt", delim = " ", col_names = TRUE)


###### Paper 29: Heuschele_Ekvall_2017 #####

# Heuschele, Jan; Ekvall, Mikael T.; Bianco, Giuseppe; , Hyl; er, Samuel; Hansson, Lars-Anders 2017	
# 	"N7X6VHZ5"
# context-dependent individual behavioral consistency in daphnia
# Ultraviolet radiation 
# longevity: 0.19 years = 69.35 days   69.35 * 0.25 = 17 days average adult age
avg_adult_age <- 17
avg_juvenile_age <- 1

# Daphnia magna
# collected in the wild
# adults measured over 3 days
# juveniles measured over 3 weeks
library(tibble)

# Pietrzak, B., Bednarska, A., & Grzesiuk, M. (2010). Longevity of Daphnia magna males and females. Hydrobiologia, 643(1), 71-75.
max_age_daphnia <- 100 

tribble(
             ~behaviour,    ~R, ~CI_lower, ~CI_upper, ~delta_t,           ~t1,              ~t2,     ~type_of_treatment,            ~treatment,
    "swimming_velocity",  0.73,      0.45,     0.92,        3, avg_adult_age,        avg_adult_age + 3,     2,            "before_UV_radiation",
    "swimming_velocity",  0.177,    -0.176,   0.614,        3, avg_adult_age,        avg_adult_age + 3,     2,                   "UV_radiation",
    "swimming_velocity",  0.576,     0.207,   0.846,        3, avg_adult_age,        avg_adult_age + 3,     2,             "after_UV_radiation",
    "swimming_velocity",  0.357,      0.097,  0.622,       21, avg_juvenile_age, avg_juvenile_age + 21,     2,            "before_UV_radiation",
    "swimming_velocity",  0.090,     -0.135,  0.387,       21, avg_juvenile_age, avg_juvenile_age + 21,     2,                   "UV_radiation",
    "swimming_velocity",  0.107,     -0.121,  0.404,       21, avg_juvenile_age, avg_juvenile_age + 21,     2,              "after_UV_radiation"
    ) %>% 
    mutate(
     Key = "N7X6VHZ5",
     species_common = "Daphnia_magna",
     species_latin = "Daphnia_magna",
     measurements_per_ind = 3, 
     sample_size = rep(c(11,22), each = 3),
     sex = 1,
     context = rep(c(2, 1), each = 3),
     life_stage = rep(c("adult", "juvenile"), each = 3),
     event = rep(c(NA, "maturation"), each = 3),
     R_se = NA,
     p_val = NA,
     max_lifespan_days = max_age_daphnia,
     remarks = "juveniles tested long term are offspring of mothers tested short term, i.e. same genes because clonal"
    ) -> meta_table

write_delim(meta_table, path = "output/Heuschele_Ekvall_2017.txt", delim = " ", col_names = TRUE)


###### Paper 30: Holtmann_Santos_2017, data is there, needs to be analysed ####
# Holtmann, Benedikt; Santos, Eduardo S. A.; Lara, Carlos E.; Nakagawa, Shinichi 2017	
# personality-matching habitat choice, rather than behavioural plasticity, is a likely driver of a phenotype-environment covariance
# V3JKY9ME

url <- "https://datadryad.org/bitstream/handle/10255/dryad.154721/Flight-initiation-distance_data.xlsx?sequence=1"
download.file(url, destfile = "data/downloaded_from_papers/holtmann_santos_2017.xlsx")

dat <- read_xlsx("data/downloaded_from_papers/holtmann_santos_2017.xlsx")



###### Paper 31: Hudson_Rangassamy_2015 #######
# Hudson, Robyn; Rangassamy, Marylin; Saldana, Amor; Banszegi, Oxana; Roedel, Heiko G.	2015	
# stable individual differences in separation calls during early development in cats and mice
# P5DGEWV5

# dat <- metaDigitise("data/to_digitise/study5_Hudson_Rangassamy_2015/", summary = FALSE)
# write_delim(dat, "data/to_digitise/study5_Hudson_Rangassamy_2015/digitised_dat.txt")
# dat <- read_delim(file = "data/to_digitise/study5_Hudson_Rangassamy_2015/digitised_dat.txt", delim = " ")
# dat

library(digitize)

mydata <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot1.png")
mydata2 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot2.png")
mydata3 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot3.png")
mydata4 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot4.png")
mydata5 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot5.png")
mydata6 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot6.png")
mydata7 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot7.png")
mydata8 <- digitize("data/to_digitise/study5_Hudson_Rangassamy_2015/plot8.png")


# saveRDS(ls(), file = "data/to_digitise/study5_Hudson_Rangassamy_2015/all_digitised_df.RData")
# obj <- readRDS(file = "data/to_digitise/study5_Hudson_Rangassamy_2015/all_digitised_df.RData")
# make big dataframe
df_list <- mget(paste0("mydata", c("",2:8)))
var_names <- rep(list(c("week1", "week2"), c("week2", "week3"), c("week3", "week4"), c("week1", "week4")), 2)

name_df <- function(df, df_names) {
    names(df) <- df_names
    df
}

dfs_proc <- map2(df_list, var_names, name_df) %>% 
             map(function(x) mutate(x, id = paste0("id", 1:nrow(x))))
            

# reshape data a litte
long_format_dfs <- map(dfs_proc, gather, key = "timepoint", value = "sep_calls", -id) 
long_format_dfs[5:8] <- map(long_format_dfs[5:8], function(x) rename(x, locomotion = "sep_calls"))

# locomotion sometimes negative (due to taking data from plot), make it 0 
long_format_dfs[5:8] <- map(long_format_dfs[5:8], function(x) mutate(x, locomotion = ifelse(locomotion < 0, 0, locomotion)))

# round calls
long_format_dfs[1:4] <- map(long_format_dfs[1:4], function(x) mutate(x, sep_calls = round(sep_calls, digits = 0)))
# repatabilities: calls: poisson, locomotion: gaussian
get_rpts_pois <- function(df) {
    res <- rptPoisson(sep_calls ~ (1|id), grname = "id", nboot = 1000, data = df)
    out <- tibble(R = res$R[2, ], R_se = res$se[2, ], 
                  t1 = unique(df$timepoint)[1], t2 = unique(df$timepoint)[2],
                  behaviour = "separation_calls")
}

get_rpts_gauss <- function(df) {
    res <- rptGaussian(locomotion ~ (1|id), grname = "id", nboot = 1000, data = df)
    out <- tibble(R = as.numeric(res$R), R_se = as.numeric(res$se), 
        t1 = unique(df$timepoint)[1], t2 = unique(df$timepoint)[2],
        behaviour = "separation_locomotor_activity")
}

all_call_rpts <- map_dfr(long_format_dfs[1:4], get_rpts_pois)
all_motor_rpts <- map_dfr(long_format_dfs[5:8], get_rpts_gauss)

meta_table_cat <- rbind(all_call_rpts, all_motor_rpts) %>% 
    mutate(sample_size = unlist(map(long_format_dfs, function(x) nrow(x) / 2)),
           Key = "P5DGEWV5",
           species_common = "domestic_cat",
           species_latin = "felis_catus",
           measurements_per_ind = 2,
           sex = 0,
           context = 1,
           type_of_treatment = 2,
           treatment = "separation_from_mother",
           life_stage = "juvenile",
           event = NA,
           CI_lower =NA,
           CI_upper = NA,
           p_val = NA,
           remarks = "not controlled for litter here, as variables were extracted from plots",
           max_lifespan_days = 10950) %>% 
    mutate(t1 = case_when(
        t1 == "week1" ~ 7,
        t1 == "week2" ~ 14,
        t1 == "week3" ~ 21,
        t1 == "week4" ~ 28
    ), 
        t2 = case_when(
            t2 == "week1" ~ 7,
            t2 == "week2" ~ 14,
            t2 == "week3" ~ 21,
            t2 == "week4" ~ 28
        ),
        delta_t = t2-t1)
    
  

meta_table_mouse <- tribble(
                     ~behaviour,  ~t1, ~t2, ~delta_t,    ~R,  ~CI_lower, ~CI_upper, ~sample_size, ~measurements_per_ind,
             "separation_calls", 13,  14,         1, 0.635,    0.220,     0.877,       18,                        2,
             "separation_calls", 14,  16,         2, 0.504,    0.332,     0.745,       59,                        2, 
"separation_locomotor_activity", 13,  14,         1, 0,        0,         0.415,       18,                        2,
"separation_locomotor_activity", 14,  16,         2, 0.248,    0.001,     0.488,       59,                        2
) %>% 
    mutate(
        Key = "P5DGEWV5",
        species_common = "mound-building_mouse",
        species_latin = "Mus_spicilegus",
        sex = 0,
        context = 1,
        type_of_treatment = 2,
        treatment = "separation_from_mother",
        life_stage = "juvenile",
        event = NA,
        R_se = NA,
        max_lifespan_days = 1278,
        remarks = "different mice for both rpts",
        pval = NA
    )
    
  
    
meta_table <- bind_rows(meta_table_cat, meta_table_mouse)

write_delim(meta_table, path = "output/Hudson_Rangassamy_2015.txt", delim = " ", col_names = TRUE)


###### Paper 32: Jablonszky_Stasz_2017 ######

# Jablonszky, Monika; Szasz, Eszter; Marko, Gabor; Torok, Janos; Herczeg, Gabor; Garamszegi, Laszlo Zsolt	2017
# escape ability and risk-taking behaviour in a hungarian population of the collared flycatcher (ficedula albicollis)
# 	ZCSXGUYW

# collared flycatcher
# 2009 - 2016
# april 11 and may 10 courtship



# may 17 and june 22 chick feeding
dmy(22062016) - dmy(11042016)

# only males
# 246 adults 134 yearlings, 80 birds more than once
# 2016: twice a day
# Escape ability

# digitize data to estimate R_se
library(digitize)
mydata <- digitize("data/to_digitise/study6/plot1.png")
mydata2 <- digitize("data/to_digitise/study6/plot2.png")
mydata3 <- digitize("data/to_digitise/study6/plot3.png")


df_list <- mget(paste0("mydata", c("",2:3)))
var_names <- list(c("day11", "day12"), c("season1", "season2"), c("year1", "year2"))

name_df <- function(df, df_names) {
    names(df) <- df_names
    df
}

dfs_proc <- map2(df_list, var_names, name_df) %>% 
    map(function(x) mutate(x, id = paste0("id", 1:nrow(x))))


# reshape data a litte
long_format_dfs <- map(dfs_proc, gather, key = "timepoint", value = "escape_ability", -id) 

# 

rpt1 <- rptGaussian(escape_ability ~ (1|id), data = long_format_dfs[[1]], grname = "id", nboot = 1000)
rpt2 <- rptGaussian(escape_ability ~ (1|id), data = long_format_dfs[[2]], grname = "id", nboot = 1000)
rpt3 <- rptGaussian(escape_ability ~ (1|id), data = long_format_dfs[[3]], grname = "id", nboot = 1000)


# R_se day: 0.122
# R_se between season: ~0.128
# R_se between years: 0.134
# within-day repeatability (2016): N = 45, R = 0.236, p = 0.022
# between periods: N = 52, R = 0.041, P = 0.356
# between years: N = 34, R = 0, p = 0.758
#
anage <- AnAgeScrapeR::scrape_AnAge("Ficedula albicollis", vars = "maximum_longevity_yrs")
avg_adult_age <- as.numeric(anage$maximum_longevity_yrs) * 0.25 * 365

tribble(
    ~R,   ~R_se,  ~p_val,          ~t1,                  ~t2,    ~delta_t, ~sample_size,
    0.236, 0.122, 0.022, avg_adult_age,  avg_adult_age + 0.5,         0.5,           45,
    0.041, 0.128, 0.356, avg_adult_age,  avg_adult_age + 72,           72,            52,
    0    , 0.134, 0.758, avg_adult_age,  avg_adult_age + 365,         365,            34
) %>% 
    mutate(Key = "ZCSXGUYW",
        species_common = "collared_flycatcher",
        species_latin = "Ficedula_albicollis",
        measurements_per_ind = 2,
        sex = 2,
        behaviour = "escape_ability_risk_taking",
        context = 3,
        type_of_treatment = 0,
        treatment = 0,
        life_stage = "adult",
        event = c(NA, "courtship_to_chick_feeding", "across_years"),
        CI_lower = NA,
        CI_upper = NA,
        max_lifespan_days = as.numeric(anage$maximum_longevity_yrs) * 365,
        remarks = "R_se estimated from recalculating rpt with a simplified model (only id as random effect)") -> meta_table

write_delim(meta_table, path = "output/Jablonszky_Stasz_2017.txt", delim = " ", col_names = TRUE)






 
###### Paper 33 Jennings_Hayden_2013 ######
# Jennings, Domhnall J.; Hayden, Thomas J.; Gammell, Martin P.	2013
# personality and predictability in fallow deer fighting behaviour: the relationship with mating success
# 	NPSM25PQ

# within rut (over two years though) ~ two weeks 14 days
# between ruts 365 days
# adult males
# fighting in the wild -> escalation rate

# 1996: 5494 interactions, 73 mature males = 75 measurements_per_ind
# 1997: 7202 interactions, 74 mature males 
# over all: 12696 int., 147 mature males, = 86 measurements per ind
# behaviour measured as number of fights / number of all interactions
# 7 blocks a 2 days, to 7 measurements per ind I guess
#### 56 over two ruts were included, minus 10 which were observed in both ###

adult_lifespan <- AnAgeScrapeR::scrape_AnAge("Dama dama", vars = c("maximum_longevity_yrs"), download_data = FALSE)
avg_adult_age <- as.numeric(adult_lifespan $maximum_longevity_yrs) * 365 * 0.25

tribble(
       ~R,    ~p_val, ~sample_size, ~measurements_per_ind,           ~t1,                ~t2, ~delta_t,
    0.113,    0.001,            46,                     7, avg_adult_age, avg_adult_age + 14,       14,
    0.221,    0.001,            10,                    14, avg_adult_age, avg_adult_age + 365,     365
) %>% 
    mutate(
        Key = "NPSM25PQ",
        species_common = "fallow_deer",
        species_latin = "Dama_dama",
        sex = 2,
        behaviour = "fighting_behaviour_during_rut",
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = c("within_rut", "between_years"),
        R_se = NA,
        CI_upper = NA,
        CI_lower = NA,
        remarks = "No se reported",
        max_lifespan_days = as.numeric(adult_lifespan$maximum_longevity_yrs) * 365
    ) -> meta_table

write_delim(meta_table, path = "output/Jennings_Hayden_2013.txt", delim = " ", col_names = TRUE)


###### Paper 34: Kluen_Brommer_2013 ####
# Kluen, Edward; Brommer, Jon E.	2013	
# context-specific repeatability of personality traits in a wild bird: a reaction-norm perspective
# 	BJKBMZRI

# neophobia, activity, escape time
# 73 individuals measured in both seasons (presumable within one year)
# blue tit breeding season: may/june
# trapping of adult birds when nestlings were between 6 and 16 days old so within 10 days
# in winter trapping around 3 month max, so probably delta_t around one and a half month ~ 45 days

# winter, say roughly 183 days between seasons

adult_lifespan <- AnAgeScrapeR::scrape_AnAge("Cyanistes caeruleus", vars = c("maximum_longevity_yrs"), download_data = FALSE)
avg_adult_age <- as.numeric(adult_lifespan $maximum_longevity_yrs) * 365 * 0.25


tribble(
    ~behaviour,         ~R,  ~p_val, ~sample_size, ~measurements_per_ind,           ~t1,                 ~t2, ~delta_t, ~event,                                 ~remarks,
    "escape_behaviour", 0.12, 0.11,           175,                     2, avg_adult_age,  avg_adult_age + 10,       10,  "within_breeding_season",                     NA,
    "escape_behaviour", 0.32, 0.001,          175,                     2, avg_adult_age,  avg_adult_age + 45,       45,  "within_winter_season",                       NA,
    "escape_behaviour",-0.013, 0.92,          73,                      2, avg_adult_age, avg_adult_age + 183,      183,  "between_breeding_and_winter_seasons", "correlation",
    "activity",         0.24, 0.002,          175,                     2, avg_adult_age,  avg_adult_age + 10,       10,  "within_breeding_season",                     NA,
    "activity",         0.18, 0.04,           175,                     2, avg_adult_age,  avg_adult_age + 45,       45,  "within_winter_season",                       NA,
    "activity",         0.424, 0.001,          73,                     2, avg_adult_age, avg_adult_age + 183,      183,  "between_breeding_and_winter_seasons", "correlation",
    "neophobia",        0,     0.5,           175,                     2, avg_adult_age,  avg_adult_age + 10,       10,  "within_breeding_season",                     NA,
    "neophobia",        0.46,  0.001,         175,                     2, avg_adult_age,  avg_adult_age + 45,       45,  "within_winter_season",                       NA,
    "neophobia",        0.021, 0.86,            73,                    2, avg_adult_age, avg_adult_age + 183,      183,  "between_breeding_and_winter_seasons", "correlation"
) %>% 
    mutate(
        Key = "BJKBMZRI",
        species_common = "blue_tit",
        species_latin  = "Cyanistes_caeruleus",
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        R_se = NA,
        CI_lower =NA,
        CI_upper =NA,
        max_lifespan_days =  as.numeric(adult_lifespan $maximum_longevity_yrs) * 365
    ) -> meta_table 

write_delim(meta_table, path = "output/Kluen_Brommer_2013.txt", delim = " ", col_names = TRUE)


###### Paper 35: Koski_2011 ######
# Koski, Sonja E.	2011	
# social personality traits in chimpanzees: temporal stability and structure of behaviourally assessed personality traits in three captive populations
# 	INKLPDFJ

library(tabulizer)
location <- "data/papers/Koski - 2011 - social personality traits in chimpanzees temporal.pdf"
location <- "data/papers/split_koski/Koski rot-rotated.pdf"
# split_pdf(location, outdir = "data/papers/split_koski")
out <- locate_areas(location)    

R_table1 <- extract_tables(location,
    output = "matrix",
    #pages = c(7), # include pages twice to extract two tables per page
    area = out,
    guess = FALSE
) 
meta_table0 <- R_table1[[1]] %>% 
    as_tibble() %>% 
    .[-c(1,2, 18:25), ] %>% 
    mutate(var_name = substr(V1, 1, str_locate(V1, "[0-9]")-2),
           numbers = substr(V1, str_locate(V1, "[0-9]"), str_length(V1))) %>% 
    select(-V1) %>% 
    rename(numbers2 = V2) %>% 
    select(var_name, numbers, numbers2) %>% 
    separate(numbers, into = c("R1", "CI_lower1", "CI_upper1", "F1", "p_val1", 
                               "R2", "CI_lower2", "CI_upper2", "F2", "p_val2", "R_irrel"), sep = " ") %>% 
    select(-R_irrel, -numbers2) %>% 
    mutate_all(funs(str_replace(., pattern = "\\(", replacement = ""))) %>% 
    mutate_all(funs(str_replace(., pattern = "\\)", replacement = ""))) %>% 
    mutate_all(funs(str_replace(., pattern = ",", replacement = ""))) %>% 
    mutate_all(funs(str_replace(., pattern = "<", replacement = ""))) %>%
    mutate_all(funs(str_replace(., pattern = "−", replacement = "-"))) %>% 
    mutate(R1 = ifelse(var_name == "Groom divers.", 0.48, R1),
           CI_lower1 = ifelse(var_name == "Being approached", 0.14, CI_lower1)) %>% 
    mutate_at(vars(R1:p_val2), as.numeric) %>% 
    select(-F1, -F2)

# reshape a little:

longevity_adult <- AnAgeScrapeR::scrape_AnAge("Pan troglodytes", vars = "maximum_longevity_yrs", 
                   download_data = FALSE) 
lifespan <- as.numeric(longevity_adult$maximum_longevity_yrs) * 365
avg_adult_age <- lifespan * 0.25

# short  : 1 year
# long: 2.5 years

meta_table1 <- meta_table0 %>% select(var_name, R1:p_val1)
meta_table2 <- meta_table0 %>% select(var_name, R2:p_val2) %>% 
               rename(R1 = R2, CI_lower1 = CI_lower2, CI_upper1 = CI_upper2, p_val1 = p_val2)

# create final meta table
rbind(meta_table1, meta_table2) %>% 
    mutate(timeframe = rep(c("long", "short"), each = 15),
           sample_size = 18,
           measurements_per_ind = rep(c(2, 6), each = 15),
           t1 = avg_adult_age,
           t2 = rep(c(avg_adult_age + 365, avg_adult_age + 365*2.5), each = 15),
           delta_t = t2 - t1) %>% 
    rename(behaviour = 1, R = 2, CI_lower = 3, CI_upper  = 4, p_val = 5) %>% 
    mutate(Key = "INKLPDFJ",
         species_common = "chimpanzee",
         species_latin = "Pan_troglodytes",
        sex = 0,
        context = 1,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = NA,
        R_se = NA,
        remarks = "zoo setting",
        max_lifespan_days = lifespan) %>% 
    select(-timeframe) -> meta_table
    
write_delim(meta_table, path = "output/Koski_2011.txt", delim = " ", col_names = TRUE)



###### Paper 36: LeCoeur_Thibault_2015 ##########
# Le Coeur, Christie; Thibault, Martin; Pisanu, Benoit; Thibault, Sophie; Chapuis, Jean-Louis; Baudry, Emmanuelle	2015
# XMJ8ZGHT		
# temporally fluctuating selection on a personality trait in a wild rodent population

dat <- read_xlsx("data/downloaded_from_papers/lecoeur_thibault_2015.xlsx")

dat %>% group_by(ID) %>% mutate(count = n()) %>% filter(count == 2) %>% arrange(ID) %>% print(n = Inf)

# filter out cases where measurements are precisely one year apart
short_term <- dat %>% 
    group_by(ID) %>% 
    mutate(count = n()) %>% 
    arrange(ID) %>% 
    filter(count == 2) %>% 
    summarise(abs_diff = diff(range(yr))) %>% 
    filter(abs_diff == 1)

dat_short <- dat[dat$ID %in% short_term$ID, ]
# check
test1 <- dat_short %>% arrange(ID) %>% group_by(ID) %>% slice(1)
test2 <- dat_short %>% arrange(ID) %>% group_by(ID) %>% slice(2)

# log transform and gaussian like described in paper with year as random and ext as fixed
rpt_short_out <- dat_short %>% 
    mutate(trappy = log(trappy)) %>% 
    rptGaussian(trappy ~ (1|ID) + (1|yr) + ext, data = ., grname = "ID")


# longer term
longer_term <- dat %>% 
    group_by(ID) %>% 
    mutate(count = n()) %>% 
    arrange(ID) %>% 
    #filter(count != 2) %>% 
    summarize(abs_diff = diff(range(yr))) %>% 
    filter(abs_diff > 1)

dat_long <- dat[dat$ID %in% longer_term$ID, ]

# all measurements further than 1 year are on average 2.3 years apart, see :
dat_long %>% arrange(ID) %>% group_by(ID) %>% summarise(yr_range = diff(range(yr))) %>% summarise(mean(yr_range))

# measurments per ind in long term data
nrow(dat_long) / length(unique(dat_long$ID))

rpt_long_out <- dat_long %>% 
    mutate(trappy = log(trappy)) %>% 
    rptGaussian(trappy ~ (1|ID) + (1|yr) + ext, data = ., grname = "ID")


adult_lifespan <- scrape_AnAge("Tamias sibiricus", vars = "maximum_longevity_yrs", download_data = FALSE)
adult_lifespan <- as.numeric(adult_lifespan$maximum_longevity_yrs) * 365
avg_adult_age <- adult_lifespan * 0.25

# trappability ~ boldness
tribble(
                 ~R,                         ~R_se,                                       ~sample_size,  ~measurements_per_ind,  ~t1,                          ~t2,    ~delta_t,
    as.numeric(rpt_short_out$R), as.numeric(rpt_short_out$se),      as.numeric(rpt_short_out$ngroups),                      2, avg_adult_age, avg_adult_age + 365,       365,
    as.numeric(rpt_long_out$R),  as.numeric(rpt_long_out$se),       as.numeric(rpt_long_out$ngroups),                       3, avg_adult_age, avg_adult_age + 365*2.3,  365*2.3     
) %>% 
    mutate(Key = "XMJ8ZGHT",
           species_common = "siberian_chipmunk",
           species_latin = "Tamias_sibiricus",
           sex = 0,
           behaviour = "trappability_boldness",
           context = 3,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = c(NA, "between_years"),
           CI_lower = NA,
           CI_upper = NA,
           p_val = NA,
           remarks = NA,
           max_lifespan_days = adult_lifespan) -> meta_table

write_delim(meta_table, path = "output/LeCoeur_Thibault_2015.txt", delim = " ", col_names = TRUE)


    
###### Paper 37: Low_Makan_2012 #####
# Low, Matthew; Makan, Troy; Castro, Isabel	2012
# TJB2RJWP 
# food availability and offspring demand influence sex-specific patterns and repeatability of parental provisioning

# some info: up to 2 broods, 30 days parental care per brood
# between year: 1 year = 21 females / 16 males ,2 years = 8 females / 12 males and 3 years = 9 females / 8 males

# avg between year time per individual 2.53*365 for 17 females (also avg measurements_per_ind)
# females
(8*2 + 3*9) / 17
# males 2.4 * 365 for 20 males
(2*12 + 3*8) / 20

# within year =  30 days

# avg age stichbird
scrape_AnAge(latin_name = "Notiomystis cincta", vars =  "maximum_longevity_yrs")
# NA
# from http://www.tiritirimatangi.org.nz/stitchbird
# 7 years

lifespan_adult <- 7*365
avg_adult_age <- lifespan_adult * 0.25

# sample size (apparently) for short term:
# 1028 supplemented
# 102 nonsupplemented
# so roughly half half males/females

## digitise figure
# dat <- metaDigitise("data/to_digitise/study7_Low_Makan_2012/")
# write_delim(dat, "data/to_digitise/study7_Low_Makan_2012/digitised.txt")
dat <- read_delim(file = "data/to_digitise/study7_Low_Makan_2012/digitised.txt", delim = " ")
dat %>% 
    select(group_id, mean, se) %>% 
    mutate(sex = rep(c("female", "male"), 3),
           type_of_treatment = 1,
           treatment = c(rep("control", 2), rep("food_supplement", 4)),
           t1 = avg_adult_age,
           t2 = c(rep(avg_adult_age + 30, 4), avg_adult_age + 2.53*365, avg_adult_age + 2.4*365),
           delta_t = t2-t1,
           measurements_per_ind = c(2,2,2,2,2.5, 2.4),
           sex = c(1,2,1,2,1,2),
           behaviour = "parental_care",
           life_stage = "adult",
           event = c(rep(NA, 4), rep("between_years", 2)),
           CI_lower = NA,
           CI_upper = NA,
           remarks = "lifespan_from_website",
           max_lifespan_days = 7*365) %>% 
    select(-group_id) %>% 
    rename(R = mean, R_se = se) %>% 
    mutate(Key = "TJB2RJWP",
           species_common = "stitchbird",
           species_latin = "Notiomystis cincta",
           sample_size = c(51,51,514,514,17,20),
           context = 3,
           p_val = NA) -> meta_table
           
write_delim(meta_table, path = "output/Low_Makan_2012.txt", delim = " ", col_names = TRUE)



###### Paper 38: Martins_Schrama_2005 #####
# Martins, CIM; Schrama, JW; Verreth, JAJ	2005	
# the consistency of individual differences in growth, feed efficiency and feeding behaviour in african catfish clarias gariepinus (burchell 1822) housed individually
# AMJQ98EZ	

# african catfish
# 48
# juveniles: 
# maturity at 1.5 years (http://www.fisheriesjournal.com/archives/2017/vol5issue2/PartD/5-1-90-634.pdf)
# 1.5*365 * 0.5 = 274 days

avg_age_juvenile <- 274 
# lifespan from african catfish Bagrus bajad not Clarias gariepinus (also AnAge)
lifespan_adult <- 17*365 

# sample size 48
# P1: 0-15
# P2: 16-30
# P3: 31-47


scrape_AnAge(latin_name = "Clarias gariepinus", download_data=FALSE, vars = c( "male_maturity_days", "maximum_longevity_yrs"))

tribble(
  ~behaviour,      ~R,         ~R_se,                ~t1,      ~t2,  
"feed_intake",   0.70,         0.06,    avg_age_juvenile, avg_age_juvenile + 47, 
"feed_intake",   0.80,         0.06,    avg_age_juvenile, avg_age_juvenile + 30,
"residual_feed_intake", 0.49,  0.10,    avg_age_juvenile, avg_age_juvenile + 47,
"residual_feed_intake", 0.61,  0.11,    avg_age_juvenile, avg_age_juvenile + 47
) %>% 
    mutate(
        Key = " AMJQ98EZ",
        species_common = "african_catfish",
        species_latin = "Clarias_gariepinus",
        sample_size = 48,
        measurements_per_ind = c(2,2,3,3),
        sex = 0,
        context = 1,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "juvenile",
        event = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        delta_t = t2-t1,
        remarks = c("residual feeding is feeding relative to weight and growth", "lifespan taken from Bagrus bajad not Clarias gariepinus", NA, NA),
        max_lifespan_days = lifespan_adult
    ) -> meta_table

write_delim(meta_table, path = "output/Martins_Schrama_2005.txt", delim = " ", col_names = TRUE)


###### Paper 39: Milligan_Radersma_2017 #######
# Milligan, Nicole D.; Radersma, Reinder; Cole, Ella F.; Sheldon, Benjamin C.	2017
# to graze or gorge: consistency and flexibility of individual foraging tactics in tits
#	JYLBGHP7

dat <- read_csv("https://datadryad.org/bitstream/handle/10255/dryad.136951/parameter%20file%20for%20PCA.csv?sequence=1")
names(dat)
table(dat$winter)
table(dat$age)

# how many years are multiple years?
# 1020 IDs across multiple years
dat_mult <- dat %>% group_by(ring, winter) %>% tally() %>% count(ring) %>% arrange(ring) %>% filter(nn > 1) 
mean(as.numeric(dat_mult$nn, na.rm = TRUE)) # ~2.3
table(dat$age)
# measurements_per_ind short term:
meas_per_ind_all <- dat %>% group_by(ring) %>% tally() %>% filter(n > 1)
mean(meas_per_ind_all$n) # 8.5

meas_per_ind_long <- dat %>% filter(ring %in% dat_mult$ring) %>% group_by(ring) %>% count()
mean(meas_per_ind_long$n) # 13.6

delta_time <- 2.3*365

lv <- scrape_AnAge(latin_name = c("Parus major", "Cyanistes caeruleus"), "maximum_longevity_yrs" )
lv_long <-lv %>% mutate(max_days = as.numeric(maximum_longevity_yrs)*365,
                  avg_adult_age = max_days * 0.25)

# only five weekends (four weeks)
dat1 <- metaDigitise("data/to_digitise/study8_Milligan_Radersma_2017/")

dat1 %>% 
    as_tibble() %>% 
    select(group_id, mean, se, n) %>% 
    rename(R = mean,
           R_se = se, 
           sample_size = n) %>% 
    mutate(species_common = ifelse(str_detect(group_id, "gt"), "great_tit", "blue_tit"),
           species_latin = ifelse(str_detect(group_id, "gt"), "Parus major", "Cyanistes caeruleus"),
           behaviour = ifelse(str_detect(group_id, "binge"), "binge_eating_feeding_behaviour", "transience_feeding_behaviour"),
           delta_t = ifelse(str_detect(group_id, "multyrs"), delta_time, 30),
           t1 = case_when(
               species_common == "great_tit" ~ as.numeric(lv_long[2, "avg_adult_age"]),
               species_common == "blue_tit" ~ as.numeric(lv_long[1, "avg_adult_age"])
           ),
           t2 = 
            case_when(
                str_detect(group_id, "multyrs") ~ t1 + delta_time,
                !str_detect(group_id, "multyrs") ~ t1 + 30
            ),
           delta_t = t2-t1) %>% 
    select(-group_id) %>% 
    mutate(
        Key = "JYLBGHP7",
        measurements_per_ind = ifelse(delta_t == 839.5, 13.6, 8.5),
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        life_stage = "both",
        treatment = NA,
        event = ifelse(delta_t == 839.5, "between_years", NA),
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        remarks = NA,
        max_lifespan_days = ifelse(species_common == "great_tit", as.numeric(lv_long[2,3]), as.numeric(lv_long[1,3]))
    ) -> meta_table
    
write_delim(meta_table, path = "output/Milligan_Radersma_2017.txt", delim = " ", col_names = TRUE)



###### Paper 40: Mitchell_Fanson_2016 ####

# Mitchell, David J.; Fanson, Benjamin G.; Beckmann, Christa; Biro, Peter A.	2016	
# towards powerful experimental and statistical approaches to study intraindividual variability in labile traits
# 	XQZPRZ9N


# short time: 3 days, mean 5 measurements
# long : 3 weeks, mean 15 measurements
# 104 individuals, 1477 obs
# activity
# male guppies (Poecilia reticulata)

# Rshort-term =0.56 [0.48, 0.64]
#Rlong-term =0.31 [0.21, 0.41]

guppie_longevity <- scrape_AnAge(latin_name = "Poecilia reticulata", download_data = FALSE, vars = "maximum_longevity_yrs" )
guppie_long_days <- as.numeric(guppie_longevity$maximum_longevity_yrs) * 365
avg_adult_age <- guppie_long_days * 0.25

tribble(
    ~R,    ~CI_lower,     ~CI_upper,     ~sample_size,    ~measurements_per_ind,         ~t1,       ~t2,            ~delta_t,
0.56,      0.48,          0.64,                 104,                         15,  avg_adult_age, avg_adult_age+3,          3,
0.31,      0.21,          0.41,                 104,                         15,  avg_adult_age, avg_adult_age + 21,      21
 ) %>% 
    mutate(Key = "XQZPRZ9N",
           species_common = "guppy",
           species_latin = "Poecilia_reticulata",
           sex = 2,
           behaviour = "activity",
           context = 1, 
           type_of_treatment = NA,
           treatment = NA,
           life_stage = "adult",
           event = NA, 
           R_se = NA,
           p_val = NA,
           remarks = NA,
           max_lifespan_days = guppie_long_days) -> meta_table 

write_delim(meta_table, path = "output/Mitchell_Fanson_2016.txt", delim = " ", col_names = TRUE)



###### Paper 41: Murphy_Sexton_2008 #####
# Murphy, M.T.; Sexton, K.; Dolan, A.C.; Redmond, L.J.2008	
# dawn song of the eastern kingbird: an honest signal of male quality?
# 	Q3AM687B	

# eastern kingbird
# Tyrannus Tyrannus
# dawn songs
# singing recorded: mid-june to third week of july: 
# 40 individuals 2003, 88 observations
# 47 individuals with 107 observation in 2004
# 83 males
# mean time between observations: 11 days, 
# 173 song bouts with 1 obs n = 22, 2 obs N=22, three obs n = 10 four N = 10 five or more N = 6 over 2 year period

# number of observation between years:
72/16
# number of observations within 2003
62/23 # 2.7, i.e 11 * 2.7 = 30 days for within season
# 2004
82 / 35

max_long <- scrape_AnAge(latin_name = "Tyrannus tyrannus", vars = "maximum_longevity_yrs", download_data = FALSE)
max_long <- as.numeric(max_long$maximum_longevity_yrs) * 365
avg_adult_age <- max_long * 0.25

tribble(
   ~behaviour,        ~R,      ~p_val,        ~sample_size, ~measurements_per_ind,           ~t1,    ~t2,
"song_start_time", 0.205,  0.019,                    16,             4.5,          avg_adult_age,  avg_adult_age + 365,
"song_30min_rate", 0.195,  0.024,                    16,             4.5,          avg_adult_age,  avg_adult_age + 365,
"actual_song_rate",0.194,  0.025,                    16,             4.5,           avg_adult_age,  avg_adult_age + 365,
"peak_song_rate",  0.376,  0.000,                    16,             4.5,          avg_adult_age,  avg_adult_age + 365,
"song_end_time",   0.040,  0.306,                    16,             4.5,          avg_adult_age,  avg_adult_age + 365,
"song_bout_length",0.028,  0.354,                    16,             4.5,          avg_adult_age,  avg_adult_age + 365,

"song_start_time", 0.134,  0.167,                    23,             2.7,          avg_adult_age,  avg_adult_age + 30,
"song_30min_rate", 0.179,  0.102,                    23,             2.7,          avg_adult_age,  avg_adult_age + 30,
"actual_song_rate", 0.203,  0.076,                    23,             2.7,          avg_adult_age,  avg_adult_age + 30,
"peak_song_rate", 0.460,  0.001,                    23,             2.7,          avg_adult_age,  avg_adult_age + 30, 
"song_end_time", 0.423,  0.001,                    23,             2.7,          avg_adult_age,  avg_adult_age + 30,  
"song_bout_length", -0.062,  0.660,                    23,             2.7,          avg_adult_age,  avg_adult_age + 30,  
    
 "song_start_time", 0.092,  0.251,                    35,             2.3,          avg_adult_age,  avg_adult_age + 30,
"song_30min_rate", 0.226,  0.051,                    35,             2.3,          avg_adult_age,  avg_adult_age + 30,
"actual_song_rate", 0.374,  0.003,                    35,             2.3,          avg_adult_age,  avg_adult_age + 30,
"peak_song_rate", 0.385,  0.002,                    35,             2.3,          avg_adult_age,  avg_adult_age + 30, 
"song_end_time", 0.023,  0.428,                    35,             2.3,          avg_adult_age,  avg_adult_age + 30,  
"song_bout_length", 0.008,  0.470,                    35,             2.3,          avg_adult_age,  avg_adult_age + 30
) %>% 
    mutate(Key = "Q3AM687B",
        species_common = "eastern_kingbird",
        species_latin = "tyrannus_tyrannus",
        sex = 2,
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = c(rep("between_years", 6), rep("within_years", 12)),
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        delta_t = t2-t1,
        remarks = "no se",
        max_lifespan_days= max_long) -> meta_table


write_delim(meta_table, path = "output/Murphy_Sexton_2008.txt", delim = " ", col_names = TRUE)


###### Paper 42: Nelson_Wilson_2008 #####
# Nelson, Ximena J.; Wilson, David R.; Evans, Christopher S. 2008
# behavioral syndromes in stable social groups: an artifact of external constraints?
# 	A2UNW3XN	

# 36 male / 36 female golden sebright bantam chickens (Gallus Gallus domesticus)
# crowing: terrority
# alarm call: antipredator
# food call: foraging

dat <- tabulizer::extract_tables("data/papers/nelson2008.pdf", pages = 6)

# extract crowing
mat1 <- dat[[1]]
mat1 <- mat1[-1, -1]
mat1[upper.tri(mat1)] <- NA
mat1[mat1 == ""] <- NA

crowing_df <- tibble(day1 = rep(NA, 28), day2 = rep(NA, 28), R = rep(NA,28))

df_row <- 1
for (row_mat in 1:8) {
    for (col_mat in 1:8) {
        if (!is.na(mat1[row_mat, col_mat])) {
            crowing_df[df_row, ] <- c(col_mat, row_mat, mat1[row_mat, col_mat])
            df_row <- df_row + 1
        }
    }
}

# extract alarm calling
mat2 <- dat[[1]]
mat2 <- mat2[-1, -1]
mat2[lower.tri(mat2)] <- NA
mat2[mat2 == ""] <- NA

alarmcalling_df <-tibble(day1 = rep(NA, 28), day2 = rep(NA, 28), R = rep(NA,28))

df_row <- 1
for (row_mat in 1:8) {
    for (col_mat in 1:8) {
        if (!is.na(mat2[row_mat, col_mat])) {
            alarmcalling_df[df_row, ] <- c(col_mat, row_mat, mat2[row_mat, col_mat])
            df_row <- df_row + 1
        }
    }
}


# extract food calling
mat3 <- dat[[2]]
mat3 <- mat3[-1, -1]
mat3[upper.tri(mat3)] <- NA
mat3[mat3 == "1.0"] <- NA

foodcalling_df <- tibble(day1 = rep(NA, 28), day2 = rep(NA, 28), R = rep(NA,28))

df_row <- 1
for (row_mat in 1:8) {
    for (col_mat in 1:8) {
        if (!is.na(mat3[row_mat, col_mat])) {
            foodcalling_df[df_row, ] <- c(col_mat, row_mat, mat3[row_mat, col_mat])
            df_row <- df_row + 1
        }
    }
}

# scrape_AnAge(latin_name = "Gallus gallus", vars =  "maximum_longevity_yrs", download_data = TRUE)
# assumed longevity: 10 years (https://en.wikipedia.org/wiki/Chicken)
# so

longevity <- 10*365
avg_adult_age <- 0.25 * longevity

# create_meta_table
rbind(crowing_df, alarmcalling_df, foodcalling_df) %>% 
                 mutate(behaviour = c(rep("crowing_territory_behav", 28), rep("alarm_calling_antipredator_behav", 28), 
                                      rep("food_calling_foraging_behav", 28))) %>%
                 mutate(p_val = case_when(
                     #str_detect(R, "\\*{1}")  ~ 0.05,
                     str_detect(R, "\\*{2}") ~ 0.01
                 )) %>% 
                mutate(p_val = case_when(
                    str_detect(R, "\\*{1}") & is.na(p_val) ~ 0.05,
                    TRUE ~ p_val
                )) %>% 
                mutate(p_val = ifelse(is.na(p_val), ">0.05", p_val),
                       sample_size = 72,
                       t1 = avg_adult_age + as.numeric(day1),
                       t2 = avg_adult_age + as.numeric(day2),
                       delta_t = t2 -t1) %>% 
                select(-day1, -day2) %>% 
                mutate(
                    Key = "A2UNW3XN",
                    species_common = "Golden_sebright_bantam_chicken",
                    species_latin = "Gallus_gallus_domesticus",
                    measurements_per_ind = 2,
                    sex = 0,
                    context = 1,
                    type_of_treatment = 0,
                    treatment = NA,
                    life_stage = "adult",
                    event = NA,
                    R_se = NA,
                    CI_lower = NA,
                    CI_upper = NA,
                    remarks = "no se",
                    max_lifespan_days = longevity
                ) -> meta_table
                

write_delim(meta_table, path = "output/Nelson_Wilson_2008.txt", delim = " ", col_names = TRUE)



###### Paper 43: Nemiroff_Deslpland_2007 ####
# Nemiroff, L.; , Despl; , E.	2007	
# consistent individual differences in the foraging behaviour of forest tent caterpillars (malacosoma disstria)
# 	9K5MKM4T

scrape_AnAge(latin_name = "Malacosoma disstria", vars = "maximum_longevity_yrs", download_data = FALSE)

# experiment started after second instar, ~15 days
# maximum lifespan around 3 month (https://en.wikipedia.org/wiki/Tent_caterpillar)

tribble(
    ~behaviour,                ~R,     ~p_val,    ~t1,     ~t2,
    "latency_to_reach_food",   0.27,     0.05,       1,      2,
    "latency_to_reach_food",   0.25,     0.05,       1,      3,
    "latency_to_reach_food",   0.21,  ">0.05",       1,      4,
    "latency_to_reach_food",   0.33,     0.01,       2,      3,
    "latency_to_reach_food",   0.20,  ">0.05",       2,      4,
    "latency_to_reach_food",   0.30,     0.01,       3,      4,
    
    "No_movement",             0.27,     0.05,       1,      2,
    "No_movement",             0.38,     0.01,       1,      3,
    "No_movement",             0.10,  ">0.05",       1,      4,
    "No_movement",             0.37,     0.01,       2,      3,
    "No_movement",             0.15,  ">0.05",       2,      4,
    "No_movement",             0.42,    0.001,       3,      4,
    
    "walking",                 0.19,  ">0.05",       1,      2,
    "walking",                 0.50,    0.001,       1,      3,
    "walking",                 0.24,     0.05,       1,      4,
    "walking",                 0.30,     0.01,       2,      3,
    "walking",                -0.04,  ">0.05",       2,      4,
    "walking",                 0.36,    0.001,       3,      4,
    
    "searching",               0.49,    0.001,       1,      2,
    "searching",               0.53,    0.001,       1,      3,
    "searching",               0.07,  ">0.05",       1,      4,
    "searching",               0.48,    0.001,       2,      3,
    "searching",               0.20,  ">0.05",       2,      4,
    "searching",               0.21,     0.05,       3,      4,
    
    "Eating",               0.49,    0.001,       1,      2,
    "Eating",                0.53,    0.001,       1,      3,
    "Eating",                0.07,  ">0.05",       1,      4,
    "Eating",                0.48,    0.001,       2,      3,
    "Eating",                0.20,  ">0.05",       2,      4,
    "Eating",                0.21,     0.05,       3,      4
    
  ) %>% 
    mutate(
        Key = "9K5MKM4T",
        species_common = "forest_tent_caterpillar",
        species_latin = "Malacosoma_disstria",
        t1 = 15 + t1,
        t2 = 15 + t2,
        delta_t = t2 - t1,
        sample_size = 78,
        context = 1,
        type_of_treatment = 0,
        measurements_per_ind = 4,
        sex = 0,
        treatment = 0,
        life_stage = "juvenile",
        event = NA,
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        remarks = "No se, spearman correlation, testing started after 2nd instar, overall 5 instars",
        max_lifespan_days = 90
    ) -> meta_table
    
write_delim(meta_table, path = "output/Nemiroff_Deslpland_2007.txt", delim = " ", col_names = TRUE)



###### Paper 44: Niemel_Vainikka_2012 #####
###### Paper 45: Olsen_Heupel_2012 ####
# Olsen, Esben Mol; Heupel, Michelle R.; Simpfendorfer, Colin A.; , Mol; , Even		2012	
# harvest selection on atlantic cod behavioral traits: implications for spatial management
# AS8P9MSC

# atlantic cod
#  all adults > 30cm long, maturity at 2-4 years
cod_data <- scrape_AnAge(latin_name = "Gadus morhua", vars = "maximum_longevity_yrs", download_data = FALSE)
max_longevity <- as.numeric(cod_data$maximum_longevity_yrs) * 365
avg_adult_age <- max_longevity * 0.25

# 3 month period
# 

dat <- metaDigitise("data/to_digitise/study9_Olsen_Heupel_2012/")
dat %>% 
    as_tibble() %>% 
    select(group_id, mean, se) %>% 
    separate(group_id, into = c("behaviour", "time"), sep = "_") %>% 
    mutate(behaviour = str_replace_all(behaviour, "-", "_")) %>% 
    mutate(t1 = avg_adult_age,
           t2 = rep(c(avg_adult_age + 30, avg_adult_age + 60, avg_adult_age + 90), 5),
           delta_t = t2 - t1, 
           sample_size = 28,
           measurements_per_ind = 2) %>% 
    select(-time) %>% 
    rename(R = mean, R_se = se) %>% 
    mutate(
        Key = "AS8P9MSC",
        species_common = "atlantic_cod",
        species_latin = "Gadus_morhua",
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        treatment = "surgical transmitter implantation for all study individuals",
        life_stage = "adult",
        event = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        max_lifespan_days = max_longevity,
        remarks = "study on natural fish movement using tagging"
    ) -> meta_table

write_delim(meta_table, path = "output/Olsen_Heupel_2012.txt", delim = " ", col_names = TRUE)



###### Paper 46: Oufiero_Garl_2009 ####
# Oufiero, C.E.; , Garl; Jr., T. 2009	
# repeatability and correlation of swimming performances and size over varying time-scales in the guppy (poecilia reticulata)
# 	TWIX36XY	

# critical swimming speed
# burst_speed

# fish at 90 days of age
fish_age_start <- 90
fish_max_age <- scrape_AnAge(latin_name = "Poecilia reticulata", vars = "maximum_longevity_yrs", download_data = FALSE)
fish_max_age <- as.numeric(fish_max_age$maximum_longevity_yrs) * 365

tribble(
    ~behaviour,              ~R,       ~p_val,        ~t1,        ~t2,       ~sample_size,    ~measurements_per_ind,
    "burst_speed",          0.14,       0.1622,       7,           7.5,                22,                        3,
    "burst_speed",          0.25,       0.027,        36,          36.5,               22,                        3,
    "burst_speed",          0.38,       0.0027,       42,          42.5,               22,                        3,
    "burst_speed",          0.27,       0.06,         510,         510.5,              22,                        3,
    "burst_speed",          0.46,       0.03,         7,            36,                22,                        6,
    "burst_speed",          0.76,       0.0001,       36,           42,                22,                        6,
    "burst_speed",          0.27,         0.4,        42,           510,                22,                       6,
    
    "swimming_speed_max",  0.35,         0.09,         7,           36,                 22,                       2,
    "swimming_speed_max",  0.54,         0.006,        36,          42,                 22,                       2,
    "swimming_speed_max",  0.42,         0.048,        42,         510,                  22,                      2,
    
    "swimming_speed_critical", 0.74,     0.001,        7,           36,                 22,                       2,
    "swimming_speed_critical", 0.58,     0.003,       36,           42,                 22,                       2,
    "swimming_speed_critical", 0.34,     0.118,       42,          510 ,               22,                        2
) %>% 
    mutate(t1 = t1 + fish_age_start,
           t2 = t2 + fish_age_start,
           delta_t = t2-t1,
          Key = "TWIX36XY",
          species_common = "trinidadian_guppy",
          species_latin = "Poecilia_reticulata",
          sex = 0,
          context = 2,
          type_of_treatment = NA,
          treatment =NA,
          life_stage = "adult",
          event = NA,
          R_se = NA,
          CI_lower = NA,
          CI_upper = NA,
          remarks = "no se",
          max_lifespan_days = fish_max_age) -> meta_table
    
write_delim(meta_table, path = "output/Oufiero_Garl_2009.txt", delim = " ", col_names = TRUE)


###### Paper 47: Petelle_McCoy_2013 #####
# Petelle, Matthew B.; McCoy, Dakota E.; , Alej; ro, Vanessa; Martin, Julien G. A.; Blumstein, Daniel T.	2013
# development of boldness and docility in yellow-bellied marmots
# 	CSPH9Z5D


# juveniles: first summer
# yearlings: second summer of life: 365:365+90

dat_marmot <- scrape_AnAge("Marmota flaviventris", vars = "maximum_longevity_yrs", download_data = FALSE)
lifespan_marmot <- as.numeric(dat_marmot$maximum_longevity_yrs) * 365
avg_adult_age <- lifespan_marmot * 0.25
# very few individuals measured at different ages so adult measuring time assumed to be 90 (one summer) too

# measurements per ind: 
# boldness 563 / 237 = 2.38
# docility 8217 / (861 + 445 + 266) = 5.23 // trapping

tribble(
    ~behaviour,     ~R,     ~R_se,      ~t1,       ~t2,              ~life_stage,   ~sample_size,
    "docility",     0.19,   0.021,        0,         90,              "juvenile",    86,
    "docility",     0.27,   0.033,       365,       455,               "yearling",   81,
    "docility",     0.29,   0.039, avg_adult_age, avg_adult_age+90,    "adult",      70,
    
      
    "boldness",     0.019,   0.027,     0,         90,                  "juvenile",  861, 
    "boldness",     0.373,   0.146,     365,       455,                 "yearling",  445,
    "boldness",     0.053,   0.057,     avg_adult_age, avg_adult_age+90,    "adult",  266
) %>% 
    mutate(
        measurements_per_ind = c(rep(5.23, 3), rep(2.38, 3)),
        Key = "CSPH9Z5D",
        species_common = "yellow_bellied_marmot",
        species_latin = "Marmota_flaviventris",
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        event = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        delta_t = t2 - t1,
        remarks = "no longer term measurement but across the lifespan",
        max_lifespan_days = lifespan_marmot) -> meta_table


write_delim(meta_table, path = "output/Petelle_McCoy_2013.txt", delim = " ", col_names = TRUE)



###### Paper 48: Polverina_Cigliano_2016 #######
# Polverino, Giovanni; Cigliano, Claudia; Nakayama, Shinnosuke; Mehner, Thomas		2016
# emergence and development of personality over the ontogeny of fish in absence of environmental stress factors
# AIP8SGJE

# 30 days old
# subadult stage event: morphogenesis of anal fin
# mosquitofish
# all within-stage tests 7 days apart
# 30 days old / 7 days juvenile test / 67 days: 7 days subadult tests / 134 days: 7 days adult experiment: final age: 141
# mosquitofish lifespan : 1 year / Haake and Dean 1983 cited in paper

tribble(
    ~behaviour,              ~R,      ~CI_lower,     ~CI_upper,      ~p_val,     ~t1,     ~t2,   ~event,
    "emergence_latency",   0.11,            NA,              NA,      0.03,      30,      141,  "ontogeny",
    "hiding_time",         0.15,            NA,              NA,      0.01,      30,      141,  "ontogeny",
    "distance_moved",      0.25,            NA ,             NA,      0.01,      30,      141,  "ontogeny",
    "freezing_time",       0.20,            NA,              NA,      0.01,      30,      141,  "ontogeny",
    
    "emergence_latency",   0.06,            0,              0.30,     0.99,      30,      37,   NA,
    "hiding_time",         0.08,            0,              0.35,     0.64,      30,      37,   NA,
    "distance_moved",      0.06,            0 ,             0.31,     0.87,      30,      37,   NA,
    "freezing_time",       0.02,            0,              0.17,     0.99,      30,      37,   NA,
    
    "emergence_latency",   0.08,            0,              0.33,     0.87,      67,      74,   NA,
    "hiding_time",         0.62,          0.36,             0.82,     0.01,      67,      74,   NA,
    "distance_moved",      0.21,            0 ,             0.50,     0.18,      67,      74,   NA,
    "freezing_time",       0.20,            0,              0.50,     0.11,      67,      74,   NA,

    "emergence_latency",   0.04,            0,              0.26,     0.87,      134,     141,   NA,
    "hiding_time",         0.44,          0.09,             0.72,     0.01,      134,     141,   NA,
    "distance_moved",      0.48,            0.13,           0.74,     0.01,      134,     141,   NA,
    "freezing_time",       0.55,           0.24,            0.78,     0.11,      134,     141,   NA
) %>% 
    mutate(
        Key = "AIP8SGJE",
        species_common = "eastern_mosquitofish",
        species_latin = "Gambusia_holbrooki",
        sample_size = 40,
        measurements_per_ind = c(rep(6, 4), rep(2, 12)),
        sex = 0,
        context = 1,
        type_of_treatment = 0,
        treatment = "add libitum food",
        life_stage = c(rep("across_ontogeny", 4), rep("juvenile", 4), rep("subadult", 4), rep("adult", 4)),
        R_se = NA,
        delta_t = t2 - t1,
        remarks = "R_se for ontogeny taking as mean R_se from across other stages",
        max_lifespan_days = 365
    ) %>% 
    mutate(R_se = (CI_upper - R) / 1.96) %>% 
    select(behaviour, R, R_se, everything()) %>% 
    # calculate ontogeny se as mean of other ses per behaviour
    group_by(behaviour) %>% 
    mutate(R_se = ifelse(is.na(R_se), mean(R_se, na.rm = TRUE), R_se)) -> meta_table

write_delim(meta_table, path = "output/Polverina_Cigliano_2016.txt", delim = " ", col_names = TRUE)


###### Paper 49: Redmond_Murphy_2009 #####
# Redmond, Lucas J.; Murphy, Michael T.; Dolan, Amy C.; Sexton, Karen 2009	
# parental investment theory and nest defense by eastern kingbirds
# 	KN8NUB8J	

# nest defense behaviour
# only for males within and between years
# between: r = 0.739, p = 0.004, n = 9


incub_to_nest <- 32
lv_dat <- scrape_AnAge(latin_name = "Tyrannus tyrannus", vars = "maximum_longevity_yrs", download_data = FALSE)
max_longevity <- as.numeric(lv_dat$maximum_longevity_yrs) * 365
avg_adult_age <- max_longevity * 0.25



tribble(
    ~R,        ~CI_lower,      ~CI_upper,   ~R_se,    ~sample_size,              ~t1,                        ~t2, 
0.284,         -0.200,         0.686,        0.21,           11,           avg_adult_age,  avg_adult_age + incub_to_nest,
0.687,          0.471,         0.828,        0.07,            44,     avg_adult_age + 365,  avg_adult_age + 365 + incub_to_nest, 
0.739,          NA,              NA,         0.21,           9,            avg_adult_age,  avg_adult_age + 365 + incub_to_nest
) %>% 
    mutate(
        Key = "KN8NUB8J",
        species_common = "eastern_kingbird",
        species_latin = "tyrannus_tyrannus",
        measurements_per_ind = 2,
        sex = 2,
        behaviour = "nest_defense",
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = c(rep("incubation_to_nestling",2), "between_years"),
        p_val = NA,
        delta_t = t2 - t1,
        remarks = NA,
        max_lifespan_days = max_longevity
    ) -> meta_table
    
write_delim(meta_table, path = "output/Redmond_Murphy_2009.txt", delim = " ", col_names = TRUE)



###### Paper 50: Rockwell_Gabriel_2012 #####
# Rockwell, Christina; Gabriel, Pia O.; Black, Jeffrey M.	2012
# bolder, older, and selective: factors of individual-specific foraging behaviors in steller's jays
# 	PHTV4Z4S

# steller's jay
# short: 19.12.2008-11.03.2009 = 82 days long
dmy(11032009) - dmy(19122008)
# second short : 16.10.2009 to 24.12.2009 = 69 days long
dmy(24122009) - dmy(16102009) 
# between winter seasons: 370 days
dmy(24122009) - dmy(19122008)

# bird age on average 1460 days
avg_adult_age <- 1460

max_life <- as.numeric(scrape_AnAge("Cyanocitta stelleri", vars = "maximum_longevity_yrs", download_data = FALSE)$maximum_longevity_yrs)*365



tribble(
             ~behaviour,  ~R,             ~t1,              ~t2,      ~sample_size, ~measurements_per_ind,
"feeder_sampling_actions", 0.35, avg_adult_age, avg_adult_age + 82,         63,                   4.9,
"feeder_sampling_actions", 0.41, avg_adult_age, avg_adult_age + 69,         57,                   7.5,
"feeder_sampling_actions", 0.38, avg_adult_age, avg_adult_age + 370,        57,                   12.4
) %>% 
    mutate(
        Key = "PHTV4Z4S",
        species_common = "Steller's_jay",
        species_latin = "Cyanocitta_stelleri",
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = c(NA, NA, "between_winter_seasons"),
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        delta_t = t2 - t1,
        remarks = "no se, nopval",
        max_lifespan_days = max_life
    ) -> meta_table 


write_delim(meta_table, path = "output/Rockwell_Gabriel_2012.txt", delim = " ", col_names = TRUE)


###### Paper 51: Sakai_2018  ######
# Sakai, Osamu	2018
# comparison of personality between juveniles and adults in clonal gecko species
# 	2UYSQ68E

# lifespan: Brown, S.G., Murphy-Walker, S., 1996. Behavioural interactions between a rare male phenotype and female unisexual Lepidodactylus lugubris. Herpetol. J. 6, 69–73
# 5 years
max_life_gecko <- 5*365

# all delta ts within 3 days, 3 measuremnts
# 41 adult, 44 juvenile
# Lepidodactylus lugubris
# mourning gecko

scrape_AnAge(latin_name = "Lepidodactylus lugubris", vars = "maximum_longevity_yrs", download_data = FALSE)

avg_adult_age <- 0.25 * 10 * 365    # maximum of 10 years, non-official source www.joshsfrogs.com/mourning-gecko-lepidodactylus-lugubris.html
avg_juvenile_age <- 0.5 * 300       # again joshsfrogs as source. ...

tribble(
 ~behaviour,    ~R,  ~p_val,    ~sample_size,   ~measurements_per_ind,     ~t1,                            ~t2,     ~life_stage,
"exploration",  0.194, 0.013,     44,                         3,            avg_juvenile_age, avg_juvenile_age + 3,  "juvenile",
"boldness",     0.436, 0.001,     44,                    3,            avg_juvenile_age, avg_juvenile_age + 3,  "juvenile",
"exploration",  0.277, 0.001,     41,                      3,             avg_adult_age,      avg_adult_age + 3,  "adult",
"boldness",     0.319, 0.001,     41,                       3,            avg_adult_age,      avg_adult_age + 3,  "adult"
) %>% 
    mutate(
        Key = "2UYSQ68E",
        species_common = "mourning_gecko",
        species_latin = "Lepidodactylus_lugubris",
        sex = 1,
        context = 2,
        type_of_treatment = 0,
        treatment = NA,
        event = NA,
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        delta_t = t2 - t1,
        remarks = "no se, all clonal females, no scientific source for longevity so used breeder website for estimate",
        max_lifespan_days = max_life_gecko
    ) -> meta_table
write_delim(meta_table, path = "output/Sakai_2018.txt", delim = " ", col_names = TRUE)


###### Paper 52: Schuster_Carl_2017 ######
# Schuster, Andrea C.; Carl, Teresa; Foerster, Katharina	2017
# repeatability and consistency of individual behaviour in juvenile and adult eurasian harvest mice
# PTJDWIY8

max_lifespan <- 11*30 # in the wild according to paper

tribble(
    ~behaviour,            ~sample_size,   ~R,      ~R_se,        ~t1,        ~t2,       ~life_stage,  
    "activity_openfield",            38,  0.45,     0.13,       42,        49,            "juvenile",   
    "boldness",                      38,  0.27,     0.13,       42,        49,            "juvenile",
    "activity_ymaze",                39,  0.05,     0.11,       42,        49,            "juvenile",
    "exploration",                   39,  0.39,     0.14,       42,        49,            "juvenile",
    "spatial_recognition",           39,  0.00,       NA,        42,       49,            "juvenile",

    "activity_openfield",            31,  0.49,     0.17,       84,        168,           "adult",
    "boldness",                      31,  0.20,     0.15,       84,        168,           "adult",
    "activity_ymaze",                31,  0.11,     0.10,       84,        168,           "adult",
    "exploration",                   31,  0.00,     0.11,       84,        168,           "adult",
    "spatial_recognition",           31,  0.00,       NA,       84,        168,           "adult",
    
    "activity_openfield",            52,  0.47,     0.11,       42,        84,            "both",   
    "boldness",                      52,  0.36,     0.13,       42,        84,            "both",
    "activity_ymaze",                51,  0.07,     0.09,       42,        84,            "both",
    "exploration",                   52,  0.13,     0.11,       42,        84,            "both",
    "spatial_recognition",           51,  0.00,       NA,       42,        84,            "both",
    
    "activity_openfield",            16,  0.77,     0.13,       84,        168,            "adult",   
    "boldness",                      16,  0.53,     0.13,       84,        168,            "adult",
    "activity_ymaze",                16,  0.00,     0.11,       84,        168,            "adult",
    "exploration",                   16,  0.00,     0.14,       84,        168,            "adult",
    "spatial_recognition",           16,  0.00,       NA,       84,        168,            "adult"
) %>% 
    filter(behaviour != "spatial_recognition") %>% 
    mutate(event = c(rep(NA, 8), rep("maturation", 4), rep("sexual_contact", 4)),
           measurements_per_ind = 2,
           sex = 0,
           Key = "PTJDWIY8",
           species_common = "eurasian_harvest_mouse",
           species_latin = "Micromys_minutus",
           context = 1,
           p_val = NA,
           type_of_treatment = 1,
           treatment = c(rep(NA, 12), rep("sexual_contact", 4)),
           CI_lower =NA,
           CI_upper = NA,
           delta_t = t2 - t1,
           remarks = NA,
           max_lifespan_days = max_lifespan) -> meta_table

write_delim(meta_table, path = "output/Schuster_Carl_2017.txt", delim = " ", col_names = TRUE)


###### Paper 53: Seltman_Ose_2012 ####
# Seltmann, Martin W.; Ost, Markus; Jaatinen, Kim; Atkinson, Shannon; Mashburn, Kendall; Hollmen, Tuula	2012
# stress responsiveness, age and body condition interactively affect flight initiation distance in breeding female eiders
# 9BYQMQXW

# females
# eider ducks
# trapping: end of incubation period: may / early june
# measurements within season per female:  
# 2009: 2.4, N = 154;  #
# 2010: 2.7, N = 166, 
# 2011: 2.6, N = 168

# Between breeding seasons: 
# 2009-2010: 3.9, N = 101, 
# 2010-2011 4.4, N = 22
# 2009-2011: 5.7, N = 40

# longevity (according to paper) : 21 years
lifespan_eider <- 21 * 365
avg_adult_age <- lifespan_eider * 0.25 
# within season delta_t ~ 30 (end of incubation period in may / early june)

dat <- tabulizer::extract_areas("data/papers/Seltmann et al. - 2012 - stress responsiveness, age and body condition inte.pdf",
                                  pages = 4)

dat %>% .[[1]] %>% 
    .[-c(1), 1:6] %>% 
    as_tibble() %>%
    mutate(R_se = str_extract(V3, "^.{4}")) %>% 
    select(-V3, -V4, -V5) %>% 
    rename(year = V1, R = V2) %>% 
    separate(V6, into = c("obs", "sample_size"), sep = " ") %>% 
    mutate(sample_size = str_replace(sample_size, "\\(", ""),
        sample_size = str_replace(sample_size, "\\)", "")) %>% 
    mutate_at(vars(R:R_se), funs(as.numeric)) %>% 
    mutate(t1 = avg_adult_age,
           t2 = c(rep(avg_adult_age + 30, 3), rep(avg_adult_age + 365, 2), avg_adult_age + 730)) %>% 
    select(-year) %>% 
    mutate(measurements_per_ind = obs/sample_size) %>% 
    select(-obs) %>% 
    mutate(
        Key = "9BYQMQXW",
        species_common = "common_eider",
        species_latin = "Somateria_mollissima",
        sex = 1,
        behaviour = "flight_initiation_distance_FID",
        context = 3,
        type_of_treatment = NA,
        treatment = NA,
        life_stage = "adult",
        event = c(rep(NA, 3), rep("between_years", 3)),
        CI_lower = NA,
        CI_upper = NA,
        p_val = c(rep(0.001, 5),0.05),
        delta_t = t2 - t1,
        remarks = NA,
        max_lifespan_days = lifespan_eider
    ) -> meta_table

write_delim(meta_table, path = "output/Seltman_Ose_2012.txt", delim = " ", col_names = TRUE)




###### Paper 54: St-Hilaire_Reale_2017 #####
# St-Hilaire, Etienne; Reale, Denis; Garant, Dany 2017
# determinants, selection and heritability of docility in wild eastern chipmunks (tamias striatus)
# 	TBSUDC6R	

all_dat <- scrape_AnAge(latin_name = "Tamias striatus", vars = c("maximum_longevity_yrs", "male_maturity_days", "female_maturity_days"),
             download_data = FALSE)

avg_adult_age <- as.numeric(all_dat$maximum_longevity_yrs) * 365 * 0.25
avg_juvenile_age <- (as.numeric(all_dat$male_maturity_days) + as.numeric(all_dat$female_maturity_days)) / 2 * 0.5

# 4690 tests on 601 individuals for site 1
4690 / 601
# 1754 tests on 311 individuals on other sites
1754 / 311
# trapping in active period: april to october ~ 183 days

tribble(
    ~R,    ~CI_lower,  ~CI_upper,      ~sample_size,     ~measurements_per_ind,      ~life_stage,     ~t1,                  ~t2,                           ~remarks,        ~sex,          
 0.41,         0.36,      0.46,           3551,                 7.8,                  "adult",        avg_adult_age,         avg_adult_age + 183,            "study_site1",     0,
 0.31,         0.22,      0.40,           1139,                 7.8,                  "juvenile",     avg_juvenile_age,       avg_juvenile_age + 183,      "study_site1",     0,
 0.37,         0.28,      0.46,            908,                 5.6,                  "adult",        avg_adult_age,        avg_adult_age + 183,            "study_site2",     0,
 0.52,         0.43,      0.63,            846,                 5.6,                  "juvenile",     avg_juvenile_age,       avg_juvenile_age + 183,       "study_site2",     0,
    
 0.43,        0.37,       0.49,          2249,                  7.8,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site1",    2,
 0.46,        0.37,       0.57,           779,                  5.6,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site2",    2,
 0.32,        0.26,       0.38,          2441,                  7.8,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site1",    1,
 0.42,        0.33,       0.51,           975,                  5.6,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site2",    1,

 0.43,        0.37,       0.49,          2249,                  7.8,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site1",    2,
 0.46,        0.37,       0.57,           779,                  5.6,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site2",    2,
 0.32,        0.26,       0.38,          2441,                  7.8,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site1",    1,
 0.42,        0.33,       0.51,           975,                  5.6,                  "both",   avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site2",    1,

 0.40,       0.35,        0.45,          3358,                  7.8,                  "both",  avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site1_reproductive_season", 0, 
 0.47,       0.39,        0.53,          1448,                  5.6,                  "both",  avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site2_reproductive_season", 0, 
 0.32,       0.27,        0.41,          1332,                  7.8,                  "both",  avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site1_non_reproductive_season", 0, 
 0.36,       0.22,        0.54,          306,                   5.6,                  "both",  avg_adult_age-avg_juvenile_age, avg_adult_age-avg_juvenile_age + 183, "study_site2_non_reproductive_season", 0
) %>% 
    mutate(Key = "TBSUDC6R",
         species_common = "eastern_chipmunk",
         species_latin = "Tamias_striatus",
         behaviour = "docility",
         context = 3,
         type_of_treatment = 0,
         treatment = c(rep(NA, 12), rep("reproductive_season", 2), rep("non_reproductive_season", 2)),
         event = c(rep(NA, 12), rep("reproductive_season", 2), rep("non_reproductive_season", 2)),
         R_se = NA,
         p_val = NA,
         delta_t = t2 - t1,
         remarks = "same data, but adult vs juv, male vs fem and reprod. vs nonreprod. season",
         max_lifespan_days = as.numeric(all_dat$maximum_longevity_yrs) * 365) -> meta_table
         
write_delim(meta_table, path = "output/St-Hilaire_Reale_2017.txt", delim = " ", col_names = TRUE)




###### Paper 55: Svartberg_Tapper_2005 #####
# Svartberg, K; Tapper, I; Temrin, H; Radesater, T; Thorman, S 2005
# 	consistency of personality traits in dogs
# TLYD3J7U


# first to second test : 30 days
# second to third: 35 days
# age: 453 days

avg_age_dogs <- 453
dog <- scrape_AnAge(latin_name = "canis familiaris", vars =c("maximum_longevity_yrs", "male_maturity_days", "female_maturity_days"),
             download_data = FALSE)

lifespan_dog <- as.numeric(dog$maximum_longevity_yrs) * 365

tribble(
    ~behaviour,     ~R,    ~t1,      ~t2,
    "boldness",    .89,      1,       30,
    "boldness",    .77,      31,      66,
    "boldness",    .83,      1,       66,
    
    "playfulness",  .77,     1,       30,
    "playfulness",  .89,    31,       66,
    "playfulness",  .76,     1,       66,
    
    "chase_proneness",  .70,     1,       30,
    "chase_proneness",  .80,    31,       66,
    "chase_proneness",  .61,     1,       66,

    "curiosity_fearlessness",  .72,     1,       30,
    "curiosity_fearlessness",  .75,    31,       66,
    "curiosity_fearlessness",  .58,     1,       66,
    
    "sociability",  .72,     1,       30,
    "sociability",  .57,    31,       66,
    "sociability",  .57,     1,       66,
    
    "aggressiveness",  .72,     1,       30,
    "aggressiveness",  .57,    31,       66,
    "aggressiveness",  .57,     1,       66
) %>% 
    mutate(t1 = t1 + avg_age_dogs,
           t2 = t2 + avg_age_dogs,
           delta_t = t2 - t1,
           Key = "TLYD3J7U",
           species_common = "dog",
           species_latin = "Canis_familiaris",
           sample_size = 40,
           measurements_per_ind = 2,
           sex = 0,
           context = 1,
           type_of_treatment  = 0,
           treatment = NA,
           life_stage = "both",
           event = "maturation?",
           R_se = NA,
          CI_lower = NA,
          CI_upper = NA,
          p_val = 0.001,
          remarks = NA,
          max_lifespan_days = lifespan_dog) -> meta_table

    
write_delim(meta_table, path = "output/Svartberg_Tapper_2005.txt", delim = " ", col_names = TRUE)

###### Paper 56: Taylor_Cooke_2014 ####
# Taylor, M. K.; Cooke, S. J.	F3ESMMAM	2014	repeatability of movement behaviour in a wild salmonid revealed by telemetry


scrape_AnAge(latin_name = "Salvelinus confluentus", vars = "maximum_longevity_yrs", download_data = FALSE)
# longevity: 12 years ( Schemmel, E. & Dunham, J. Environ Biol Fish (2010) 89: 161. https://doi.org/10.1007/s10641-010-9708-8)

max_lifespan <- 12 * 365
avg_adult_age <- max_lifespan * 0.25

library(digitize)
dat <- digitize("data/to_digitise/study10_Taylor_Cooke_2014/Screenshot 2019-03-15 11.18.06.png")
dat2 <- digitize("data/to_digitise/study10_Taylor_Cooke_2014/Screenshot 2019-03-15 11.18.12.png")

dat_df <- dat %>% 
    rename(autumn = x, spring = y) %>% 
    gather(key = "time", value = "mean_movement") %>% 
    mutate(ID = rep(1:(n() / 2), 2)) 

dat_df2 <- dat2 %>% 
    rename(pm = x, am = y) %>% 
    gather(key = "time", value = "mean_movement") %>% 
    mutate(ID = rep(1:(n() / 2), 2)) 


# write_delim(rbind(dat_df, dat_df2), path = "data/to_digitise/study10_Taylor_Cooke_2014/digitised.txt")
 
rpt1 <- rptGaussian(mean_movement ~ (1|ID), data = dat_df, grname = "ID")
rpt2 <- rptGaussian(mean_movement ~ (1|ID), data = dat_df2, grname = "ID")

# season delta
mid_spring <- dmy(21042009) + (dmy(15062009) - dmy(21042009)) / 2
mid_autumn <- dmy(20102008) + (dmy(08122008) -  dmy(20102008)) / 2

delta_season <- ymd(mid_spring) - ymd(mid_autumn) # 187

# Salvelinus confluentus

tribble(
    ~R,                 ~R_se,                 ~t1,       ~t2,
    as.numeric(rpt1$R), as.numeric(rpt1$se),  avg_adult_age, avg_adult_age + as.numeric(delta_season),
    as.numeric(rpt2$R), as.numeric(rpt2$se),  avg_adult_age, avg_adult_age + 0.5
) %>% 
    mutate(Key = "F3ESMMAM",
          species_common = "bull_trout",
          species_latin = "Salvelinus_confluentus",
          sample_size = 17,
          measurements_per_ind = 2,
          sex = 0,
          behaviour = "movement_telemetry",
          context = 3,
          type_of_treatment = 0,
          treatment = NA,
          life_stage = "adult",
          event = c(NA, "autumn_to_spring"),
          CI_lower = NA,
          CI_upper = NA,
          delta_t = t2 - t1,
          remarks = NA,
          p_val = NA,
          max_lifespan_days = max_lifespan) -> meta_table
    
write_delim(meta_table, path = "output/Taylor_Cooke_2014.txt", delim = " ", col_names = TRUE)

    

###### Paper 57: Thys_Eens_2017 ####
# Thys, Bert; Eens, Marcel; Aerts, Silke; Delory, Am; , ine; Iserbyt, Arne; Pinxten, Rianne	2017
# exploration and sociability in a highly gregarious bird are repeatable across seasons and in the long term but are unrelated
# A5QEWIB6

dat <- extract_areas("data/papers/Thys et al. - 2017 - exploration and sociability in a highly gregarious.pdf",
                     pages = 7)

# starlings
# 30 juvenile males captured from the wild in October 2008, experiments started in spring 2011 and
# ended in spring 2013
dmy(06042011) - dmy(30102008)  # 888 days

starling_data <- scrape_AnAge(latin_name = "Sturnus vulgaris", vars = c("male_maturity_days", "maximum_longevity_yrs"),
             download_data = FALSE)

avg_age_juvenile <- as.numeric(starling_data$male_maturity_days) * 0.5
avg_age_thisstudy <- avg_age_juvenile + 888
max_lifespan_days <- as.numeric(starling_data$maximum_longevity_yrs) * 365
# calculating delta_ts
dmy(16042013) - dmy(06042011)
dmy(16042013) - dmy(27042011)

dat %>% 
    .[[1]] %>% 
    as_tibble() %>% 
    .[-c(1,2), ] %>% 
    filter(!(V1 %in% c("Overall", "Short term", "Across season" , "Across year"))) %>% 
    mutate(time = c(rep("overall", 4), "short_term", rep("across_season", 4), rep("across_year", 4))) %>% 
    select(V1, V4, time) %>% 
    separate(V4, into = c("R", "CI_lower", "CI_upper"), sep = " ") %>% 
    mutate(CI_lower = str_replace(CI_lower, "\\(", ""),
           CI_lower = str_replace(CI_lower, ";", ""),
           CI_upper = str_replace(CI_upper, "\\)", "")) %>% 
    rename(behaviour = V1) %>% 
    mutate(delta_t = c(741, 720, 720, 720, 11, 182, rep(174, 3), 724, rep(716, 3)),
           t1 = avg_age_thisstudy + c(0,7,7,7,0,0,7,7,7,0,7,7,7),
           t2 = t1 + delta_t) %>% 
    mutate(behaviour = case_when(
        str_detect(behaviour, "Expl") ~ "exploratory_behaviour",
        str_detect(behaviour, "TR") ~ "sociability_combined",
        str_detect(behaviour, "NB") ~ "sociability_nestbox",
        str_detect(behaviour, "FE") ~ "sociability_near_female"
    )) %>% 
    mutate(Key = "A5QEWIB6",
           species_common = "european_starling",
           species_latin = "Sturnus_vulgaris",
           sample_size = 30,
           measurements_per_ind = c(4,3,3,3,2,2,2,2,2,2,2,2,2),
           sex = 2,
           context = 2,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = time,
           R_se = NA,
           p_val = NA,
           remarks = "one sociability measure is the sum of the other two",
           max_lifespan_days = max_lifespan_days) %>% 
    select(-time) -> meta_table

write_delim(meta_table, path = "output/Thys_Eens_2017.txt", delim = " ", col_names = TRUE)



###### Paper 58: 







###### Paper 58: Twiss_Cairns_2012 #####
# Twiss, Sean D.; Cairns, Charlotte; Culloch, Ross M.; Richards, Shane A.; Pomeroy, Patrick P.		2012
# variation in female grey seal (halichoerus grypus) reproductive performance correlates to proactive-reactive behavioural types
# THVJK6B7

# data within season:
dmy(311009) - dmy(290909)   # 32 days
dmy(011110) - dmy(290910)    # 33 days

# RCV tests:
# 2009: N = 19 (2 test)
# 2010: N = 17
# 2009-2010: N = 7
# only adult females

dat_seal <- scrape_AnAge(latin_name = "halichoerus grypus", vars = "maximum_longevity_yrs", download_data = FALSE)
avg_adult_age <- as.numeric(dat_seal$maximum_longevity_yrs) * 365 * 0.25

avg_age_female <- 1

tribble(
    ~R,   ~CI_lower,    ~CI_upper,             ~t1,            ~t2,    ~sample_size,
    0.81,  0.58,          0.92,     avg_age_female, avg_age_female + 32, 19,
    0.70,  0.34,          0.88,     avg_age_female, avg_age_female + 33, 17,
    0.72,  0.011,         0.95,     avg_age_female, avg_age_female + 365, 7
) %>% 
    mutate(
        Key = "THVJK6B7",
        species_common = "grey_seal",
        species_latin = "halichoerus_grypus",
        measurements_per_ind = 2,
        sex = 1,
        behaviour = "pup_checking_by_mothers",
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        p_val = NA,
        event = NA,
        R_se = NA,
        delta_t = t2 - t1,
        remarks = "wolf_sound_during_all_trials",
        max_lifespan_days = as.numeric(dat_seal$maximum_longevity_yrs) * 365
    ) -> meta_table

write_delim(meta_table, path = "output/Twiss_Cairns_2012.txt", delim = " ", col_names = TRUE)




###### Paper 59: Wang_Brennan_2015 ####
# Wang, Mu-Yun; Brennan, Caroline H.; Lachlan, Robert F.; Chittka, Lars	2015
# speed-accuracy trade-offs and individually consistent decision making by individuals and dyads of zebrafish in a colour discrimination task
# KZDWYY3T

dat_day_1 <- digitize("data/to_digitise/study11_Wang_Brennan_2015/l7Q1lLBm.png")
#write_delim(dat_day_1, path = "data/to_digitise/study11_Wang_Brennan_2015/digitised.txt")

dat <- dat_day_1
dat %>% 
    mutate(day = rep(1:3, 17),
           ID = rep(1:17, each = 3)) %>% 
    rename(accuracy = x, decision_time = y) %>% 
    select(-accuracy)%>% 
    .[-c(13:15, 19:21), ] -> dat_final

# del: 13-15
# 19-21

# age around 360 days
avg_adult_age <- 360
# max lifespan?
lifespan <- scrape_AnAge("danio rerio", vars = "maximum_longevity_yrs",
        download_data = FALSE)


dat12 <- filter(dat_final, day != 3)
rpt12 <- rptGaussian(decision_time ~ (1|ID), grname = "ID", data = dat12)
dat23 <- filter(dat_final, day != 1)
rpt23 <- rptGaussian(decision_time ~ (1|ID), grname = "ID", data = dat23)
dat13 <- filter(dat_final, day != 2)
rpt13 <- rptGaussian(decision_time ~ (1|ID), grname = "ID", data = dat13)

tribble(
    ~R,           ~R_se,       ~t1,          ~t2,
0.397,           0.205,         1,             2,
0.442,           0.204,         2,             3,
0.705,           0.144,         1,             3
) %>% 
    mutate(
        behaviour = "decision_time",
        t1 = avg_adult_age+t1,
        t2 = avg_adult_age+t2,
        delta_t = t2-t1,
        Key = "KZDWYY3T",
        species_common = "zebrafish",
        species_latin = "danio_rerio",
        sample_size = 15,
        measurements_per_ind = 2,
        sex = 0,
        context = 1,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        remarks = NA,
        max_lifespan_days = as.numeric(lifespan$maximum_longevity_yrs) * 365
    ) -> meta_table

write_delim(meta_table, path = "output/Wang_Brennan_2015.txt", delim = " ", col_names = TRUE)


# zebrafish
# danio rerio

###### Paper 60: Watkins_1997 #####
# Watkins, TB 1997	
# the effect of metamorphosis on the repeatability of maximal locomotor performance in the pacific tree frog hyla regilla
# 	AZY734PE	

# timeline:
# jump_distance: 38 days after metamorphosis
# trials adult frogs within a few hours
# hyla regilla
# larval stage ~ 2 month
# stage 37: premetamorphic: 2 month - 15 days
# stage 42: start of metamorphosis around 2 month 
# from 42 to frog jump distance measurements ~ 2 month

scrape_AnAge(latin_name = "Pseudacris regilla", vars = "maximum_longevity_yrs")
# lifespan in the wild: 3 years
# https://pages.uoregon.edu/titus/herp_old/regillahistory.htm


tribble(
    ~R,     ~p_val,    ~sample_size,    ~t1,     ~t2,   ~event,    ~life_stage,
0.65,         0.01,     29,               45,     46,       NA,     "juvenile",
-0.029,       0.88,     29,               45,     60,    "hindlimps_grew", "juvenile",
-0.009,       -0.009,   29,               45,    105,    "metamorphosis",  "both",
-0.157,       0.42,     29,               60,    105,     "metamorphosis",  "both",
0.791,       0.001,     29,               105,   105.5,   NA,         "adult"
) %>% 
    mutate(
        Key = "AZY734PE",
        species_common = "pacific_tree_frog",
        species_latin = "hyla_regilla/Pseudacris_regilla",
        measurements_per_ind = c(rep(2, 4), 5),
        sex  = 0,
        behaviour = "burst_locomotor_performance",
        context = 2,
        type_of_treatment = 0,
        treatment = NA,
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        delta_t = t2 - t1,
        remarks = c(rep("pearson", 4), "kendall"),
        max_lifespan_days = 365 * 3
    ) -> meta_table

write_delim(meta_table, path = "output/Watkins_1997.txt", delim = " ", col_names = TRUE)



###### Paper 61: Wexler_Subach_2016 ####

#Wexler, Yonatan; Subach, Aziz; Pruitt, Jonathan N.; Scharf, Inon	2016
# behavioral repeatability of flour beetles before and after metamorphosis and throughout aging
# RKZ9H5QA


# Islam, W. (2017). Eco-Friendly Approaches for the Management of Red Flour Beetle: Tribolium castaneum (Herbst). Science Letters, 5(2), 105-114.
# life cycle:
# 9 days egg
# 7 weeks larvae
# 8 days pupae
# 2 years adult

# red flour beetle 
# larvae 2-4 days before population, t
# movement activity and adge prefewrence twice, with a 1-day interval
# adults 4 days after eclosion, twice with a 1-day interval
# adults 18 days after eclosion, twice with a 1 day interval
# 27 females / 29 males

# second set of beetles:
# tested after populaiton for 4 month, 9 pairs of 2 measurements
# 25 male / 25 female / starting from adult age onwards

age_larv <- 56
age_adult <- 56 + 14

tribble(
~behaviour,       ~R,      ~p_val,    ~t1,     ~t2,    ~sex,    ~life_stage, ~event,
"activity",   0.63,    0.004,    age_larv, age_larv + 1, 1,       "larva",    NA,
"activity",   0.53,    0.003,    age_larv, age_larv + 1, 2,       "larva",    NA,
"activity",   0.45,    0.019,    age_adult, age_adult + 1, 1,     "adult",    NA,
"activity",   0.77,    0.001,    age_adult, age_adult + 1, 2,     "adult",    NA,
"activity",  -0.15,    0.45,    age_larv,  age_adult,      1,     "across_ontogeny", "metamorphosis",
"activity",  -0.29,    0.13,    age_larv,  age_adult,      2,     "across_ontogeny","metamorphosis",
"activity",  0.61,     0.001,   age_adult, age_adult + 14, 1,     "adult", NA,
"activity",  0.59,     0.001,   age_adult, age_adult + 14, 2,     "adult", NA,   
    
"edge_preference",   0.44,    0.023,    age_larv, age_larv + 1, 1,    "larva",   NA,
"edge_preference",   0.49,    0.007,    age_larv, age_larv + 1, 2,    "larva",   NA,
"edge_preference",   0.48,    0.011,    age_adult, age_adult + 1, 1,  "adult",   NA,
"edge_preference",   0.32,    0.086,    age_adult, age_adult + 1, 2,  "adult",   NA,
"edge_preference",   0.10,    0.61,    age_larv,  age_adult,      1,  "across_ontogeny",   "metamorphosis",
"edge_preference",  -0.003,   0.99,    age_larv,  age_adult,      2,  "across_ontogeny",   "metamorphosis",
"edge_preference",  0.16,     0.43,   age_adult, age_adult + 14, 1,   "adult",  NA,
"edge_preference",  0.27,     0.18,   age_adult, age_adult + 14, 2,   "adult",  NA
             
) %>% 
    mutate(
        Key = "RKZ9H5QA",
        species_common = "red_flour_beetle",
        species_latin = "Tribolium_castaneum",
        sample_size = 58,
        measurements_per_ind = 2,
        sex = 0,
        context = 1,
        type_of_treatment = 0,
        treatment = NA,
        R_se = NA,
        p_val = NA,
        CI_lower = NA,
        CI_upper =NA,
        delta_t = t2 - t1,
        remarks = "some more data across the lifespan, but not sure about it, check with holger",
        max_lifespan_days = 787
    ) -> meta_table


write_delim(meta_table, path = "output/Wexler_Subach_2016.txt", delim = " ", col_names = TRUE)


###### Paper 62: White_Meekan_2015 #########
# White, James R.; Meekan, Mark G.; McCormick, Mark I.	2015	
# individual consistency in the behaviors of newly-settled reef fish
# ES4ET2BP

# Ambon damsel
# Pomacentrus amboinensis

# newly metamorphosed juveniles

# field:
# short term: 3 minutes, n = 18
# days: 9 observations over 3 days, n = 21
avg_juvenile_age <- 21

# scrape_AnAge(latin_name = "Pomacentrus amboinensis", vars = "maximum_longevity_yrs", download_data = FALSE)
adult_lifespan <- 6.5*365 # McCormick, M. I. (2016). Protogyny in a tropical damselfish: females queue for future benefit. PeerJ, 4, e2198.

tribble(
    ~behaviour,              ~R,    ~CI_lower,   ~CI_upper,   ~t1,                 ~t2,
    "bite_rate",            .64,      .39,       .83,        avg_juvenile_age,  avg_juvenile_age + 3/ (24*60),
    "distance_ventured",    .69,      .46,       .86,        avg_juvenile_age,  avg_juvenile_age + 3/ (24*60),
    "reef_height",          .52,      .24,       .76,        avg_juvenile_age,  avg_juvenile_age + 3/ (24*60),
    
    "bite_rate",            .77,      .64,       .88,        avg_juvenile_age,  avg_juvenile_age + 3,
    "distance_ventured",    .62,      .45,       .79,        avg_juvenile_age,  avg_juvenile_age + 3,
    "reef_height",          .33,      .16,       .55,        avg_juvenile_age,  avg_juvenile_age + 3
) %>% 
    mutate(Key = "ES4ET2BP",
        species_common = "ambon_damselfish",
        species_latin = "Pomacentrus_amboinensis",
        sample_size = rep(c(18,21), each = 3),
        measurements_per_ind = rep(c(2,9), each = 3),
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "juvenile",
        event = NA,
        R_se = NA,
        p_val = NA,
        delta_t = t2 - t1,
        remarks = NA,
        max_lifespan_days = adult_lifespan) -> meta_table

write_delim(meta_table, path = "output/White_Meekan_2015.txt", delim = " ", col_names = TRUE)

    
###### Paper 63: White_Briffa_2017 ####
# White, Stephen J.; Briffa, Mark	2017
# how do anthropogenic contaminants (acs) affect behaviour? multi-level analysis of the effects of copper on boldness in hermit crabs
# BSARAJJ5	

# females, males
# NN 17, 15
# NC 11, 21
# sample size n = 32 in each group
# 10 observations / 5 per period per individual

avg_adult_age <- 1

scrape_AnAge(latin_name = "Pagurus bernhardus", vars = "maximum_longevity_yrs", download_data = FALSE) # no

# 3-10 years / Bridger, D., Bonner, S. J., & Briffa, M. (2015). Individual quality and personality: bolder males are less fecund in the hermit crab Pagurus bernhardus. Proceedings of the Royal Society B: Biological Sciences, 282(1803), 20142492.

avg_adult_age <- 3.5*365 * 0.25

tribble(
    ~R,       ~CI_lower,        ~CI_upper,              ~t1,              ~t2,   ~type_of_treatment,     ~treatment, 
    .54,            .39,             .71,     avg_adult_age, avg_adult_age + 5,                 2,               NA,
    .70,            .56,             .81,     avg_adult_age + 8, avg_adult_age + 13,            2,               NA,
    
    .59,            .43,             .72,     avg_adult_age, avg_adult_age + 5,                 2,               "copper",
    .59,            .46,             .75,     avg_adult_age + 8, avg_adult_age + 13,            2,               "copper"
) %>% 
    mutate(behaviour = "boldness",
        species_common = "hermit_crab",
        species_latin = "Pagurus_bernhardus",
        t1 = t1 + avg_adult_age,
        t2 = t2 + avg_adult_age,
        delta_t = t2 - t1,
        sample_size = rep(32, 4),
        measurements_per_ind = 5,
        sex = 0,
        context = 2,
        life_stage = "adult",
        event = c(NA, NA, "copper", "copper"),
        R_se = NA,
        p_val = NA,
        remarks = NA,
        max_lifespan_days = 3.5*365,
        Key = "BSARAJJ5") -> meta_table

write_delim(meta_table, path = "output/White_Briffa_2017.txt", delim = " ", col_names = TRUE)





###### Paper 64: Winney_Schroeder_2018 #####
## Winney, I. S.; Schroeder, J.; Nakagawa, S.; Hsu, Y. -H.; Simons, M. J. P.; Sanchez-Tojar, A.; Mannarelli, M. -E.; Burke, T.	2018	
# DBYUTD2X
# heritability and social brood effects on personality in juvenile and adult life-history stages in a wild passerine

# breeding season april-august 
# boldness: 2011 - 2014 during breeding season delta_t ~5 month max
# exploration: 2010-2014, during non-breeding season delta_t ~ 3 month
# between year always just first measure in a season

# boldness: 90 inds more than once (between two and 9) (adults
# exploration: 117 between two and six times (adults
# activiity: 
# 

delta_t_breed <- 5*30
delta_t_nonbreed <- 3*30

library(tabulizer)

location <- "data/papers/Winney et al. - 2018 - heritability and social brood effects on personali.pdf"

R_table <- extract_areas(location, pages = 6)

# figure out average age

sparrow_dat <- scrape_AnAge("Passer domesticus", vars = "maximum_longevity_yrs", download_data = FALSE)
sparrow_lv <- as.numeric(sparrow_dat$maximum_longevity_yrs) * 365
avg_adult_age <- sparrow_lv * 0.25


R_table %>% .[[1]] %>% as_tibble() %>% 
    select(V1, V2, V3, V4, V11, V12) %>% 
    rename(group = V1, Nind = V2, sample_size = V3, Nobs = V4, R = V11, R_se = V12) %>% 
    filter(group != "" & Nind != "" & R_se != "") %>% 
    mutate_at(vars(Nind:R_se), as.numeric) %>% 
    mutate(measurements_per_ind = (Nobs - (Nind- sample_size)) / sample_size) %>% 
    select(-Nind, -Nobs) %>% 
    .[-c(1,2,3,8,9,10,15), ] %>% 
    mutate(behaviour = c(rep("boldness", 4), rep("exploration", 4), rep("activity", 3))) %>% 
    mutate(t1 = c(rep(avg_adult_age, 8), 10, 12, 10),
           t2 = case_when(
                str_detect(group, pattern = "Within") & (behaviour == "boldness") ~ t1 + delta_t_breed,
                str_detect(group, pattern = "Within") & (behaviour == "exploration") ~ t1 + delta_t_nonbreed,
                str_detect(group, pattern = "Across year") ~ t1 + (measurements_per_ind - 1) * 365,
                str_detect(group, pattern = "Across days") ~ t1 + 3,
                str_detect(group, pattern = "Day 10") ~ t1 + 0.5,
                str_detect(group, pattern = "Day 12") ~ t1 + 0.5)
    ) %>% 
    select(-group) %>% 
    mutate(
        Key = "DBYUTD2X",
        species_common = "house_sparrow",
        species_latin = "Passer domesticus",
        sex = 0,
        context = 3,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = c(rep("adult", 8), rep("juvenile", 3)),
        event = c(rep("within_year", 3), "between_years", rep("within_year", 3), "between_years", "within_day", "within_day", "between_days"),
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        delta_t = t2-t1,
        remarks = NA,
        max_lifespan_days = sparrow_lv
    ) -> meta_table
    

write_delim(meta_table, path = "output/Winney_Schroeder_2018.txt", delim = " ", col_names = TRUE)


    
###### Paper 65: Wuerz_Krueger_2015 ####
# Wuerz, Yvonne; Krueger, Oliver 2015	
# 3T3APIHX	
# personality over ontogeny in zebra finches: long-term repeatable traits but unstable behavioural syndromes

# 
# test days
t1 <- c(56, 103, 367)
t2 <- c(73, 121, 381)

filepath <- "data/papers/Wuerz and Krueger - 2015 - personality over ontogeny in zebra finches long-t.pdf"
R_table <- tabulizer::extract_areas(file = filepath, pages = 5)

R_table %>% .[[1]] %>% as_tibble() %>% 
    select(-V12, -V8, -V11) %>% 
    filter(V2 != "") %>% 
    mutate(behaviour = c("", rep(c("fearlessness", "exploration", "boldness", "activity", "aggression"), each = 3))) %>% 
    filter(behaviour != "") %>% 
   # separate(V4, into = c("CI_lower", "CI_upper"), sep = "-") %>% 
    separate(V7, into = c("CI", "p_val"), sep = " ", extra = "merge") %>% 
    separate(V10, into = c("CI2", "p_val2"),  sep = " ", extra = "merge") %>% 
    select(-V1) %>% 
    separate(V13, into = c("R4", "CI4", "p_val4"), sep = " ", extra = "merge") %>% 
    rename(sex = V2, R1 = V3, CI1 = V4, p_val1 = V5, R2 = V6, CI2 = CI, p_val2 = p_val,
           R3 = V9, CI3 = CI2, p_val3 = p_val2 ) -> meta_table_raw


meta_table_raw2 <- map(list(meta_table_raw[c(1:4, 14)], meta_table_raw[c(1, 5:7, 14)], meta_table_raw[c(1, 8:10, 14)],
               meta_table_raw[c(1, 11:13, 14)]), function(x) {
                   names(x) <- c("sex", "R", "CI", "p_val", "behaviour")
                   x
               }) %>% 
               bind_rows() %>% 
               mutate(timeperiod = rep(rep(c("subadult", "young_adult", "mature_adult", "long_interval"), each = 3), each = 5))

# lifespan from paper
maximum_lifespan_wild <- 5*365

# check
ggplot(meta_table_raw2, aes(timeperiod, as.numeric(R), fill = behaviour)) +
    geom_boxplot()

meta_table_raw2 %>% 
    separate(CI, into = c("CI_lower", "CI_upper"), sep = "-") %>% 
    mutate(t1 = c(rep(56, 15), rep(103, 15), rep(367, 15), rep(56, 15)),
           t2 = c(rep(73, 15), rep(121, 15), rep(381, 15), rep(381, 15))) %>% 
    mutate(sex = case_when(
        sex == "both" ~ 0,
        sex == "males" ~ 2,
        sex == "females" ~ 1
    )) %>% 
    mutate(
        Key = "3T3APIHX",
        species_common = "zebra_finch",
        species_latin = "Taeniopygia_guttata",
        sample_size = rep(c(52, 22, 30), 20),
        measurements_per_ind = 2,
        context = 1,
        type_of_treatment = 0,
        treatment = 0,
        event = c(rep(NA, 45), rep("maturation", 15)),
        R_se = NA,
        delta_t = t2 - t1,
        remarks = NA,
        max_lifespan_days = maximum_lifespan_wild,
        life_stage = "adult" # merging early and late adulthood here
    ) %>% 
    select(-timeperiod) -> meta_table


write_delim(meta_table, path = "output/Wuerz_Krueger_2015.txt", delim = " ", col_names = TRUE)



###### Paper 66: Zsebk_Herczeg_2017 ######
# Zsebk, S.; Herczeg, G.; Blzi, G.; Laczi, M.; Nagy, G.; Szsz, E.; Mark, G.; Trk, J.; Garamszegi, L.Z.	2017	
# short- and long-term repeatability and pseudo-repeatability of bird song: sensitivity of signals to varying environments
# IM7M86ND	

library(metaDigitise)

# sample sizes: 54 within-day, 29 between days, 16 between year
# male collared flycatchers (Ficedula albicollis)

# within day : 6 minutes
# between day: 2.92 days
# between years: 1.36 years

# all mature males
flycatcher_data <- scrape_AnAge("Ficedula albicollis", vars = "maximum_longevity_yrs", download_data = FALSE)
maximum_long_flycatcher <- as.numeric(flycatcher_data$maximum_longevity_yrs) * 365
avg_adult_age <- maximum_long_flycatcher * 0.25

dat <- metaDigitise("data/to_digitise/study13_Zsebk_Herczeg_2017/")
# write_delim(dat, "data/to_digitise/study13_Zsebk_Herczeg_2017/digitised.txt")
dat <- read_delim("data/to_digitise/study13_Zsebk_Herczeg_2017/digitised.txt", delim = " ")
dat %>% 
    as_tibble() %>% 
    select(group_id, mean, se) %>% 
    rename(R = mean, R_se = se) %>% 
    mutate(R = ifelse(R<0, 0, R)) %>% 
    mutate(timeperiod = rep(c("within_day", "between_days", "between_years"), 7)) %>% 
    mutate(group_id = str_replace(group_id, paste0("_", timeperiod), "")) %>% 
    mutate(t1 = avg_adult_age,
           t2 = case_when(
               timeperiod == "within_day" ~ avg_adult_age + 0.004,
               timeperiod == "between_days" ~ avg_adult_age + 2.92,
               timeperiod == "between_years" ~avg_adult_age + 1.36*365
           )) %>% 
    select(-timeperiod) %>% 
    rename(behaviour = group_id) %>% 
    mutate(
        Key = "IM7M86ND",
        species_common = "collared_flycatcher",
        species_latin = "ficedula_albicollis",
        sample_size = rep(c(54, 29, 16), 7),
        measurements_per_ind = 2,
        sex = 2,
        context = 3,
        type_of_treatment = NA,
        treatment = NA,
        life_stage = "adult",
        event = rep(c(NA, NA, "between_years"), 7),
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        delta_t = t2 - t1,
        remarks = "20 song recordings averaged per measurement",
        max_lifespan_days = maximum_long_flycatcher
    ) -> meta_table


write_delim(meta_table, path = "output/Zsebk_Herczeg_2017.txt", delim = " ", col_names = TRUE)

    


###### Paper 67: Amy_Ung_2017 ####
# Amy, Mathieu; Ung, Davy; Beguin, Nathalie; Leboucher, Gerard 2017	
# personality traits and behavioural profiles in the domestic canary are affected by sex and photoperiod
# 	GUF3YCJ4	

# this data I already extracted, so it's loaded here from an excel sheet

meta_table_raw <- read_xlsx(path = "data/amy_ung_2017.xlsx") %>% 
    select(species_common:remarks)

canary_longevity <- scrape_AnAge("Serinus canaria", vars = "maximum_longevity_yrs", download_data = FALSE)
canary_max_life <- as.numeric(canary_longevity$maximum_longevity_yrs) * 365


names(meta_table)[which(!(names(meta_table) %in% names(meta_table_raw) ))]

meta_table_raw %>% 
    mutate(measurements_per_ind = ifelse(delta_t == 190, 4, 2),
           Key = "GUF3YCJ4",
           species_latin = "Serinus_canaria",
           max_lifespan_days = canary_max_life) -> meta_table

write_delim(meta_table, path = "output/Amy_Ung_2017.txt", delim = " ", col_names = TRUE)

###### Paper 68: Boulton_Grimmer_2014 ####
# Boulton, Kay; Grimmer, Andrew J.; Rosenthal, Gil G.; Walling, Craig A.; Wilson, Alastair J.	2014	
# how stable are personalities? a multivariate view of behavioural variation over long and short timescales in the sheepshead swordtail, xiphophorus birchmanni
# PYZL2E8A

## 373 for long term
## 32 for short term
## Open field trial (OFT), emergence and exploration (EET) (both sort of boldness)
## LT: 30 weeks
## four OFT and four EET per fish
## ST: 4 days interval, 5 times

# mean longevity lab = 450 days reported
max_lifespan_fish <- 450 + 2*8.10 # mean and se

dat <- metaDigitise("data/to_digitise/study14_boulton_grimmer/")
#write_delim(dat, path = "data/to_digitise/study14_boulton_grimmer/digitised.txt")
scrape_AnAge(latin_name = "Xiphophorus birchmanni", vars = "maximum_longevity_yrs", download_data = FALSE)

dat %>% 
    as_tibble() %>% 
    select(group_id, mean, se) %>% 
    rename(R = mean, R_se = se) %>% 
    mutate(
        sample_size = rep(c(373, 32), 5),
        measurements_per_ind = rep(c(4,5), 5),
        sex = 0,
        behaviour = group_id,
        context = 1,
        type_of_treatment =0,
        treatment = NA,
        life_stage = "adult",
        event = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_val = NA,
        t1 = rep(c(203, 715), 5),
        t2 = rep(c(427, 732), 5),
        delta_t = t2-t1,
        remarks = c("open field test and emergence and exploration test, both should measure boldness"),
        max_lifespan_days = max_lifespan_fish,
        Key = "PYZL2E8A",
        species_common = "sheepshead_swordtail",
        species_latin = "Xiphophorus_birchmanni"
    ) %>% 
    select(-group_id) %>% 
    mutate(behaviour = str_remove(behaviour, pattern = "st_"),
           behaviour = str_remove(behaviour, pattern = "lt_"),
           behaviour = ifelse(behaviour == "emergence", paste0("eet_", behaviour), paste0("oft_", behaviour))) -> meta_table

write_delim(meta_table, path = "output/Boulton_Grimmer_2014.txt", delim = " ", col_names = TRUE)


###### Paper 69: English_Nakagawa_2010 ####
# English, S.; Nakagawa, S.; Clutton-Brock, T. H.	2010	
# consistent individual differences in cooperative behaviour in meerkats (suricata suricatta)
# 	XR32CK8A

# max_age = 4 years (paper)
# 2-19 measures per individual

dat <- tabulizer::extract_areas("data/papers/English et al. - 2010 - consistent individual differences in cooperative b.pdf",
                                pages = 4)
dat <- dat[[1]]
dat_df <- rbind(dat[,1:3], dat[,4:6])

# age classes: <1, 1-2, >2. avg_ages simply as minimum here and delta_t as spanning the whole period
dat_df %>% 
    as_tibble() %>% 
    separate(V1, into = c("R", "R_se"), sep = " \\(") %>% 
    mutate(R_se = str_replace(R_se, "\\)", "")) %>% 
    separate(V3, into = c("nobs", "sample_size")) %>% 
    rename(p_val = V2) %>% 
    mutate(behaviour = c(rep("babysitting", 3), rep("provisioning_pup", 3))) %>% 
    mutate(t1 = rep(c(182.5, 365, 730), 2),
           t2 = rep(c(365, 730, 1460), 2),
           delta_t = t2-t1) -> meta_table_raw

tribble(
    ~R,       ~R_se,     ~nobs,     ~sample_size,      ~behaviour,     ~t1,       ~t2,   
    0.218,     0.046,     896+948+738, 562+521+250,    "babysitting",  182.5,     1460,  
    0.513,     0.020,     896+948+738, 562+521+250,    "provisioning_pup", 182.5,  1460
) %>% 
    mutate(delta_t = t2-t1,
            p_val =  NA) -> meta_table_raw2

rbind(meta_table_raw, meta_table_raw2) %>% 
    mutate(Key = "XR32CK8A",
           species_common = "meerkat",
           species_latin = "Suricata_suricatta",
           measurements_per_ind = as.numeric(nobs)/as.numeric(sample_size),
           sex = 0,
           context = 3,
           type_of_treatment = NA,
           treatment = NA,
           life_stage = c("juvenile", "adult", "adult", "juvenile", "adult", "adult", "both", "both"),
           event = c(rep(NA, 6), "maturation_growth", "maturation_growth"),
           CI_lower = NA,
           CI_upper = NA,
           remarks = NA,
           max_lifespan_days = 4*365) %>% 
    select(-nobs) -> meta_table
           
write_delim(meta_table, path = "output/English_Nakagawa_2010 .txt", delim = " ", col_names = TRUE)

    

###### Paper 70: Hirata_Taketomi_2013 ####
# Hirata, Masahiko; Taketomi, Ikuko; Matsumoto, Yuka; Kubo, Shotaro	2013
# trade-offs between feeding and social companionship in cattle: intra-animal consistency over short and extended periods
# 76WVUJ9V

# 3 and 6 November 2010 (Experiment 1), 27–30
# November 2010 (Experiment 2) and 22–23 November 2011 (Experiment 3)
# 1,2,3: 8,16,12     1and2: 8, 1and3 7, 2and3: 12
#maximum (Dmax) and mean (Dmean) distance from the group, number of total (Ntotal) and different (Ndiff) tub visits, and
# proportion of time eating concentrate (Peatconc) and grazing sward (Pgraze)

# behaviours: willingness to trade sociability for feeding

# mean age see:
# T. Gotoh, H. Takahashi, T. Nishimura, K. Kuchida, H. Mannen, Meat produced by Japanese Black cattle and Wagyu, Animal Frontiers, Volume 4, Issue 4, October 2014, Pages 46–54, https://doi.org/10.2527/af.2014-0033
# slaughter age ~ 28-30 month
# minus 1 year ~18 max
# so probably around one year

avg_adult_age <- 365
# pval of 0.5 translated from NS in the paper
dat_bos <- scrape_AnAge("bos taurus", vars = "maximum_longevity_yrs", download_data = FALSE)
max_long_bos <- as.numeric(dat_bos$maximum_longevity_yrs) * 365

tribble(
    ~behaviour,              ~R,   ~p_val,   ~sample_size,   ~t1,   ~t2,
    "max_dist_from_group",  0.86,   0.01,             8, avg_adult_age, avg_adult_age + 3,  
    "mean_dist_from_group", 0.75,   0.05,             8, avg_adult_age, avg_adult_age + 3,
    "total_tubs_visited",   0.92,   0.01,             8, avg_adult_age, avg_adult_age + 3,
    "different_tubs_visited", 0.83, 0.05,             8, avg_adult_age, avg_adult_age + 3,
    "prop_eating_tub",      0.53,   0.5,              8, avg_adult_age, avg_adult_age + 3,
    "prop_grazing",         0.63,   0.5,              8, avg_adult_age, avg_adult_age + 3,
    
    "max_dist_from_group",  0.72,   0.01,             16, avg_adult_age, avg_adult_age + 3,  
    "mean_dist_from_group", 0.63,   0.01,             16, avg_adult_age, avg_adult_age + 3,
    "total_tubs_visited",   0.82,   0.001,             16, avg_adult_age, avg_adult_age + 3,
    "different_tubs_visited", 0.60, 0.05,             16, avg_adult_age, avg_adult_age + 3,
    "prop_eating_tub",      0.30,   0.5,              16, avg_adult_age, avg_adult_age + 3,
    "prop_grazing",         0.34,   0.5,              16, avg_adult_age, avg_adult_age + 3,
    
    "max_dist_from_group",  0.43,   0.5,             12, avg_adult_age, avg_adult_age +  2,  
    "mean_dist_from_group", 0.41,   0.5,             12, avg_adult_age, avg_adult_age +  2,
    "total_tubs_visited",   0.84,   0.001,            12, avg_adult_age, avg_adult_age + 2,
    "different_tubs_visited", 0.67, 0.05,             12, avg_adult_age, avg_adult_age + 2,
    "prop_eating_tub",      0.68,   0.05,              12, avg_adult_age, avg_adult_age +2,
    "prop_grazing",         0.59,   0.05,              12, avg_adult_age, avg_adult_age +2,
    
    "max_dist_from_group",  0.89,   0.01,             8, avg_adult_age, avg_adult_age + 21,  
    "mean_dist_from_group", 0.91,   0.01,             8, avg_adult_age, avg_adult_age + 21,
    "total_tubs_visited",   0.88,   0.01,             8, avg_adult_age, avg_adult_age + 21,
    "different_tubs_visited", 0.81, 0.05,             8, avg_adult_age, avg_adult_age + 21,
    "prop_eating_tub",      0.75,   0.05,              8, avg_adult_age, avg_adult_age + 21,
    "prop_grazing",         0.83,   0.05,              8, avg_adult_age, avg_adult_age + 21,
    
    "max_dist_from_group",  0.72,   0.01,            12,  avg_adult_age, avg_adult_age +  365,  
    "mean_dist_from_group", 0.57,   0.5,             12, avg_adult_age, avg_adult_age +  365,
    "total_tubs_visited",   0.68,   0.05,            12,  avg_adult_age, avg_adult_age +  365,
    "different_tubs_visited", 0.52, 0.5,             12, avg_adult_age, avg_adult_age +  365,
    "prop_eating_tub",      0.30,   0.5,             12,  avg_adult_age, avg_adult_age + 365,
    "prop_grazing",         0.31,   0.5,             12,  avg_adult_age, avg_adult_age + 365,

    "max_dist_from_group",  0.72,   0.05,            7, avg_adult_age, avg_adult_age + 386,  
    "mean_dist_from_group", 0.57,   0.5,             7,avg_adult_age, avg_adult_age +  386,
    "total_tubs_visited",   0.68,   0.01,            7, avg_adult_age, avg_adult_age + 386,
    "different_tubs_visited", 0.52, 0.05,             7,avg_adult_age, avg_adult_age +  386,
    "prop_eating_tub",      0.30,   0.5,             7, avg_adult_age, avg_adult_age + 386,
    "prop_grazing",         0.31,   0.5,             7, avg_adult_age, avg_adult_age + 386

) %>% 
    mutate(
        Key = "76WVUJ9V",
        species_common = "japanese_black_cattle",
        species_latin = "bos_taurus",
        measurements_per_ind = 2,
        sex = 1,
        context = 1,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = NA,
        R_se = NA,
        CI_lower =NA,
        CI_upper = NA,
        delta_t = t2 - t1,
        remarks = "all pearson correlation",
        max_lifespan_days = max_long_bos) -> meta_table
    

write_delim(meta_table, path = "output/Hirata_Taketomi_2013.txt", delim = " ", col_names = TRUE)



###### Paper 71: Spinka_Stehulova_2002 #####
# Spinka, M; Stehulova, I; Zacharova, J; Maletinska, J; Illmann, G	2002
# nursing behaviour and nursing vocalisations in domestic sows: repeatability and relationship with maternal investment
# P8R6KJHZ

# 2 lactations: 6 month
# 2-5: 6 month to 3 years, avg: 365+365*0.25 = 456.25


scrape_AnAge(latin_name = "Sus domesticus", vars = c("maximum_longevity_yrs", "female_maturity_days"), 
            download_data = FALSE)

# assumed average age: one and a half years
avg_adult_age <- 547.5
max_adult_age <- 3650 # needs citation, from wikipedia

tribble(
    ~behaviour,                                ~R,     ~p_val,    ~t1,                                      ~t2,  ~sample_size,               
    "behaviour_nursing_interval",             0.04,     "NS",  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_interval_nutritive",   0.09,     "NS",  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_prop_non_nutrative",   0.09,     "NS",  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_udder_massage_duration",   0.08, "NS",  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_udder_massage_duration2", 0.07,  "NS",  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_prop_terminated",       0.23,     "NS", avg_adult_age + 11,  avg_adult_age + 11 + 456.25,   11,
    "behaviour_nursing_prop_initiated",        0.28,     0.1,  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_preference_left_side",  0.27,     0.1,  avg_adult_age + 11, avg_adult_age +  11 + 456.25,   11,
    "behaviour_nursing_preference_one_side",   0.35,     0.05, avg_adult_age + 11,  avg_adult_age + 11 + 456.25,   11,
    "vocalisation_nursing_rate",               0.91,    0.001, avg_adult_age + 11,  avg_adult_age + 11 + 456.25,   11,
    "vocalisation_nursing_total",              0.88,    0.001, avg_adult_age + 11,  avg_adult_age + 11 + 456.25,   11,
    "vocalisation_grunt_rate_increase",        0.73,    0.001, avg_adult_age + 11,  avg_adult_age + 11 + 456.25,   11,
    
    "behaviour_nursing_interval",             0.05,     "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_interval_nutritive",   0.04,     "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_prop_non_nutrative",   0.01,     "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_udder_massage_duration",   0.23, "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_udder_massage_duration2", 0.37,  "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_prop_terminated",       0.40,     0.1, avg_adult_age + 7,  avg_adult_age + 7 + 182.5,   14,
    "behaviour_nursing_prop_initiated",        0.04,     "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_preference_left_side",  0.15,     "NS",  avg_adult_age + 7, avg_adult_age +  7 + 182.5,   14,
    "behaviour_nursing_preference_one_side",   0.03,     "NS", avg_adult_age + 7,  avg_adult_age + 7 + 182.5,   14,
    "vocalisation_nursing_rate",               0.80,    0.01, avg_adult_age + 7,  avg_adult_age + 7 + 182.5,   14,
    "vocalisation_nursing_total",              0.63,    0.01, avg_adult_age + 7,  avg_adult_age + 7 + 182.5,   14,
    "vocalisation_grunt_rate_increase",        0.66,    0.05, avg_adult_age + 7,  avg_adult_age + 7 + 182.5,   14,
    
    "behaviour_nursing_interval",             0.02,     "NS",  avg_adult_age  + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_interval_nutritive",   0.01,     "NS",  avg_adult_age  + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_prop_non_nutrative",   0.10,     "NS",  avg_adult_age  + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_udder_massage_duration",   0.24, "NS",  avg_adult_age  + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_udder_massage_duration2", 0.00,  "NS",  avg_adult_age  + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_prop_terminated",       0.04,     0.1, avg_adult_age  + 28,  avg_adult_age + 28 + 182.5,   14,
    "behaviour_nursing_prop_initiated",        0.09,     "NS",  avg_adult_age + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_preference_left_side",  0.58,     0.05,  avg_adult_age + 28, avg_adult_age + 28  + 182.5,   14,
    "behaviour_nursing_preference_one_side",   0.05,     "NS", avg_adult_age  + 28, avg_adult_age + 28 + 182.5,   14,
    "vocalisation_nursing_rate",               0.92,    0.001, avg_adult_age  + 28,  avg_adult_age + 28 + 182.5,   14,
    "vocalisation_nursing_total",              0.88,    0.001, avg_adult_age  + 28,  avg_adult_age + 28 + 182.5,   14,
    "vocalisation_grunt_rate_increase",        0.65,    0.1, avg_adult_age  + 28,  avg_adult_age + 28 + 182.5,   14,
    
    "behaviour_nursing_interval",             0.35,      0.1,  avg_adult_age  + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_interval_nutritive",   -0.26,     "NS",  avg_adult_age  + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_prop_non_nutrative",   -0.03,     "NS",  avg_adult_age  + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_udder_massage_duration",   0.24, "NS",  avg_adult_age  + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_udder_massage_duration2", -0.05,  "NS",  avg_adult_age  + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_prop_terminated",       0.06,     0.1, avg_adult_age   + 11,  avg_adult_age + 28  , 14,
    "behaviour_nursing_prop_initiated",        0.19,     "NS",  avg_adult_age + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_preference_left_side",  0.22,     0.05,  avg_adult_age + 11,  avg_adult_age + 28  ,  14,
    "behaviour_nursing_preference_one_side",   0.01,     "NS", avg_adult_age  + 11,  avg_adult_age + 28  , 14,
    "vocalisation_nursing_rate",               0.83,    0.001, avg_adult_age  + 11,   avg_adult_age + 28 ,  14,
    "vocalisation_nursing_total",              0.83,    0.001, avg_adult_age  + 11,   avg_adult_age + 28 ,  14,
    "vocalisation_grunt_rate_increase",        0.65,    0.01, avg_adult_age  + + 11,  avg_adult_age + 28  , 14
   
) %>% 
    mutate(Key = "P8R6KJHZ",
        species_common = "domestic_pig",
        species_latin = "Sus_domesticus",
        measurements_per_ind = 2,
        sex = 1,
        context = 1,
        type_of_treatment = 0,
        treatment = NA,
        life_stage = "adult",
        event = c(rep("between_lactations", 36), rep("within_lactation", 12)),
        R_se = NA,
        CI_lower = NA,
        CI_upper = NA,
        delta_t = t2 - t1,
        remarks = "pearson correlation",
        max_lifespan_days = max_adult_age) -> meta_table

write_delim(meta_table, path = "output/Spinka_Stehulova_2002.txt", delim = " ", col_names = TRUE)



### meta tables merged in different file

###### Paper 72: Goold_Newberry_2017 #####
# Goold, Conor; Newberry, Ruth C.	2017
# aggressiveness as a latent personality trait of domestic dogs: testing local independence and measurement invariance
# C993NZT8

url <- "https://github.com/ConorGoold/GooldNewberry_aggression_shelter_dogs/blob/master/raw_data.csv?raw=true"

dat <- read_csv(url)
names(dat)
hist(as.numeric(dat$mean_age))

dat$id  

dat[dat$id == 18, ]
dat %>% 
    select(site, neutered, total_days, source_type)

###### Paper 73: Ferrari_Millot_2015 #####
# Ferrari, Sebastien; Millot, S; , ie; Leguay, Didier; Chatain, Beatrice; Begout, Marie-Laure	2015
# consistency in european seabass coping styles: a life-history approach
# ZAKGTDQK	

seabass_dat <- scrape_AnAge("Dicentrarchus labrax", vars = "maximum_longevity_yrs")

seabass_max <- as.numeric(seabass_dat$maximum_longevity_yrs) * 365

tribble(
  ~behaviour,                          ~t1,            ~t2,         ~sample_size,   ~R,     ~p_val,
  "feeding_recovery_behaviour_pca",         129,       283,           29,           0.33,  0.09, 
  "feeding_recovery_behaviour_pca",         129,       548,           21,           0.05,  0.83,
  "feeding_recovery_behaviour_pca",         283,       548,           29,           0.12,  0.59,
  "restraint_test_escape_behaviour",        557,       739,           22,           0.19,  0.39,
  "restraint_test_escape_behaviour",        557,       758,           17,           0.36,  0.15,
  "restraint_test_escape_behaviour",        739,       758,           17,         -0.06,   0.8,
  "risk_taking",                            187,       202,           30,           0.49,  0.006,
  "risk_taking",                            187,       216,           30,           0.72,  0.001,
  "risk_taking",                            202,       216,           30,            0.33,  0.06
) %>% 
  mutate(Key = "ZAKGTDQK",
         species_common = "seabass",
         species_latin = "Dicentrarchus_labrax",
         measurements_per_ind = 2,
         sex = 0,
         context = 1,
         type_of_treatment = 0,
         treatment = NA,
         life_stage = "juvenile",
         event = NA,
         R_se = NA,
         CI_lower = NA,
         CI_upper = NA,
         delta_t = t2-t1,
         remarks = "spearman correlation",
         max_lifespan_days = seabass_max) -> meta_table
         
write_delim(meta_table, path = "output/Ferrari_Millot_2015.txt", delim = " ", col_names = TRUE)


###### Paper 74: Araya_Ajoy_2017 ######
# Araya-Ajoy, Yimen G.; Dingemanse, Niels J. 2017
# 	2UJ6RD3M
# repeatability, heritability, and age-dependence of seasonal plasticity in aggressiveness in a wild passerine bird

# nest stage 1: one model, nest stage 2: another model
# random: brood id, male id, nest box

# : 5 week
# zwischen Bruten innerhalb von Jahr: Mai/Juni 6 * 7 Tage

# 2 measurements in egg laying phase
# 2 measurements in incubation
# ~7 days

# between years ~ 365 days
library(AnAgeScrapeR)
tit_dat <- scrape_AnAge("Parus major", vars = "maximum_longevity_yrs")
tit_max <- as.numeric(tit_dat$maximum_longevity_yrs) * 365

library(rptR)
md <- read_csv("data/downloaded_from_papers/dingemanse_great_tit.csv")
md %>% 
  group_by(MaleID, AggressionYear) %>% 
  filter( NestStage==1) %>% 
  count()

# short term repeatability with broodID (only one brood per year)
mod_st_1 <- rptGaussian(TMaleMinDistance ~ (1|BroodID), grname = "BroodID", data = subset(md, NestStage==1))
mod_st_2 <- rptGaussian(TMaleMinDistance ~ (1|BroodID), grname = "BroodID", data = subset(md, NestStage==2))

# age / sample size short term
age_df <- md %>% 
  group_by(NestStage) %>% 
  summarise(age = mean(minAge) * 365)
md %>% 
  group_by(BroodID) %>% 
  count()

# long term
mod_lt_1 <- rptGaussian(TMaleMinDistance ~ (1|MaleID), grname = "MaleID", data = subset(md, NestStage==1))
mod_lt_2 <- rptGaussian(TMaleMinDistance ~ (1|MaleID), grname = "MaleID", data = subset(md, NestStage==2))

mult_yrs <- md %>% 
  group_by(NestStage, MaleID, AggressionYear) %>% 
  count() %>% 
  filter(NestStage == 2)
sum(duplicated(mult_yrs$MaleID))
sum(duplicated(mult_yrs$MaleID))

tribble(
  ~behaviour,                ~t1,             ~t2,     ~sample_size,                  ~R,   ~CI_lower,                 ~CI_upper,                         ~p_val,
  "aggressiveness",   age_df$age[1], age_df$age[1] + 7,        1027,   unlist(mod_st_1$R),  unlist(mod_st_1$CI_emp)[[1]], unlist(mod_st_1$CI_emp)[[2]], unlist(mod_st_1$P)[[1]],
  "aggressiveness",   age_df$age[1]+8, age_df$age[1] + 15,     1027,   unlist(mod_st_2$R),  unlist(mod_st_2$CI_emp)[[1]], unlist(mod_st_2$CI_emp)[[2]], unlist(mod_st_2$P)[[1]],
  "aggressiveness",   age_df$age[1], age_df$age[1] + 365,      246,    unlist(mod_lt_1$R),  unlist(mod_lt_1$CI_emp)[[1]], unlist(mod_lt_1$CI_emp)[[2]], unlist(mod_lt_1$P)[[1]],
  "aggressiveness",   age_df$age[1]+8, age_df$age[1] + 373,    298,    unlist(mod_lt_1$R),  unlist(mod_lt_1$CI_emp)[[1]], unlist(mod_lt_1$CI_emp)[[2]], unlist(mod_lt_1$P)[[1]],
) %>% 
mutate(Key = "2UJ6RD3M",
       species_common = "great_tit",
       species_latin = "parus_major",
       measurements_per_ind = c(2,2,4,4),
       sex = 0,
       context = 3,
       type_of_treatment = 0,
       treatment = NA,
       life_stage = "adult",
       event = NA,
       R_se = NA,
       delta_t = t2-t1,
       remarks = NA,
       max_lifespan_days = tit_max) -> meta_table

write_delim(meta_table, path = "output/ArayaAjoy_Dingemanse_2017.txt", delim = " ", col_names = TRUE)

###### Paper 75: Freasneau_Kluen_2014 ##############
# Fresneau, Nolwenn; Kluen, Edward; Brommer, Jon E. 2014
# 2CG7EZYD	
# a sex-specific behavioral syndrome in a wild passerine

#  blue tits Cyanistes caeruleus
# 
tit_dat <- scrape_AnAge("Cyanistes caeruleus", vars = "maximum_longevity_yrs")
tit_max <- as.numeric(tit_dat$maximum_longevity_yrs) * 365

avg_age <- tit_max * 0.25
years_avg <- (((73*2 + 41*3 + 10 * 4) / (73 + 41 + 10)) - 1) * 365

tribble(
  ~behaviour,                ~t1,             ~t2,     ~sample_size,  ~ measurements_per_ind,      ~R,   ~R_se,   ~p_val,
 "hatching_defense",       avg_age,     avg_age + 3,     73,         278/73,           0.46,   0.07,   0.0001,
 "hatching_defense",       avg_age,     avg_age + 3,     91,         218 /91,            0.31,   0.08,   0.0001,
 "hatching_defense",       avg_age,     avg_age + 3,     116,        275 /116,            0.39,   0.07,   0.0001,
 "hatching_defense",       avg_age,     avg_age + 3,     84,         204 /84,            0.25,   0.09,   0.0054,
 "hatching_defense",       avg_age,     avg_age + 3,     105,        276 /105,            0.33,   0.07,   0.0001,
 "hatching_defense",       avg_age,     avg_age + 3,     87,         199 /87,            0.22,   0.09,   0.0029,
 "hatching_defense",       avg_age,     avg_age + years_avg, 124,    years_avg/365 + 1,        0.078,  0.066,  0.23
) %>% 
  mutate(Key = "2CG7EZYD",
         species_common = "blue_tit",
         species_latin = "Cyanistes_caeruleus",
         sex = 1,
         context = 3,
         type_of_treatment = 0,
         treatment = NA,
         life_stage = "adult",
         event = NA,
         CI_upper = NA,
         CI_lower = NA,
         delta_t = t2-t1,
         remarks = NA,
         max_lifespan_days = tit_max) -> meta_table

write_delim(meta_table, path = "output/ Freasneau_Kluen_2014.txt", delim = " ", col_names = TRUE)


