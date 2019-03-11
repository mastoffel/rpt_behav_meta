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
# variables from the meta table
meta_table_template <- tibble("Key" = NA,                # identifier in meta-table "data/meta_table_filled.xlsx"
                         "species_common" = NA, 
                         "species_latin" = NA,
                         "sample_size" = NA, 
                         "measurements_per_ind" = NA, # new, check papers 1-6 again
                         "sex" = NA,                # 0 = both, 1 = females, 2 = males
                         "behaviour" = NA,          # measured behaviour as stated by authors
                         "context" = NA,            # 1 = lab exp. / lab-reared, 2 = lab exp. / wild-caught, 3 field exp / maybe another category: 4 field behaviour?
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
                         "remarks"= NA,
                         "max_lifespan_days" = NA)


###### Paper 1: Baker_Goodman 2018 ##########
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


###### Paper 3, Bierbach_Laskowski 2017 incomplete #######
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
write_delim(meta_table, path = "output/Aplin_Firth_2015.txt", delim = " ", col_names = TRUE)


###### Paper 5, Arroyo_Mougeot 2017 incomplete #######
# Arroyo, Beatriz; Mougeot, Francois; Bretagnolle, Vincent
# individual variation in behavioural responsiveness to humans leads to differences in breeding success and long-term population phenotypic changes
# Key: V5R8XCKE
# DOI: 10.1111/ele.12729

# a bit tricky: Sample sizes are a bit unclear (number of measurements are reported only, see table)
# between years consists of 30 individuals with 2-7 years 

meta_table <- tribble(
    ~behaviour,                   ~R,  ~R_se, ~sample_size, ~t1, ~t2, ~delta_t, ~remarks, ~sex,
    "nest_departure_distance", 0.001,  0.023,         1662,  NA,  NA,      365, 'within-year-repeatability', 1,
    "fleeing probability",     0.541,  0.049,         2399,  NA,  NA,      365, 'within-year-repeatability', 1,    
    "passive_if_present",      0.520,  0.059,         1844,  NA,  NA,      365, 'within-year-repeatability', 1, 
    "defence_intensity_PCA",   0.518,  0.033,         1014,  NA,  NA,      365, 'within-year-repeatability', 1,        
    "nest_departure_distance", 0.001,  0.017,          197,  NA,  NA,      365, 'between-year-repeatability (range 2-7 years)', 1,
    "fleeing probability",     0.341,  0.135,          299,  NA,  NA,      365, 'between-year-repeatability (range 2-7 years)', 1,    
    "passive_if_present",      0.522,  0.194,          252,  NA,  NA,      365, 'between-year-repeatability (range 2-7 years)', 1, 
    "defence_intensity_PCA",   0.131,  0.073,          189,  NA,  NA,      365, 'between-year-repeatability (range 2-7 years)', 1
)

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
    
meta_table <- meta_table %>% 
    mutate(Key = "RR968ED6",
           species_common = "largemouth_bass",
           context = 2,
           type_of_treatment = 0,
           treatment = NA,
           p_val = NA,
           event = NA)

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
            
write_delim(meta_table, path = "output/Bosco_Riechert_2016.txt", delim = " ", col_names = TRUE)



###### Paper 8 Boulton et al 2014 To be discussed ####
# Boulton, Kay; Grimmer, Andrew J.; Rosenthal, Gil G.; Walling, Craig A.; Wilson, Alastair J. 2014		
# PYZL2E8A	
# how stable are personalities? a multivariate view of behavioural variation over long and short 
#timescales in the sheepshead swordtail, xiphophorus birchmanni


### talk to holger, not sure this is ok, as individuals for short and long timescales are different! ##





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


meta_table <- tribble(
    ~behaviour,              ~R,  ~R_se, ~CI_lower, ~CI_upper, ~sample_size,  ~t1,   ~t2, ~delta_t, ~measurements_per_ind, ~remarks,
    "boldness",            0.61,     NA,      0.57,      0.66,          458, 6753,   8760,    2008,                    5.5, "one measurement per year and average of 5.5 equals 5.5 years span of rpt",
    "boldness",            0.82,     NA,      0.54,      0.92,           35, 6753,  7118,      365,                      3, "measurement times start based on mean adult life time" 
)

meta_table <- meta_table %>% 
    mutate(Key = "AKG6LAL3",
           species_common = "grey_seal",
           sex = 1,
           context = 3,
           type_of_treatment = NA,
           treatment = NA,
           life_stage = "adult",
           event = NA,
           p_val = NA)

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


write_delim(meta_table, path = "output/Burtka_Grindstaff_2013.txt", delim = " ", col_names = TRUE)




###### Paper 11: Cabrera_Andres_2017 #######
 #Cabrera, Doreen; Andres, Daniel; McLoughlin, Philip D.; Debeffe, Lucie; Medill, Sarah A.; Wilson, Alastair J.; Poissant, Jocelyn	2017	
# ESF6PD88	
# island tameness and the repeatability of flight initiation distance in a large herbivore

location <- "data/papers/Cabrera et al. - 2017 - island tameness and the repeatability of flight in-rotated.pdf"

# average time of measurements among days: 
delta_t_among <- 10.2
# all data from foles (maximum of one year old, but rather a few days to a few months)
# assumed here: 365/2 = 182.5 days
age_fole <- 365/2

tribble(
    ~timespan,         ~sample_size, ~measurements_per_ind,    ~R, ~R_se,     ~t1,                           ~t2,      ~delta_t,
  #  "within_among_days",          103,               376/103,   0.42,  0.06,  age_fole, age_fole + delta_t_among, delta_t_among,
    "within_days",                156,               375/156,   0.55,  0.05,  age_fole,           age_fole + 0.5,           0.5,
    "among_days",                  45,                 99/45,   0.39,  0.12,  age_fole, age_fole + delta_t_among, delta_t_among
) %>% 
    mutate(Key = "ESF6PD88",
           species_common = "horse",
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
           remarks = "foles under 1 year") %>% 
    select(-timespan) -> meta_table

write_delim(meta_table, path = "output/Cabrera_Andres_2017.txt", delim = " ", col_names = TRUE)


###### Paper 12: D'Eath_2004 #######
# D'Eath, RB 2004	
# 2WGN2QQC		
# consistency of aggressive temperament in domestic pigs: the effects of social experience and social disruption

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
           sex = 0,
           behaviour = "aggressiveness",
           context = 1,
           type_of_treatment = c(0,0,0,2),
           teatment = c(NA, NA, NA, "litter_change"),
           life_stage = c("juvenile", "adult", "adult", "both"),
           event = NA,
           CI_lower = NA,
           CI_upper = NA,
           remarks = c("correlation", "correlation", "correlation", "repeatability")
           ) -> meta_table

write_delim(meta_table, path = "output/D'eath_2004.txt", delim = " ", col_names = TRUE)

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
           sex = 1,
           context = 1,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = NA,
           R_se = NA,
           remarks = NA
    ) -> meta_table
    
write_delim(meta_table, path = "output/David_Auclair_2012.txt", delim = " ", col_names = TRUE)

###### Paper 14: DeWitt_Sih_1999 ####
# DeWitt, TJ; Sih, A; Hucko, JA	1999
# 7QAACTY6	
# trait compensation and cospecialization in a freshwater snail: size, shape and antipredator behaviour

# snail were wild caught

# To quantify the repeatability of individual variation in antipredator behaviour, we repeated our
# behavioural assays on 6 days spread over a 13-day period.
tribble(
    ~behaviour,   ~sample_size, ~measurements_per_ind,    ~R, ~CI_lower, ~CI_upper,  ~p_val,    ~t1,   ~t2, ~delta_t,
    "antipredator_behaviour",         96,           3,   0.33,      NA,        NA,   0.0001,     NA,    NA,         3,  
    "antipredator_behaviour",         96,           6,   0.27,      NA,        NA,   0.0001,     NA,    NA,        13     
) %>% 
    mutate(Key = "7QAACTY6",
           species_common = "physid_snail",
           context = 2,
           remarks = "No SE/CI, also data was obtained by pooling across video observation, i.e. measurement per ind not reliable",
           sex = 0,
           context = 2, 
           type_of_treatment = 0,
           treatment = NA, 
           lifes_stage = NA,
           event = NA,
           R_se = NA) -> meta_table
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

delta_t_between <- (1.7*365)+185

tribble(
    ~behaviour,   ~sample_size, ~measurements_per_ind,    ~R, ~CI_lower, ~CI_upper,  ~p_val,    ~t1,                    ~t2,        ~delta_t, ~remarks,
    "docility",         116,           4,               0.41,      0.30,      0.51,      NA,   1132,              1132+185,             185, "delta_t is full season",
    "docility",          65,           2.7,             0.31,      0.13,      0.49,      NA,   1132,  1132+delta_t_between, delta_t_between, "delta_t is full season plus an 2.7*365 days"   
) %>% 
    mutate(Key = "4ABCSFR6",
           species_common = "roe_deer",
           sex = 0,
           context = 3,
           type_of_treatment = 0,
           life_stage = "both",
           event = NA,
           treatment = NA,
           R_se = NA) -> meta_table

write_delim(meta_table, path = "output/Debeffe_Lemaitre_2015.txt", delim = " ", col_names = TRUE)








###### Paper 16: English_Nakagawa_2010 incomplete, longitudinal? #######
# English, S.; Nakagawa, S.; Clutton-Brock, T. H. 2010	
# XR32CK8A	
# consistent individual differences in cooperative behaviour in meerkats (suricata suricatta)

# 2-19 measures per individual (average 10 ?)
# few individuals survive beyond 4 years of age
# tribble(
#     ~behaviour,     ~R,   ~R_se, ~p_val, ~N_o, N_i,   ~t1,                    ~t2,        ~delta_t, ~remarks,
#   #  "babysitting",  0.218, 0.046, 0.0001, 6460, 646,                              
# )


###### Paper 17: Erhard_Mendl_1997 #######

# Erhard, HW; Mendl, M		1997	
# SX3CYX5C
# measuring aggressiveness in growing pigs in a resident-intruder situation

# 1 and 2: 11 weeks of age and 1 day later

tribble(
    ~behavior,         ~R, ~R_se, ~p_val, ~sample_size, ~measurements_per_ind, ~t1, ~t2, ~delta_t,
    "aggressiveness",  0.56, NA,     0.01,   85,                              2,  77,  78,        1,
    "aggressiveness",  0.73, NA,     0.01,   78,                              2,  77,  78,        1,
    "aggressiveness",  0.57, NA,     0.01,   53,                              2,  49,  77,        28
) %>% 
    mutate(Key = "SX3CYX5C",
           species_common = "pig",
           sex = 0,
           context = 1,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "juvenile",
           event = NA,
           CI_lower = NA,
           CI_upper = NA,
           remarks = "All spearman rank correlations, no SE, not the same individuals for short and long term rpt") -> meta_table

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



tribble(
      ~R, ~CI_lower, ~CI_upper, ~sample_size, ~measurements_per_ind,        ~t1,       ~t2,   ~delta_t,  ~p_val, 
    0.49,      0.33,      0.64,          44,                   3.4,        1825,    1825+365,      365,  0.0001,
    0.41,      0.28,      0.54,          65,                   3.4,  1825+2*365,  1825+3*365,      365,  0.0001,
    0.74,      0.49,      0.88,          20,                     2,        1825,  1825+3*365,    3*365,  0.0001   
) %>% 
    mutate(Key = "NWYGHZWI",
           species_common = "stellers_jay",
           behaviour = "risk_taking",
           sex = 0,
           context = 3,
           type_of_treatment = 0,
           treatment = NA,
           life_stage = "adult",
           event = "between_season",
           R_se = NA,
           remarks = "sample sizes and measurement were a bit of guesswork") -> meta_table

write_delim(meta_table, path = "output/Gabriel_Black_2010.txt", delim = " ", col_names = TRUE)



###### Paper 20: Garamszegi_Mark_2015 ######

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

write_delim(meta_table, path = "output/Garamszegi_Mark_2015.txt", delim = " ", col_names = TRUE)




###### Paper 21: Gifford_Clay_2014 : FROM here also species_latin variable and avg_adult var when age not reported ######

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
           remarks = NA) -> meta_table

write_delim(meta_table, path = "output/Gifford_Clay_2014.txt", delim = " ", col_names = TRUE)









###### Paper 22: Goold_Newberry_2017 incomplete, data_to_be_analysed #######

url <- "https://github.com/ConorGoold/GooldNewberry_aggression_shelter_dogs/blob/master/raw_data.csv"

dat <- read_csv(url)

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
todigit <- metaDigitise("to_digitise/study3_grace2014/")
 # write_delim(todigit, path = "to_digitise/study3_grace2014/digitized_figure.txt")
todigit <- read_delim("to_digitise/study3_grace2014/digitized_figure.txt", delim = " ")
# to be figured out
avg_adult <- NA


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
           t1 = NA,
           t2 = NA,
           remarks = "age_unknown",
           Key = "X75NVMBP",
           species_common = "nazca_booby",
           species_latin = "sula_granti") -> meta_table

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
           life_stage = "adult",
           R_se = NA,
           p_val = NA) -> meta_table

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
        remarks = "no se or ci or pval, maybe simulation?") -> meta_table
    
write_delim(meta_table, path = "output/Haage_Bergvall_2013.txt", delim = " ", col_names = TRUE)


###### Paper 27: Hammond-Tooke_Cally_2012, metadigitise here #######   
### metadigitise
# Hammond-Tooke, Cally A.; Nakagawa, Shinichi; Poulin, Robert 2012	
# 6TGIGAET	
# parasitism and behavioural syndromes in the fish gobiomorphus cotidianus

# common_bullies
# Gobiomorphus cotidianus
# "Young bullies" between 25 and 55 mm in length

# fish from wild to lab
 # three sessions, three weeks between each
# one session : fish tested twice, 7 days apart
# second session with predator cue added

# all young fish, so guess: half a year old (they mature at 1 year)
avg_juvenile <- 183

# digitalise
# todigit <- metaDigitise("to_digitise/study4_hammond-tooke_cally2012/")
# write_delim(x = todigit, path = "to_digitise/study4_hammond-tooke_cally2012/rpt_table.txt")
todigit <- read_delim(path = "to_digitise/study4_hammond-tooke_cally2012/rpt_table.txt")

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
           CI_lower = NA,
           CI_upper = NA,
           remarks = "age_guessed") -> meta_table
    
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

# common_voles 
# Microtus_arvalis
# 1) before and after maturation over three month 62 +- 20 days old and 90 days later after mat. lab-born , 17 voles 9 males 8 females
# 2) adult during one week / different voles, adult voles (avg life-span 4.5 month), 88 males, 80 females, wild trapped and lab tested
# 3) adult over three month, 48 adults males and females 2.5 month apart wild trapped and lab tested
# 
# assumed average adult age: 0.25 * 4.8 year = 1.2 years = 438 days

avg_adult_age <- 438

meta_table <- meta_table_raw %>% 
    mutate(sample_size = as.numeric(sample_size),
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

# dat <- metaDigitise("to_digitise/study5_Hudson_Rangassamy_2015/", summary = FALSE)
# write_delim(dat, "to_digitise/study5_Hudson_Rangassamy_2015/digitised_dat.txt")
# dat <- read_delim(file = "to_digitise/study5_Hudson_Rangassamy_2015/digitised_dat.txt", delim = " ")
# dat

library(digitize)

mydata <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot1.png")
mydata2 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot2.png")
mydata3 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot3.png")
mydata4 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot4.png")
mydata5 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot5.png")
mydata6 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot6.png")
mydata7 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot7.png")
mydata8 <- digitize("to_digitise/study5_Hudson_Rangassamy_2015/plot8.png")


# saveRDS(ls(), file = "to_digitise/study5_Hudson_Rangassamy_2015/all_digitised_df.RData")
# obj <- readRDS(file = "to_digitise/study5_Hudson_Rangassamy_2015/all_digitised_df.RData")
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
mydata <- digitize("to_digitise/study6/plot1.png")
mydata2 <- digitize("to_digitise/study6/plot2.png")
mydata3 <- digitize("to_digitise/study6/plot3.png")


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
    0.041, 0.128, 0.356, avg_adult_age,  avg_adult_age + 72,          72,            52,
    0    , 0.134, 0.758, avg_adult_age,  avg_adult_age + 365,        365,            34
) %>% 
    mutate(Key = "ZCSXGUYW",
        species_common = "collared_flycatcher",
        species_latin = "Ficedula albicollis",
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
        remarks = "R_se estimated from recalculating rpt with 
        a simplified model (only id as random effect)") -> meta_table

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


    
###### Paper 37: