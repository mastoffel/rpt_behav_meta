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
# variables from the meta table
meta_table_template <- tibble("Key" = NA,                # identifier in meta-table "data/meta_table_filled.xlsx"
                         "species_common" = NA, 
                         "species_latin" = NA,
                         "sample_size" = NA, 
                         "measurements_per_ind" = NA, # new, check papers 1-6 again
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


