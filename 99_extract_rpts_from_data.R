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
# variables from the meta table
meta_table_template <- tibble("Key" = NA,                # identifier in meta-table "data/meta_table_filled.xlsx"
                         "species_common" = NA, 
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

##### Paper 13: David_Auclair_2012 ######
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
