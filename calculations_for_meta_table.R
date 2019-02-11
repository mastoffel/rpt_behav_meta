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
average_time_diff_winter1 <- as.numeric(winter1_end - winter1_start) / 11 # see paper, median of 11 measurements

winter2_start <- dmy(01122012)
winter2_end <- dmy(03032013)
average_time_diff_winter2 <- as.numeric(winter2_end - winter2_start) / 10

winter3_start <- dmy(30112013)
winter3_end <- dmy(02032014)
average_time_diff_winter3 <- as.numeric(winter3_end - winter3_start) / 12

# repeatability across all 3 winters (~ 365*2 = 730 days in between)
meta_table <- data.frame("Key" = NA,                # identifier in meta-table "data/meta_table_filled.xlsx"
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


