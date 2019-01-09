# constructing the meta-table.

# What happened so far:
# (1) WoS and Scopus search (see file 1_)
# (2) Rayyan abstract screening: 218 articles tagged with longer period / more than 2 measurements
# (3) exported to bibtex, changed "url" to "doi" in bibtex file, imported in zotero
# (4) downloaded pdfs in zotero and not exported csv file with standard utf encoding and "notes" 

# What will happen here: file preparation for the meta-table

library(tidyverse)
library(writexl)
# load csv with all articles
zotero_articles <- read_csv("data/rayyan_longer_period_articles_from_zotero_04_12_2018.csv")
names(zotero_articles)

meta_table <- zotero_articles %>% 
                # select relevant fields to identify study / leave the names unchanged for 
                # easier merging later
                select(Key, `Publication Year`, Author, Title, `Publication Title`, ISBN, ISSN,
                    DOI, Pages, Volume) %>% 
                # add variables to be collected for metaanalysis
                mutate(species_common = "", 
                       sample_size = "",
                       sex = "", # 0 = both, 1 = females, 2 = males
                      # species_latin = "",
                       behaviour = "", # broader behaviour as described by the authors (boldness, etc..)
                                       # see Bell 2009 for a potential classification
                       context = "", # 3 levels: 1: lab experiment, lab-reared animals
                                     #           2: lab experiment, wild-caught animals
                                     #           3: Field experiment
                       type_of_treatment = "", # 1: Between-subject treatment
                                               # 2: Within-subject treatment
                                               # NA: No treatment
                       treatment = "", # verbal description of treatment (hormones etc..)
                       life_stage = "", # to be build: potentially: juvenile, adult, both (maybe more 
                                       # specific, larvae or so)
                       event = "", # any big life event between measurements?
                       R = "", # point estimate repeatability
                       R_se = "", # standard error or r
                       CI_lower = "", # lower CI
                       CI_upper = "", # upper CI
                       p_val = "",
                       t1 = "", # timepoint of first measurement in days old#
                       t2 = "", # timepoint of second measurement in days old#
                       delta_t = "", # difference between timepoints 
                       remarks = ""
                    )

write_xlsx(meta_table, path = "data/meta_table.xlsx")
