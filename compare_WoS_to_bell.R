# check whether web of science search terms retrieve the literature from the Bell paper

library(revtools)
library(bib2df)

# used word supplementary from bell paper and parsed it with 'anystyle'. Saved as bib.
bell <- read_bibliography("data/bell_paper.bib")

# make a df for later analysis
bell_df <- bib2df("data/bell_paper.bib")

# Web Of Science search, core collection, All years 1900-2018
# TS = ((repeatab* OR intraclass) AND behav* AND (personalit* OR temperament* OR (behav* NEAR/2 syndrom*) OR (*individ* NEAR/2 difference*))) 740

WoS1 <- read_bibliography("data/WoS1.bib")
WoS2 <- read_bibliography("data/WoS2.bib")[, -ncol(WoS2)]
WoS <- rbind(WoS1, WoS2)

# all names to lowercase for comparability
names(bell_df) <- tolower(names(bell_df))
# figure out which fields are in both bibtex files
congruent_fields <- names(bell_df)[which(names(bell_df) %in% names(WoS1))]

# filter for congruent fields
WoS_short <- WoS[congruent_fields]
bell_short <- bell_df[congruent_fields]
all_lit <- rbind(WoS_short, bell_short)

?find_duplicates

dups <- find_duplicates(all_lit)

ref_unique <- extract_unique_references(all_lit, dups)

screen_duplicates(all_lit)
