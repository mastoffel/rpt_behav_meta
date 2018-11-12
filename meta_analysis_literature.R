# potential tools
# from https://docs.google.com/spreadsheets/d/1n_CmT1kvFmaHUJ3YIzSMUbiREAfIiXCKERub0HZz4vg/edit#gid=0
# library(metagear)# variety of tools
# library(tabulizer) # extract tables from pdfs
# library(pdftools) # keyword search in pdf files
# library(revtools) # github version: literature and screening
# library(digitizeR) # extract data from plots and images


##### Searching Web of Science ##### ---------------------------------------------------------------
# Wildcards:
# * Zero of more characters
# $ Zero or one character
# ? One character. Can be repeated.

# Boolean Operators (capitalization doesnt matter)
# AND  e.g. "stem cell" AND lymphoma
# OR 
# NOT 

# Proximity Operators:
# Phrase searching: "stem cell", wildcard characters are allowed to assess variation in spelling
# NEAR/x: Terms occur within a number of words. Example: canine NEAR/10 virus 

# Search types (Field next to where search terms are entered in Web of Science)
# Topic: searches article titles, abstracts and keywords

# Parentheses:
# group Boolean operators, e.g.
# (personality OR temperament) AND (repeatability or intra$class$coeff*)

# Hyphens
# speech-impairment finds records containing speech-impairment and speech impairment


###### Web of Science search ######
### Try 1 with the following search terms:
# Web of Science Advanced Search, Core Collection, All years 1900-2018, 31.October 2018,10:44 CEST, 781 results
# TS = ((repeatab* OR intraclass) AND behav* AND (personalit* OR temperament* OR (behav* NEAR/10 syndrom*) OR (*individ* NEAR/10 difference*) OR "coping style*")) 

library(revtools) 
library(dplyr)
library(readr)
library(tidyr)
library(wordcloud)

WoS1 <- as.data.frame(read_bibliography("data/WoS_try_1.bib")) %>% 
                dplyr::select(label, title, author, journal, #issn, 
                              volume, pages, year, doi, abstract)
WoS2 <- as.data.frame(read_bibliography("data/WoS_try_2.bib")) %>% 
                dplyr::select(label, title, author, journal, #issn, 
                            volume, pages, year, doi, abstract)

WoS <- rbind(WoS1, WoS2) %>% 
        mutate(source = "WoS")
names(WoS)
table(WoS$year)

# Scopus advanced search, 6.11.2018, 14:00, 771 results
# TITLE-ABS-KEY((repeatab* OR intraclass) AND behav* AND (personalit* OR temperament* OR (behav* W/10 syndrom*) OR (*individ* W/10 difference*) OR "coping style*")) 

scopus <-  as.data.frame(read_bibliography("data/scopus/scopus.bib")) %>% 
    dplyr::select(label, title, author, journal, #issn, 
        volume, pages, year, doi, abstract) %>% 
    mutate(source = "scopus")

all_papers <- rbind(WoS, scopus)

# are any duplicates in the database
any_duplicated <- find_duplicates(all_papers) 

# how many duplicates?
nrow(any_duplicated) - length(unique(any_duplicated$duplicate_group))

# extract unique references
unique_refs <- extract_unique_references(any_duplicated, show_source = TRUE)

# how many were new in scopus?
new_in_scopus <- unique_refs %>% 
    filter(source == "scopus" & n_duplicates == 1)

# compare word frequencies in titles
# wordcloud(paste(new_in_scopus$title, collapse = " "), colors=brewer.pal(6,"Dark2"),random.order=FALSE)
# wordcloud(paste(unique_refs$title, collapse = " "), colors=brewer.pal(6,"Dark2"),random.order=FALSE)
# 
# compare word frequencies in journals
# wordcloud(paste(unique_refs$abstract, collapse = " "), colors=brewer.pal(6,"Dark2"),random.order=FALSE)

# some duplicates still seem to be there, potentially problems with different publication years
# in WoS and scopus ------------------------------------------

# figure out duplicates only by title and remove too 

# all to lowercase for better comparison
unique_refs_lower <- unique_refs %>% 
    mutate(title = tolower(title))

# this will be the vector containing the groups, some of which can be duplicates
dup_groups <- as.numeric(rep(NA, nrow(unique_refs_lower)))

find_duplicates2 <- function(i, unique_refs_lower, dup_group){
    all_stringsims <- unique_refs_lower %>% 
        .[i:nrow(.), ] %>% 
        apply(1, function(x) stringsim(a = unique_refs_lower$title[i], b = x[["title"]]))
    dups <- i-1 + which(as.numeric(all_stringsims) > 0.90, arr.ind = TRUE)
    
    if (!any(is.na(dup_groups[dups]))) return()
    dup_groups[dups] <<- i
}

# fills up dup_groups vector
lapply(1:length(dup_groups), find_duplicates2, unique_refs_lower, dup_group)

group_num_dup <- as.numeric(names(table(dup_groups)[which(as.logical(table(dup_groups) > 1), arr.ind = TRUE)]))

new_unique_refs <- unique_refs_lower %>% 
    mutate(dup_groups) %>% 
    arrange(dup_groups) %>% 
    filter(dup_groups %in% group_num_dup)
new_unique_refs$title

unique_refs_lower <- mutate(unique_refs_lower, duplicate_group = dup_groups)

final_refs <- extract_unique_references(unique_refs_lower) # 30 more removed

# save as csv for rayyan
# rayyan example csv
rayyan <- read_csv("data/rayyan_example.csv")

# check differences in formatting for rayyan input
setdiff(names(rayyan), names(final_refs))
setdiff(names(final_refs), names(rayyan))

# add and rename columns
final_refs_upd <- final_refs %>% 
    dplyr::rename(key = label,
                  authors = author,
                  url = doi) %>% 
    mutate(issue = NA,
           publisher = NA,
           issn = NA)

# check differences again
setdiff(names(rayyan), names(final_refs_upd))
setdiff(names(final_refs_upd), names(rayyan))

# reorder according to ryan
final_refs_rayyan <- final_refs_upd[match(names(rayyan), names(final_refs_upd ))]

# save as csv
write_csv(final_refs_rayyan, path = "data/complete_search_rayyan_formatted.csv")








# 
first_screen <- read_csv("data/first_screening.csv")
names(first_screen)
head(first_screen)
first_screen$notes


WoS_screen <- readRDS("WoS_screen.rds")
revtools::screen_abstracts(WoS)
WoS_screen[WoS_screen$order_random == 1, ]
# load a bibliographic dataset with the authors, titles, and abstracts of multiple study refer
data(example_references_metagear)
names(example_references_metagear)
example_references_metagear["JOURNAL"]

effort_distribute(example_references_metagear, initialize = TRUE, reviewers = "martin", save_split = TRUE)
abstract_screener("effort_martin.csv", aReviewer = "martin")

theBiblio <- scrape_bibliography(as.character(example_references_metagear$DOI)[4])


install.packages("tabulizer")
