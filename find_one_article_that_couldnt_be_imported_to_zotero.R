library(readr)
library(stringr)

library("bib2df")

test <- readLines("data/rayyan_longer_period_articles.bib")

sum(str_detect(test, "@article"))

zotero <- bib2df(file = "data/zotero_import.bib", separate_names = TRUE)

zotero <- read_csv(file = "data/zotero_import.csv")

before_zotero <- bib2df(file = "data/rayyan_longer_period_articles.bib", separate_names = FALSE)

before_zotero_titles <- tolower(before_zotero$TITLE)
zotero_titles <- tolower(zotero$Title)

which(!(before_zotero_titles %in% zotero_titles))

before_zotero_titles[30]
before_zotero_titles[69] # its this one, probably because of a lack of doi or so
