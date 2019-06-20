# dingemanse paper mit data (1)
# nest stage 1: one model, nest stage 2: another model
# random: brood id, male id, nest box

# innerhalb von Brut: 5 week
# zwischen Bruten innerhalb von Jahr: Mai/Juni 6 * 7 Tage

# kurzfristig
summary(lmer(TMaleMinDistance ~ (1|BroodID), data=subset(md, NestStage==1)))
# langfristig
summary(lmer(TMaleMinDistance ~ (1|MaleID) + (1|BroodID), data=subset(md, NestStage==1)))

library(tidyverse)
dat <- read_csv("data/downloaded_from_papers/dingemanse_great_tit.csv")
table(dat$MaleID)
table(dat$AggressionYear)

dat %>% 
    group_by(MaleID, AggressionYear) %>% 
    count() -> counts_n


# (2) fresneau paper: 

# hatching defense 
# within year vs across year rpt is there

# (3) fisher

url1 <- "https://ore.exeter.ac.uk/repository/bitstream/handle/10871/16930/bae.13.txt?sequence=3&isAllowed=y"
dat <- read_delim(url1, delim = "\t")

