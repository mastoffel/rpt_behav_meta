# merge all tables saved in output
library(tidyverse)
source("theme_simple.R")

# path to files
files <- list.files("output", full.names = TRUE)

# how many studies?
length(files)

# check that all metatables have the same format
col_names <- files %>% 
    .[1] %>% 
    read_delim(" ") %>% 
    names(.)

col_types <- files %>% 
    .[1] %>% 
    read_delim(" ") %>% 
    map_chr(function(x) class(x))

# read all meta_tables
all_tables <- map(files, read_delim, 
                  delim = " ")

# check whether all files have all columns and are mergeable
check_table <- function(i, col_names) {
    missing_cols <- col_names[!(col_names %in% names(all_tables[[i]]))]
    if (length(missing_cols != 0)) print(paste0(files[i], " lacks column(s) " , missing_cols))
}
# if this doesn't return anything, things are fine. 
walk(1:length(all_tables), check_table, col_names)

# convert everything to character for merging
all_tables_char <- map(all_tables, ~mutate_if(., is.numeric, as.character))
# give each table a unique study ID
all_tables_ided <- map(1:length(all_tables_char), function(x) {
                                          current_table <- all_tables_char[[x]]
                                          current_table$study_id <- x
                                          current_table
                                          })

# create meta_table for all
meta_table_full <- bind_rows(all_tables_ided)

# some R values still have ** for significances, let's remove those, as they
# are converted to NA otherwise when made numeric
meta_table_full$R
meta_table_corrected <- meta_table_full %>% 
    mutate(R = str_replace_all(R, "\\*", ""),
           R = str_replace_all(R, "\\(", ""),
           R = str_replace_all(R, "\\)", "")) 

# check whether R converts properly
which(is.na(as.numeric(meta_table_corrected$R)))

# one minus seems to be weird
meta_table_corrected$R[149] <- "-0.083"

# convert a few columns to numeric and create standardised delta
# for blue and great tit from study 35, max lifespan hasn'T been transformed to days
mt_all <- meta_table_corrected %>% 
    mutate_at(vars(sample_size, R, R_se, CI_lower, CI_upper, t1, t2, delta_t,
                   measurements_per_ind, max_lifespan_days), as.numeric) %>% 
    select(study_id, Key, species_common, species_latin, behaviour, R, R_se, CI_lower, CI_upper, p_val,
           t1, t2, delta_t, everything()) %>% 
    mutate(delta_stand_by_lifespan = delta_t/max_lifespan_days) 

# fill in missing SE's \ CI's where only one of them is given
# note: some CI's become negative here. I didn't want to set them to 0
# because sometimes the estimate is not the repeatability but something
# like spearman rank correlation
mt_all <- mt_all %>% 
    mutate(R_se = ifelse(is.na(R_se), (CI_upper-CI_lower)/3.92, R_se),
           CI_lower = ifelse(is.na(CI_lower), R - 1.96*R_se, CI_lower),
           CI_upper = ifelse(is.na(CI_upper), R + 1.96*R_se, CI_upper))

mt_all %>% 
    filter(is.na(R_se)) 
# fishers r to z transformation
# mt_all <- mt_all %>% 
#     mutate(R_stand = (0.5*log(1 +(measurements_per_ind - 1)*R/(1-R))))

WriteXLS::WriteXLS(mt_all, "data/meta_table_prelim.xls") 


# some plots
source("theme_simple.R")
mt_all %>% 
    mutate(R = ifelse(R<0, 0, R)) %>% 
    ggplot(aes(delta_stand_by_lifespan, R)) +
    geom_point(aes(size = sample_size), alpha = 0.2) +
    geom_smooth(method = "lm") +
    scale_x_sqrt() +
    theme_martin()

mt_all

mt_all %>% 
    group_by(species_common) %>% 
    count() %>% 
    arrange(n) %>% 
    .$species_common -> species_names_ordered

mt_all %>% filter(species_common == "desert_funnel_web_spider") %>% .$study_id


mt_all %>% 
    mutate(species_common = forcats::fct_relevel(species_common ,species_names_ordered)) %>% 
    group_by(study_id, species_common) %>% 
   # count() %>% 
    #arrange(species_common) %>% 
    ggplot(aes(species_common)) +
    geom_bar(aes(fill = as.factor(study_id))) +
    theme_martin() +
    coord_flip() +
    ylab("number of effect sizes") + 
    xlab("species") +
    ggtitle("Number of effect sizes per species.\n different color = different study") +
    theme(legend.position = "none")
    
    
as.data.frame(table(mt_all$species_common)) %>% 
    rename("species" = Var1,
           "n_effect_sizes" = Freq) %>% 
    arrange(n_effect_sizes) %>% 
    ggplot()

library(wordcloud)
wordcloud(words = mt_wc$species_common, freq = mt_wc$n, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))


