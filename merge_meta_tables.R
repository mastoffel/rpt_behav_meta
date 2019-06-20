# merge all tables saved in output
library(tidyverse)
source("martin.R")
library("wordcloud") # word-cloud generator 

length(list.files("output/"))
# read all meta_tables
all_tables <- map(paste0("output/", list.files("output/")), read_delim, delim = " ")
?bind_rows

# convert everything to character
all_tables_char <- map(all_tables, ~mutate_if(., is.numeric, as.character))
# give each table a unique study ID
all_tables_ided <- map(1:length(all_tables_char), function(x) {
                                          current_table <- all_tables_char[[x]]
                                          current_table$study_id <- x
                                          current_table
                                          })
all_tables_ided[[63]]
# some names in some tables are wrong so correct
# names(all_tables_ided[[61]])[names(all_tables_ided[[61]]) == "timeperiod"] <- "life_stage"
#all_tables_ided[[23]]$p_val <- NA
#all_tables_ided[[23]]$session <- NULL

correct_meta_table_names <- function(meta_table){
    if (any(names(meta_table) == "teatment")){
       names(meta_table)[which(names(meta_table) == "teatment")] <- "treatment"
    }
    if (any(names(meta_table) == "lifes_stage")){
        names(meta_table)[which(names(meta_table) == "lifes_stage")] <- "life_stage"
    }
    
    if (any(names(meta_table) == "behavior")){
        names(meta_table)[which(names(meta_table) == "behavior")] <- "behaviour"
    }
    
    if (any(names(meta_table) == "pval")){
        names(meta_table)[which(names(meta_table) == "pval")] <- "p_val"
    }
    
    meta_table
}

meta_tables_corrected <- map(all_tables_ided, correct_meta_table_names)

# create meta_table for all
meta_table_full <- bind_rows(meta_tables_corrected)

# some R values still have ** for significances, let's remove those
meta_table_full$R
as.numeric(meta_table_full$R)

meta_table_corrected2 <- meta_table_full %>% 
    mutate(R = str_replace_all(R, "\\*", ""),
           R = str_replace_all(R, "\\(", ""),
           R = str_replace_all(R, "\\)", "")) 

# one minus seems to be weird
meta_table_corrected2$R[145] <- "-0.083"

# convert a few columns to numeric and create standardised delta
# for blue and great tit from study 35, max lifespan hasn'T been transformed to days
mt_all <- meta_table_corrected2 %>% 
    mutate_at(vars(sample_size, R, R_se, CI_lower, CI_upper, t1, t2, delta_t,
                   measurements_per_ind, max_lifespan_days), as.numeric) %>% 
    select(study_id, Key, species_common, species_latin, behaviour, R, R_se, CI_lower, CI_upper, p_val,
           t1, t2, delta_t, everything()) %>% 
    mutate(delta_stand_by_lifespan = delta_t/max_lifespan_days) 
 

test <- mt_all %>% filter(study_id == 29) %>% as.matrix() %>% as_tibble()
# fishers r to z transformation
# mt_all <- mt_all %>% 
#     mutate(R_stand = (0.5*log(1 +(measurements_per_ind - 1)*R/(1-R))))

WriteXLS::WriteXLS(mt_all, "meta_table_prelim.xls") 

source("martin.R")
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

wordcloud(words = mt_wc$species_common, freq = mt_wc$n, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))


