# explore accepted literature using some wordclouds

library(tidyverse)
library(stringr)
library(wordcloud)

ray_data <- read_csv("data/rayyan_annotated_literature.csv")
names(ray_data)
notes <- ray_data$notes

inclusion_labels <- str_split(notes, "RAYYAN-LABELS: ") %>% 
        lapply(function(x) x[2]) %>% 
        unlist() %>% 
        str_replace(" \\| RAYYAN-EXCLUSION-REASONS: ", ",") %>% # only twice or so in the data
        str_split(",") %>% 
        unlist() %>% 
        .[!is.na(.)] %>% 
        .[!(. == "longer period")] %>% 
        wordcloud(colors=brewer.pal(15,"Dark2"), scale = c(5, 0.8), random.order=FALSE, min.freq = 1)
inclusion_labels

exclusion_labels <- str_split(notes, "RAYYAN-EXCLUSION-REASONS: ") %>% 
    lapply(function(x) x[2]) %>% 
    unlist() %>% 
    str_split(",") %>% 
    unlist() %>% 
    .[!is.na(.)] %>% 
    tibble(reasons = .) %>% 
    ggplot(aes(x = reasons)) + 
    geom_histogram( stat="count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust=0.5))
exclusion_labels 


# wordcloud of longer period studies
abstract_wordcloud_included <- ray_data %>% 
    filter(str_detect(notes, "true")) %>% 
    filter(str_detect(notes, "longer period"))

wordcloud(paste(abstract_wordcloud_included$abstract, collapse = " "),  
    colors=brewer.pal(6,"Dark2"),random.order=FALSE)
wordcloud(abstract_wordcloud_included$abstract,   
    colors=brewer.pal(6,"Dark2"),random.order=FALSE)

included <- notes[str_detect(string = notes, "true")]

included_lables <- str_locate(included, "RAYYAN-LABELS:")
included_lables <- substring(included, c(included:40))
