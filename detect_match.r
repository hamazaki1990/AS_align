library(tidyverse)
setwd("~/github/AS_align/")
(match <-read_csv("SF3_5mer_exactmatch20.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match %>% group_by(fragment)%>%mutate(.,"count"=n())%>%distinct()%>%ungroup()%>%filter(count>2)%>%View()

(match <-read_csv("detect_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match%>%View()
match %>% group_by(fragment)%>%mutate(.,"count"=n())%>%distinct()%>%ungroup()%>%filter(count>2)%>%View()

(match <-read_csv("SF1_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match %>% group_by(fragment)%>%mutate(.,"count"=n())%>%distinct()%>%ungroup()%>%filter(count>2)%>%View()

(match <-read_csv("SF2_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match %>% group_by(fragment)%>%mutate(.,"count"=n())%>%distinct()%>%ungroup()%>%filter(count>2)%>%arrange( hit_id, hit_start, fragment)%>%View()

(match <-read_csv("SF3_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match %>% group_by(fragment)%>%mutate(.,"count"=n())%>%distinct()%>%ungroup()%>%filter(count>2)%>%arrange(fragment, hit_id)%>%View()


(match <-read_csv("SF3_2mer_junc_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match %>% group_by(fragment)%>%mutate(.,"count"=n())%>%distinct()%>%ungroup()%>%filter(count>2)%>%arrange(fragment, hit_id)%>%View()

(match <-read_csv("SF1_SF3_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match%>%View()

(match <-read_csv("SF2_SF3_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match%>%View()

(match <-read_csv("SF2_SF1_exactmatch.csv", col_names = c("fragment", "hit_id", "hit_start", "seq")))
match%>%View()
