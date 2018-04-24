library(tidyverse)
getwd()
setwd("~/github/AS_align/")

chr1a <-as.tibble(c('A', 'A'))
chr3a <-as.tibble(c('A', 'A', 'A', 'AA', 'A', 'A', 'A', 'A'))
chr5a <-as.tibble(c('A', 'A', 'A', 'A'))
chr6a <-as.tibble(c('A', 'A', 'AA', 'A', 'A', 'A', 'A', 'A'))
chr7a <-as.tibble(c('A'))
chr10a <-as.tibble(c('A', 'A'))
chr12a <-as.tibble(c('A', 'A'))
chr16a <-as.tibble(c('A', 'A'))

(SF1a<-chr1a%>%bind_rows(chr3a, chr5a, chr6a, chr7a, chr10a, chr12a, chr16a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(SF="1"))

chr1b <-as.tibble(c('B', 'B'))
chr3b <-as.tibble(c('B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'))
chr5b <-as.tibble(c('B', 'B', 'B', 'B'))
chr6b <-as.tibble(c('BBB', 'B', 'B', 'B', 'B', 'B', 'B', 'B'))
chr7b <-as.tibble(c('B'))
chr10b <-as.tibble(c('BB', 'B'))
chr12b <-as.tibble(c('B', 'B'))
chr16b <-as.tibble(c('B', 'B'))

(SF1b<-chr1b%>%bind_rows(chr3b, chr5b, chr6b, chr7b, chr10b, chr12b, chr16b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(SF="1")) 

chr2a <-as.tibble(c('A', 'A', 'A', 'A'))
chr4a <-as.tibble(c('A', 'A', 'A'))
chr8a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A', 'A'))
chr9a <-as.tibble(c('A'))
chr13a <-as.tibble(c('AA', 'A', 'A', 'A', 'A'))
chr15a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A'))
chr18a <-as.tibble(c('A', 'A', 'A', 'A', 'A'))
chr20a <-as.tibble(c('A', 'A'))

(SF2a<-chr2a%>%bind_rows(chr4a, chr8a, chr9a, chr13a, chr15a, chr18a, chr20a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(SF="2")) 

chr2b <-as.tibble(c('B', 'B', 'B', 'B'))
chr4b <-as.tibble(c('B', 'BB', 'B'))
chr8b <-as.tibble(c('BB', 'B', 'B', 'B', 'B', 'B', 'B'))
chr9b <-as.tibble(c('B'))
chr13b <-as.tibble(c('B', 'B', 'B', 'B', 'B'))
chr15b <-as.tibble(c('B', 'BB', 'B', 'B', 'BB', 'B'))
chr18b <-as.tibble(c('BB', 'B', 'B', 'B', 'B'))
chr20b <-as.tibble(c('BB', 'B'))

(SF2b<-chr2b%>%bind_rows(chr4b, chr8b, chr9b, chr13b, chr15b, chr18b, chr20b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(SF="2")) 

chr11a <-as.tibble(c('AA'))
chr17a <-as.tibble(c('AA', 'A', 'AA', 'A'))
chr17ba <-as.tibble(c('AAA', 'AA'))
chrXa <-as.tibble(c('AA', 'AA', 'A'))

(SF3a<-chr11a%>%bind_rows(chr17a, chr17ba, chrXa)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(SF="3")) 

chr11b <-as.tibble(c('BBB'))
chr17b <-as.tibble(c('BBB', 'BBB', 'B', 'BBB'))
chr17bb <-as.tibble(c('BBBBB', 'B'))
chrXb <-as.tibble(c('B', 'BBB', 'BBB'))

(SF3b<-chr11b%>%bind_rows(chr17b, chr17bb, chrXb)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(SF="3")) 

(E<-as.tibble(c(0:5))%>%mutate(A=c("A","AA","AAA","AAAA","AAAAA","AAAAAA"))%>%mutate(B=c("B","BB","BBB","BBBB","BBBBB","BBBBBB"))%>%mutate(dg=dgeom(value,0.5))%>%mutate(qua="E"))

(E1<-E%>%mutate(cnt=dg*29)%>%mutate(SF="1"))
E2<-E%>%mutate(cnt=dg*33)%>%mutate(SF="2")
E3<-E%>%mutate(cnt=dg*10)%>%mutate(SF="3")
Ea<-E1%>%bind_rows(E2,E3)%>%select(A,cnt,SF,qua)
Ea

Eb<-E1%>%bind_rows(E2,E3)%>%select(B,cnt,SF,qua)
Eb

SFa<-SF1a%>%bind_rows(SF2a,SF3a)%>%mutate(qua="O")%>%rename(A=value)
SFa<-SFa%>%bind_rows(Ea)%>%bind_rows(tibble(A=c("AAA","AAAA","AAAAA","AAAAAA"),cnt=c(0,0,0,0),SF=c("1","1","1","1"),qua=c("O","O","O","O")))%>%bind_rows(tibble(A=c("AAA","AAAA","AAAAA","AAAAAA"),cnt=c(0,0,0,0),SF=c("2","2","2","2"),qua=c("O","O","O","O")))%>%bind_rows(tibble(A=c("AAAA","AAAAA","AAAAAA"),cnt=c(0,0,0),SF=c("3","3","3"),qua=c("O","O","O")))
SFa
g_a<-SFa%>%ggplot(aes(x=A, y=cnt, group=qua, position=qua, fill=qua))
g_a<-g_a+geom_col(position = "dodge")+facet_grid(SF~.,margins = TRUE)
print(g_a)

SFb<-SF1b%>%bind_rows(SF2b,SF3b)%>%mutate(qua="O")%>%rename(B=value)
SFb<-SFb%>%bind_rows(Eb)%>%bind_rows(tibble(B=c("BBBB","BBBBB","BBBBBB"),cnt=c(0,0,0),SF=c("1","1","1"),qua=c("O","O","O")))%>%bind_rows(tibble(B=c("BBB","BBBB","BBBBB","BBBBBB"),cnt=c(0,0,0,0),SF=c("2","2","2","2"),qua=c("O","O","O","O")))%>%bind_rows(tibble(B=c("BB","BBBB","BBBBBB"),cnt=c(0,0,0),SF=c("3","3","3"),qua=c("O","O","O")))
g_b<-SFb%>%ggplot(aes(x=B, y=cnt, group=qua, position=qua, fill=qua))
g_b<-g_b+geom_col(position = "dodge")+facet_grid(SF~.,margins = TRUE)
print(g_b)






