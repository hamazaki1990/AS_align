library(tidyverse)
install.packages("viridis")
library(viridisLite)

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

(SF1a<-chr1a%>%bind_rows(chr3a, chr5a, chr6a, chr7a, chr10a, chr12a, chr16a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF1"))

chr1b <-as.tibble(c('B', 'B'))
chr3b <-as.tibble(c('B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'))
chr5b <-as.tibble(c('B', 'B', 'B', 'B'))
chr6b <-as.tibble(c('BBB', 'B', 'B', 'B', 'B', 'B', 'B', 'B'))
chr7b <-as.tibble(c('B'))
chr10b <-as.tibble(c('BB', 'B'))
chr12b <-as.tibble(c('B', 'B'))
chr16b <-as.tibble(c('B', 'B'))

(SF1b<-chr1b%>%bind_rows(chr3b, chr5b, chr6b, chr7b, chr10b, chr12b, chr16b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF1")) 

chr2a <-as.tibble(c('A', 'A', 'A', 'A'))
chr4a <-as.tibble(c('A', 'A', 'A'))
chr8a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A', 'A'))
chr9a <-as.tibble(c('A'))
chr13a <-as.tibble(c('AA', 'A', 'A', 'A', 'A'))
chr15a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A'))
chr18a <-as.tibble(c('A', 'A', 'A', 'A', 'A'))
chr20a <-as.tibble(c('A', 'A'))

(SF2a<-chr2a%>%bind_rows(chr4a, chr8a, chr9a, chr13a, chr15a, chr18a, chr20a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF2")) 

chr2b <-as.tibble(c('B', 'B', 'B', 'B'))
chr4b <-as.tibble(c('B', 'BB', 'B'))
chr8b <-as.tibble(c('BB', 'B', 'B', 'B', 'B', 'B', 'B'))
chr9b <-as.tibble(c('B'))
chr13b <-as.tibble(c('B', 'B', 'B', 'B', 'B'))
chr15b <-as.tibble(c('B', 'BB', 'B', 'B', 'BB', 'B'))
chr18b <-as.tibble(c('BB', 'B', 'B', 'B', 'B'))
chr20b <-as.tibble(c('BB', 'B'))

(SF2b<-chr2b%>%bind_rows(chr4b, chr8b, chr9b, chr13b, chr15b, chr18b, chr20b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF2")) 

chr11a <-as.tibble(c('AA'))
chr17a <-as.tibble(c('AA', 'A', 'AA', 'A'))
chr17ba <-as.tibble(c('A','AA', 'AA'))
chrXa <-as.tibble(c('AA', 'AA', 'A'))

(SF3a<-chr11a%>%bind_rows(chr17a, chr17ba, chrXa)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF3")) 

chr11b <-as.tibble(c('BBB'))
chr17b <-as.tibble(c('BBB', 'BBB', 'B', 'BBB'))
chr17bb <-as.tibble(c('BBB','BBBBB', 'B'))
chrXb <-as.tibble(c('B', 'BBB', 'BBB'))

(SF3b<-chr11b%>%bind_rows(chr17b, chr17bb, chrXb)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF3")) 

(E1a<-as.tibble(c(0:4))%>%mutate(A=c("A","AA","AAA","AAAA","AAAAA"))%>%mutate(frequency=dgeom(value,31/59))%>%mutate(qua="Expected")%>%mutate(SF="SF1"))
(E1b<-as.tibble(c(0:4))%>%mutate(B=c("B","BB","BBB","BBBB","BBBBB"))%>%mutate(frequency=dgeom(value,28/59))%>%mutate(qua="Expected")%>%mutate(SF="SF1"))
(E2a<-as.tibble(c(0:4))%>%mutate(A=c("A","AA","AAA","AAAA","AAAAA"))%>%mutate(frequency=dgeom(value,34/70))%>%mutate(qua="Expected")%>%mutate(SF="SF2"))
(E2b<-as.tibble(c(0:4))%>%mutate(B=c("B","BB","BBB","BBBB","BBBBB"))%>%mutate(frequency=dgeom(value,36/70))%>%mutate(qua="Expected")%>%mutate(SF="SF2"))
(E3a<-as.tibble(c(0:4))%>%mutate(A=c("A","AA","AAA","AAAA","AAAAA"))%>%mutate(frequency=dgeom(value,18/36))%>%mutate(qua="Expected")%>%mutate(SF="SF3"))
(E3b<-as.tibble(c(0:4))%>%mutate(B=c("B","BB","BBB","BBBB","BBBBB"))%>%mutate(frequency=dgeom(value,18/36))%>%mutate(qua="Expected")%>%mutate(SF="SF3"))
Ea<-E1a%>%bind_rows(E2a,E3a)%>%select(A,frequency,SF,qua)
Ea

Eb<-E1b%>%bind_rows(E2b,E3b)%>%select(B,frequency,SF,qua)
Eb

SFa<-SF1a%>%bind_rows(SF2a,SF3a)%>%mutate(qua="Observed")%>%rename(A=value)
SFa<-SFa%>%bind_rows(tibble(A=c("AAA","AAAA","AAAAA"),frequency=c(0,0,0),SF=c("SF1","SF1","SF1"),qua=c("Observed","Observed","Observed")))%>%bind_rows(tibble(A=c("AAA","AAAA","AAAAA"),frequency=c(0,0,0),SF=c("SF2","SF2","SF2"),qua=c("Observed","Observed","Observed")))%>%bind_rows(tibble(A=c("AAAA","AAAAA"),frequency=c(0,0),SF=c("SF3","SF3"),qua=c("Observed","Observed")))
SFa<-SFa%>%bind_rows(Ea)
g_a<-SFa%>%ggplot(aes(x=A, y=frequency, group=qua, position=qua, fill=qua))
g_a<-g_a+geom_col(position = "dodge")+facet_grid(SF~.)+theme(strip.text.y = element_text(angle = 0))+scale_fill_manual(values=c("#cccccc","#377eb8"))+guides(fill=guide_legend(title=NULL))+theme_classic(base_size=11, base_family='')
print(g_a)
ggsave("continuousAtype.png",g_a, width=4,height=4)

(SFb<-SF1b%>%bind_rows(SF2b,SF3b)%>%mutate(qua="Observed")%>%rename(B=value))
SFb<-SFb%>%bind_rows(tibble(B=c("BBBB","BBBBB"),frequency=c(0,0),SF=c("SF1","SF1"),qua=c("Observed","Observed")))%>%bind_rows(tibble(B=c("BBB","BBBB","BBBBB"),frequency=c(0,0,0),SF=c("SF2","SF2","SF2"),qua=c("Observed","Observed","Observed")))%>%bind_rows(tibble(B=c("BB","BBBB"),frequency=c(0,0),SF=c("SF3","SF3"),qua=c("Observed","Observed")))
SFb<-SFb%>%bind_rows(Eb)
g_b<-SFb%>%ggplot(aes(x=B, y=frequency, group=qua, position=qua, fill=qua))
g_b<-g_b+geom_col(position = "dodge")+facet_grid(SF~.)+theme(strip.text.y = element_text(angle = 0))+scale_fill_manual(values=c("#cccccc","#e41a1c"))+guides(fill=guide_legend(title=NULL))+theme_classic(base_size=11, base_family='')
print(g_b)
ggsave("continuousBtype.png",g_b, width=4,height=4)

#conserved A/B
chr1a <-as.tibble(c())
chr3a <-as.tibble(c('A'))
chr5a <-as.tibble(c())
chr6a <-as.tibble(c())
chr7a <-as.tibble(c())
chr10a <-as.tibble(c())
chr12a <-as.tibble(c())
chr16a <-as.tibble(c())

(SF1a<-chr1a%>%bind_rows(chr3a, chr5a, chr6a, chr7a, chr10a, chr12a, chr16a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF1"))

chr1b <-as.tibble(c('B', 'B'))
chr3b <-as.tibble(c('B', 'B', 'B', 'B', 'B', 'B', 'B'))
chr5b <-as.tibble(c('B', 'B', 'B'))
chr6b <-as.tibble(c('B', 'B', 'B', 'B', 'B', 'B'))
chr7b <-as.tibble(c('B'))
chr10b <-as.tibble(c('BB', 'B'))
chr12b <-as.tibble(c('B', 'B'))
chr16b <-as.tibble(c('B', 'B'))

(SF1b<-chr1b%>%bind_rows(chr3b, chr5b, chr6b, chr7b, chr10b, chr12b, chr16b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF1")) 

chr2a <-as.tibble(c('A', 'A', 'A'))
chr4a <-as.tibble(c('A'))
chr8a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A'))
chr9a <-as.tibble(c('A'))
chr13a <-as.tibble(c('A', 'A', 'A'))
chr15a <-as.tibble(c('A', 'A', 'A', 'A'))
chr18a <-as.tibble(c('A', 'A'))
chr20a <-as.tibble(c('A', 'A'))

(SF2a<-chr2a%>%bind_rows(chr4a, chr8a, chr9a, chr13a, chr15a, chr18a, chr20a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF2")) 

chr2b <-as.tibble(c('B', 'B', 'B', 'B'))
chr4b <-as.tibble(c('BB', 'B'))
chr8b <-as.tibble(c('BB', 'B', 'B', 'B', 'B', 'B'))
chr9b <-as.tibble(c('B'))
chr13b <-as.tibble(c('B', 'B', 'B', 'B', 'B'))
chr15b <-as.tibble(c('B', 'BB', 'B', 'B', 'BB', 'B'))
chr18b <-as.tibble(c('B', 'B', 'B', 'B'))
chr20b <-as.tibble(c('B'))

(SF2b<-chr2b%>%bind_rows(chr4b, chr8b, chr9b, chr13b, chr15b, chr18b, chr20b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF2")) 

chr11a <-as.tibble(c())
chr17a <-as.tibble(c('A', 'A'))
chr17ba <-as.tibble(c('A', 'A'))
chrXa <-as.tibble(c('A', 'A'))

(SF3a<-chr11a%>%bind_rows(chr17a, chr17ba, chrXa)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF3")) 

chr11b <-as.tibble(c())
chr17b <-as.tibble(c())
chr17bb <-as.tibble(c())
chrXb <-as.tibble(c('B', 'B', 'BBB'))

(SF3b<-chr11b%>%bind_rows(chr17b, chr17bb, chrXb)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF3")) 


SFa<-SF1a%>%bind_rows(SF2a,SF3a)%>%mutate(qua="Observed")%>%rename(A=value)
SFa<-SFa%>%bind_rows(tibble(A=c("AA","AAA","AAAA","AAAAA"),frequency=c(0,0,0,0),SF=c("SF1","SF1","SF1","SF1"),qua=c("Observed","Observed","Observed","Observed")))%>%bind_rows(tibble(A=c("AA","AAA","AAAA","AAAAA"),frequency=c(0,0,0,0),SF=c("SF2","SF2","SF2","SF2"),qua=c("Observed","Observed","Observed","Observed")))%>%bind_rows(tibble(A=c("AA","AAA","AAAA","AAAAA"),frequency=c(0,0,0,0),SF=c("SF3","SF3","SF3","SF3"),qua=c("Observed","Observed","Observed","Observed")))
SFa<-SFa%>%bind_rows(Ea)
g_a<-SFa%>%ggplot(aes(x=A, y=frequency, group=qua, position=qua, fill=qua))
g_a<-g_a+geom_col(position = "dodge")+facet_grid(SF~.)+theme(strip.text.y = element_text(angle = 0))+scale_fill_manual(values=c("#cccccc","#377eb8"))+guides(fill=guide_legend(title=NULL))+theme_classic(base_size=11, base_family='')
print(g_a)
ggsave("continuousAbox.png",g_a, width=4,height=4)

SFb<-SF1b%>%bind_rows(SF2b,SF3b)%>%mutate(qua="Observed")%>%rename(B=value)
SFb<-SFb%>%bind_rows(tibble(B=c("BBB","BBBB","BBBBB"),frequency=c(0,0,0),SF=c("SF1","SF1","SF1"),qua=c("Observed","Observed","Observed")))%>%bind_rows(tibble(B=c("BBB","BBBB","BBBBB"),frequency=c(0,0,0),SF=c("SF2","SF2","SF2"),qua=c("Observed","Observed","Observed")))%>%bind_rows(tibble(B=c("BB","BBBB","BBBBB"),frequency=c(0,0,0),SF=c("SF3","SF3","SF3"),qua=c("Observed","Observed","Observed")))
SFb<-SFb%>%bind_rows(Eb)
g_b<-SFb%>%ggplot(aes(x=B, y=frequency, group=qua, position=qua, fill=qua))
g_b<-g_b+geom_col(position = "dodge")+facet_grid(SF~.)+theme(strip.text.y = element_text(angle = 0))+scale_fill_manual(values=c("#cccccc","#e41a1c"))+guides(fill=guide_legend(title=NULL))+theme_classic(base_size=11, base_family='')
print(g_b)

ggsave("continuousBbox.png",g_b, width=4,height=4)
getwd()


#<2 mutation in A/B
chr1a <-as.tibble(c('A', 'A'))
chr3a <-as.tibble(c('A', 'A', 'A', 'AA', 'A', 'A', 'A', 'A'))
chr5a <-as.tibble(c('A', 'A', 'A', 'A'))
chr6a <-as.tibble(c('A', 'A', 'AA', 'A', 'A', 'A', 'A', 'A'))
chr7a <-as.tibble(c('A'))
chr10a <-as.tibble(c('A', 'A'))
chr12a <-as.tibble(c('A', 'A'))
chr16a <-as.tibble(c('A', 'A'))

(SF1a<-chr1a%>%bind_rows(chr3a, chr5a, chr6a, chr7a, chr10a, chr12a, chr16a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF1"))

chr1b <-as.tibble(c('B', 'B'))
chr3b <-as.tibble(c('B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'))
chr5b <-as.tibble(c('B'))
chr6b <-as.tibble(c('BBB', 'B', 'B', 'B', 'B', 'B', 'B'))
chr7b <-as.tibble(c('B'))
chr10b <-as.tibble(c('BB', 'B'))
chr12b <-as.tibble(c('B', 'B'))
chr16b <-as.tibble(c('B', 'B'))

(SF1b<-chr1b%>%bind_rows(chr3b, chr5b, chr6b, chr7b, chr10b, chr12b, chr16b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF1")) 

chr2a <-as.tibble(c('A', 'A', 'A', 'A'))
chr4a <-as.tibble(c('A', 'A', 'A'))
chr8a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A', 'A'))
chr9a <-as.tibble(c('A'))
chr13a <-as.tibble(c('AA', 'A', 'A', 'A', 'A'))
chr15a <-as.tibble(c('A', 'A', 'A', 'A', 'A', 'A'))
chr18a <-as.tibble(c('A', 'A', 'A', 'A', 'A'))
chr20a <-as.tibble(c('A', 'A'))

(SF2a<-chr2a%>%bind_rows(chr4a, chr8a, chr9a, chr13a, chr15a, chr18a, chr20a)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF2")) 

chr2b <-as.tibble(c('B', 'B', 'B', 'B'))
chr4b <-as.tibble(c('B', 'BB', 'B'))
chr8b <-as.tibble(c('BB', 'B', 'B', 'B', 'B', 'B', 'B'))
chr9b <-as.tibble(c('B'))
chr13b <-as.tibble(c('B', 'B', 'B', 'B', 'B'))
chr15b <-as.tibble(c('B', 'BB', 'B', 'B', 'BB', 'B'))
chr18b <-as.tibble(c('B', 'B', 'B', 'B'))
chr20b <-as.tibble(c('BB'))

(SF2b<-chr2b%>%bind_rows(chr4b, chr8b, chr9b, chr13b, chr15b, chr18b, chr20b)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF2")) 

chr11a <-as.tibble(c('AA'))
chr17a <-as.tibble(c('AA', 'A', 'AA', 'A'))
chr17ba <-as.tibble(c('AAA', 'AA'))
chrXa <-as.tibble(c('AA', 'AA', 'A'))

(SF3a<-chr11a%>%bind_rows(chr17a, chr17ba, chrXa)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF3")) 

chr11b <-as.tibble(c('BB'))
chr17b <-as.tibble(c('BB', 'BB', 'B', 'BBB'))
chr17bb <-as.tibble(c('B', 'B'))
chrXb <-as.tibble(c('B', 'B', 'B', 'BBB'))

(SF3b<-chr11b%>%bind_rows(chr17b, chr17bb, chrXb)%>%mutate(c=1)%>%group_by(value)%>%summarise(cnt=sum(c))%>%ungroup()%>%mutate(frequency=cnt/sum(cnt))%>%mutate(SF="SF3")) 


SFa<-SF1a%>%bind_rows(SF2a,SF3a)%>%mutate(qua="Observed")%>%rename(A=value)
SFa<-SFa%>%bind_rows(tibble(A=c("AAA","AAAA","AAAAA","AAAAAA"),cnt=c(0,0,0,0),SF=c("SF1","SF1","SF1","SF1"),qua=c("Observed","Observed","Observed","Observed")))%>%bind_rows(tibble(A=c("AAA","AAAA","AAAAA","AAAAAA"),cnt=c(0,0,0,0),SF=c("2","2","2","2"),qua=c("Observed","Observed","Observed","Observed")))%>%bind_rows(tibble(A=c("AAAA","AAAAA","AAAAAA"),cnt=c(0,0,0),SF=c("3","3","3"),qua=c("Observed","Observed","Observed")))
SFa<-SFa%>%bind_rows(Ea)
SFa
g_a<-SFa%>%ggplot(aes(x=A, y=cnt, group=qua, position=qua, fill=qua))
g_a<-g_a+geom_col(position = "dodge")+facet_grid(SF~.,margins = TRUE)
print(g_a)

SFb<-SF1b%>%bind_rows(SF2b,SF3b)%>%mutate(qua="Observed")%>%rename(B=value)
SFb<-SFb%>%bind_rows(tibble(B=c("BBBB","BBBBB","BBBBBB"),cnt=c(0,0,0),SF=c("1","1","1"),qua=c("Observed","Observed","Observed")))%>%bind_rows(tibble(B=c("BBB","BBBB","BBBBB","BBBBBB"),cnt=c(0,0,0,0),SF=c("2","2","2","2"),qua=c("Observed","Observed","Observed","Observed")))%>%bind_rows(tibble(B=c("BB","BBBB","BBBBBB"),cnt=c(0,0,0),SF=c("3","3","3"),qua=c("Observed","Observed","Observed")))
SFb<-SFb%>%bind_rows(Eb)
g_b<-SFb%>%ggplot(aes(x=B, y=cnt, group=qua, position=qua, fill=qua))
g_b<-g_b+geom_col(position = "dodge")+facet_grid(SF~.,margins = TRUE)
print(g_b)

ggsave("Bbox_save.png", g_b, width = 4, height = 4, dpi = 100)
getwd()
























