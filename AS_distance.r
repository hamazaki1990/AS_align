library(tidyverse)

setwd("~/github/AS_align//")

(pairdist<-read_csv("AS_full.csv")%>%rename("AS"=X1))
gp<-(pairdist%>%select(AS, W3_DXZ1_monomer1_X02418.1))%>%ggplot(.,aes(x=W3_DXZ1_monomer1_X02418.1))+geom_bar()

gp<-pairdist%>%gather(vs, distance, -AS)%>%arrange(AS)%>%ggplot(.,aes(x=AS, y=vs, fill=distance))+geom_tile()
gp<-gp + theme(axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))
print(gp)


chr1<-read_csv("chr1_ASs_labels.csv")%>%unlist()%>%c(.)
chr1
pairdist%>%filter(AS=="W3_DXZ1_monomer1_X02418.1")%>%select(AS,chr1)%>%gather(vs, distance, -AS)

for (i in c(1:22, "X", "Y")){
  csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
  chr<-read_csv(csv)%>%unlist()%>%c(.)
  assign(paste("table",i,sep = ""), pairdist%>%filter(AS=="W3_DXZ1_monomer1_X02418.1")%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i))
}

df<-table1%>%bind_rows(table2, table3, table4, table5, table6, table7, table8)
df
