library(tidyverse)

setwd("~/github/AS_align//")

(pairdist<-read_csv("AS_full.csv")%>%rename("AS"=X1))
gp<-(pairdist%>%select(AS, W3_DXZ1_monomer1_X02418.1))%>%ggplot(.,aes(x=W3_DXZ1_monomer1_X02418.1))+geom_bar()

(pairdist%>%select(`J2_D1Z7/D5Z2/D19Z3mono1_Acc_No._AJ295044.1`))


(pairdist<-pairdist%>%gather(vs, distance, -AS))
(df<-pairdist%>%arrange(AS)%>%arrange(distance))
gp<-df%>%ggplot(.,aes(x=AS, y=vs, fill=distance)) + geom_raster()
gp<-gp+ scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))  + theme(axis.ticks = element_blank(), axis.text.x = element_text(size=5, angle = 280, hjust = 0),axis.text.y = element_text(size=5))
print(gp)

(pairdist<-read_csv("AS_full.csv")%>%rename("AS"=X1))

chr1<-read_csv("chr1_ASs_labels.csv")%>%unlist()%>%c(.)
chr1
pairdist%>%filter(AS=="W3_DXZ1_monomer1_X02418.1")%>%select(AS,chr1)%>%gather(vs, distance, -AS)

chromosomes<-c(1:17, "17b", 18:22,  "X", "Y")

chromosomes

for (i in chromosomes){
  csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
  chr<-read_csv(csv)%>%unlist()%>%c(.)
  assign(paste("table",i,sep = ""), pairdist%>%filter(AS=="W3_DXZ1_monomer1_X02418.1")%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
}
table1

gp<-table1%>%ggplot(.,aes(x=reorder(vs,mono), y=AS, fill=distance))+geom_tile()+coord_fixed(0.7)
print(gp)

df<-table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY)%>%mutate(AS_lab="m1_W3")
df
df$chr_f = factor(df$chromosome, levels=chromosomes)

Xmono<-read_csv("chrX_ASs_labels.csv")%>%unlist()%>%c(.)
Xmono

for (j in Xmono){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfW3_DXZ1_monomer1_X02418.1%>%bind_rows(dfW4_DXZ1_monomer2_X02418.1, dfW5_DXZ1_monomer3_X02418.1, dfW1_DXZ1_monomer4_X02418.1, dfW2_DXZ1_monomer5_X02418.1,dfW3_DXZ1_monomer6_X02418.1, dfW3_DXZ1_monomer6_X02418.1,dfW4_DXZ1_monomer7_X02418.1, dfW5_DXZ1_monomer8_X02418.1,dfW1_DXZ1_monomer9_X02418.1, dfW2_DXZ1_monomer10_X02418.1,dfW3_DXZ1_monomer11_X02418.1,dfW4_DXZ1_monomer12_X02418.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))%>%filter(chromosome=="11")
df1$value[df1$distance_per>0.5]<-""
df1%>%select(value)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=Xmono)
AS_labs<-c("W3_DXZ1_monomer1_X02418.1"="W3_m1","W4_DXZ1_monomer2_X02418.1"="W4_m2","W5_DXZ1_monomer3_X02418.1"="W5_m3","W1_DXZ1_monomer4_X02418.1"="W1_m4","W2_DXZ1_monomer5_X02418.1"="W2_m5","W3_DXZ1_monomer6_X02418.1"="W3_m6","W4_DXZ1_monomer7_X02418.1"="W4_m7","W5_DXZ1_monomer8_X02418.1"="W5_m8","W1_DXZ1_monomer9_X02418.1"="W1_m9","W2_DXZ1_monomer10_X02418.1"="W2_m10","W3_DXZ1_monomer11_X02418.1"="W3_m11","W4_DXZ1_monomer12_X02418.1"="W4_m12")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F) +coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.3))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)


c1mono<-read_csv("chr1_ASs_labels.csv")%>%unlist()%>%c(.)

for (j in c1mono){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-`dfJ2_D1Z7/D5Z2/D19Z3mono1_Acc_No._AJ295044.1` %>%bind_rows(`dfJ1_D1Z7/D5Z2/D19Z3mono2_Acc_No._AJ295044.1` , `dfJ2_D1Z7/D5Z2/D19Z3mono3_Acc_No._AJ295044.1` , `dfJ1_D1Z7/D5Z2/D19Z3mono4_Acc_No._AJ295044.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c1mono)
AS_labs<-c("J2_D1Z7/D5Z2/D19Z3mono1_Acc_No._AJ295044.1"="J2_m1","J1_D1Z7/D5Z2/D19Z3mono2_Acc_No._AJ295044.1"="J1_m2","J2_D1Z7/D5Z2/D19Z3mono3_Acc_No._AJ295044.1"="J2_m3","J1_D1Z7/D5Z2/D19Z3mono4_Acc_No._AJ295044.1"="J1_m4")

gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
print(gp1)


(c2<-read_csv("chr2_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c2){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-dfD2_D2Z1mono1_Acc_No._J04773.1%>%bind_rows(dfD1_D2Z1mono2_Acc_No._J04773.1, dfD2_D2Z1mono3_Acc_No._J04773.1, dfD1_D2Z1mono4_Acc_No._J04773.1, dfD2_D2Z1mono5_Acc_No._J04773.1,dfD1_D2Z1mono6_Acc_No._J04773.1,dfD2_D2Z1mono7_Acc_No._J04773.1, dfD1_D2Z1mono8_Acc_No._J04773.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c2)
AS_labs<-c("D2_D2Z1mono1_Acc_No._J04773.1"="D2_m1","D1_D2Z1mono2_Acc_No._J04773.1"="D1_m2","D2_D2Z1mono3_Acc_No._J04773.1"="D2_m3","D1_D2Z1mono4_Acc_No._J04773.1"="D1_m4","D2_D2Z1mono5_Acc_No._J04773.1"="D2_m5","D1_D2Z1mono6_Acc_No._J04773.1"="D1_m6","D2_D2Z1mono7_Acc_No._J04773.1"="D1_m7","D1_D2Z1mono8_Acc_No._J04773.1"="D1_m8")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

c3<-read_csv("chr3_ASs_labels.csv")%>%unlist()%>%c(.)
c3

for (j in c3){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfJ2_D3Z1mono1_Acc_No._Z12006.1%>%bind_rows(dfJ1_D3Z1mono2_Acc_No._Z12006.1, dfJ2_D3Z1mono3_Acc_No._Z12006.1, dfJ1_D3Z1mono4_Acc_No._Z12006.1, dfJ2_D3Z1mono5_Acc_No._Z12006.1, dfJ1_D3Z1mono6_Acc_No._Z12006.1, dfJ2_D3Z1mono7_Acc_No._Z12006.1, dfJ1_D3Z1mono8_Acc_No._Z12006.1, `dfJ1-74W1_D3Z1mono9_Acc_No._Z12006.1`, dfJ2_D3Z1mono10_Acc_No._Z12006.1,dfJ1_D3Z1mono11_Acc_No._Z12006.1,dfJ2_D3Z1mono12_Acc_No._Z12006.1,dfJ1_D3Z1mono13_Acc_No._Z12006.1,dfJ2_D3Z1mono14_Acc_No._Z12006.1,dfJ1_D3Z1mono15_Acc_No._Z12006.1,dfJ2_D3Z1mono16_Acc_No._Z12006.1,dfJ1_D3Z1mono17_Acc_No._Z12006.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c3)
AS_labs<-c("J2_D3Z1mono1_Acc_No._Z12006.1"="J2_m1","J1_D3Z1mono2_Acc_No._Z12006.1"="J1_m2","J2_D3Z1mono3_Acc_No._Z12006.1"="J2_m3","J1_D3Z1mono4_Acc_No._Z12006.1"="J1_m4","J2_D3Z1mono5_Acc_No._Z12006.1"="J2_m5","J1_D3Z1mono6_Acc_No._Z12006.1"="J1_m6","J2_D3Z1mono7_Acc_No._Z12006.1"="J2_m7","J1_D3Z1mono8_Acc_No._Z12006.1"="J1_m8","J1-74W1_D3Z1mono9_Acc_No._Z12006.1"="J1_m9","J2_D3Z1mono10_Acc_No._Z12006.1"="J2_m10","J1_D3Z1mono11_Acc_No._Z12006.1"="m11","J2_D3Z1mono12_Acc_No._Z12006.1"="m12", "J1_D3Z1mono13_Acc_No._Z12006.1"="m13","J2_D3Z1mono14_Acc_No._Z12006.1"="m14","J1_D3Z1mono15_Acc_No._Z12006.1"="m15","J2_D3Z1mono16_Acc_No._Z12006.1"="m16","J1_D3Z1mono17_Acc_No._Z12006.1"="m17")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)


c4<-read_csv("chr4_ASs_labels.csv")%>%unlist()%>%c(.)
c4

for (j in c4){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-dfD1_D4Z1mono1_Acc_No._Z12011.1%>%bind_rows(dfD2_D4Z1mono2_Acc_No._Z12011.1, dfD1_D4Z1mono3_Acc_No._Z12011.1, dfD1_D4Z1mono4_Acc_No._Z12011.1, dfD2_D4Z1mono5_Acc_No._Z12011.1,dfD1_D4Z1mono6_Acc_No._Z12011.1,dfD2_D4Z1mono7_Acc_No._Z12011.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c4)
AS_labs<-c("D1_D4Z1mono1_Acc_No._Z12011.1"="D1_m1","D2_D4Z1mono2_Acc_No._Z12011.1"="D2_m2","D1_D4Z1mono3_Acc_No._Z12011.1"="D3_m3","D1_D4Z1mono4_Acc_No._Z12011.1"="m4","D2_D4Z1mono5_Acc_No._Z12011.1"="m5","D1_D4Z1mono6_Acc_No._Z12011.1"="m6","D2_D4Z1mono7_Acc_No._Z12011.1"="m7")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)


(c5<-read_csv("chr5_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c5){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-`dfJ2_D5Z2/D19Z3mono1__Acc_No._M26919.1`%>%bind_rows(`dfJ1_D5Z2/D19Z3mono2_Acc_No._M26919.1`)

df1
df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c("J2_D5Z2/D19Z3mono1__Acc_No._M26919.1","J1_D5Z2/D19Z3mono2_Acc_No._M26919.1"))
AS_labs<-c("J2_D5Z2/D19Z3mono1__Acc_No._M26919.1"="m1","J1_D5Z2/D19Z3mono2_Acc_No._M26919.1"="m2")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)


(c6<-read_csv("chr6_ASs_labels.csv")%>%unlist()%>%c(.))
for (j in c6){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfD1_D6Z1mono1_Acc_No_GJ211907.1%>%bind_rows(`dfD1(D2)_D6Z1mono2_Acc_No_GJ211907.1`, `dfD2_D6Z1mono3_Acc_No_GJ211907.1`, `dfD1(D2)_D6Z1mono4_Acc_No_GJ211907.1`, `dfD2(D1)_D6Z1mono5_Acc_No_GJ211907.1`, `dfD1(D2)_D6Z1mono6_Acc_No_GJ211907.1`, `dfD2_D6Z1mono7_Acc_No_GJ211907.1`, `dfD2(D1)_D6Z1mono8_Acc_No_GJ211907.1`, `dfD1(D2)_D6Z1mono9_Acc_No_GJ211907.1`, dfD2_D6Z1mono10_Acc_No_GJ211907.1,`dfD1(D2)_D6Z1mono11_Acc_No_GJ211907.1`,dfD2_D6Z1mono12_Acc_No_GJ211907.1,dfD1_D6Z1mono13_Acc_No_GJ211907.1,dfD2_D6Z1mono14_Acc_No_GJ211907.1,dfD1_D6Z1mono15_Acc_No_GJ211907.1,dfD2_D6Z1mono16_Acc_No_GJ211907.1,dfD1_D6Z1mono17_Acc_No_GJ211907.1, dfD2_D6Z1mono18_Acc_No_GJ211907.1,dfD1_D6Z1mono19_Acc_No_GJ211907.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c6)
AS_labs<-c("D1_D6Z1mono1_Acc_No_GJ211907.1"="m1","D1(D2)_D6Z1mono2_Acc_No_GJ211907.1"="m2","D2_D6Z1mono3_Acc_No_GJ211907.1"="m3","D1(D2)_D6Z1mono4_Acc_No_GJ211907.1"="m4","D2(D1)_D6Z1mono5_Acc_No_GJ211907.1"="m5","D1(D2)_D6Z1mono6_Acc_No_GJ211907.1"="m6","D2_D6Z1mono7_Acc_No_GJ211907.1"="m7","D2(D1)_D6Z1mono8_Acc_No_GJ211907.1"="m8","D1(D2)_D6Z1mono9_Acc_No_GJ211907.1"="m9","D2_D6Z1mono10_Acc_No_GJ211907.1"="m10","D1(D2)_D6Z1mono11_Acc_No_GJ211907.1"="m11","D2_D6Z1mono12_Acc_No_GJ211907.1"="m12", "D1_D6Z1mono13_Acc_No_GJ211907.1"="m13","D2_D6Z1mono14_Acc_No_GJ211907.1"="m14","D1_D6Z1mono15_Acc_No_GJ211907.1"="m15","D2_D6Z1mono16_Acc_No_GJ211907.1"="m16","D1_D6Z1mono17_Acc_No_GJ211907.1"="m17", "D2_D6Z1mono18_Acc_No_GJ211907.1"="m18","D1_D6Z1mono19_Acc_No_GJ211907.1"="m19")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c7<-read_csv("chr7_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c7){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-`dfJ1_D7Z1mono1_Acc_No._AC142529.3`%>%bind_rows(`dfJ2_D7Z1mono2_Acc_No._AC142529.3`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c7)
AS_labs<-c("J1_D7Z1mono1_Acc_No._AC142529.3"="m1","J2_D7Z1mono2_Acc_No._AC142529.3"="m2")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c8<-read_csv("chr8_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c8){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfD1_D8Z2mono1_Acc_No._M64779.1%>%bind_rows(`dfD2_D8Z2mono2_Acc_No._M64779.1`, `dfD1_D8Z2mono3_Acc_No._M64779.1`, `dfD2_D8Z2mono4_Acc_No._M64779.1`, `dfD1_D8Z2mono5_Acc_No._M64779.1`, `dfD2_D8Z2mono6_Acc_No._M64779.1`, `dfD1_D8Z2mono7_Acc_No._M64779.1`, `dfD2_D8Z2mono8_Acc_No._M64779.1`, `dfD1_D8Z2mono9_Acc_No._M64779.1`, dfD2_D8Z2mono10_Acc_No._M64779.1,`dfD1_D8Z2mono11_Acc_No._M64779.1`,dfD2_D8Z2mono12_Acc_No._M64779.1,dfD1_D8Z2mono13_Acc_No._M64779.1,dfD2_D8Z2mono14_Acc_No._M64779.1,dfD1_D8Z2mono15_Acc_No._M64779.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c8)
AS_labs<-c("D1_D8Z2mono1_Acc_No._M64779.1"="m1","D2_D8Z2mono2_Acc_No._M64779.1"="m2","D1_D8Z2mono3_Acc_No._M64779.1"="m3","D2_D8Z2mono4_Acc_No._M64779.1"="m4","D1_D8Z2mono5_Acc_No._M64779.1"="m5","D2_D8Z2mono6_Acc_No._M64779.1"="m6","D1_D8Z2mono7_Acc_No._M64779.1"="m7","D2_D8Z2mono8_Acc_No._M64779.1"="m8","D1_D8Z2mono9_Acc_No._M64779.1"="m9","D2_D8Z2mono10_Acc_No._M64779.1"="m10","D1_D8Z2mono11_Acc_No._M64779.1"="m11","D2_D8Z2mono12_Acc_No._M64779.1"="m12", "D1_D8Z2mono13_Acc_No._M64779.1"="m13","D2_D8Z2mono14_Acc_No._M64779.1"="m14","D1_D8Z2mono15_Acc_No._M64779.1"="m15")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

  
(c9<-read_csv("chr9_ASs_labels.csv")%>%unlist()%>%c(.))
for (j in c9){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-`dfD2_D9Z4mono1_Acc_No._M64320.1`%>%bind_rows(`dfD1_D9Z4mono2_Acc_No._M64320.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c9)
AS_labs<-c("D2_D9Z4mono1_Acc_No._M64320.1"="m1","D1_D9Z4mono2_Acc_No._M64320.1"="m2")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)


(c10<-read_csv("chr10_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c10){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-dfJ2_D10Z1mono1_Acc_No._X63622.1%>%bind_rows(dfJ1_D10Z1mono2_Acc_No._X63622.1, dfJ2_D10Z1mono3_Acc_No._X63622.1, dfJ1_D10Z1mono4_Acc_No._X63622.1, dfJ2_D10Z1mono5_Acc_No._X63622.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c10)
AS_labs<-c("J2_D10Z1mono1_Acc_No._X63622.1"="m1","J1_D10Z1mono2_Acc_No._X63622.1"="m2","J2_D10Z1mono3_Acc_No._X63622.1"="m3","J1_D10Z1mono4_Acc_No._X63622.1"="m4","J2_D10Z1mono5_Acc_No._X63622.1"="m5")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c11<-read_csv("chr11_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c11){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-dfW4_D11Z1mono1_Acc_No._M21452.1%>%bind_rows(dfW5_D11Z1mono2_Acc_No._M21452.1, dfW1_D11Z1mono3_Acc_No._M21452.1, dfW2_D11Z1mono4_Acc_No._M21452.1, dfW3_D11Z1mono5_Acc_No._M21452.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))
df1$value[df1$distance_per>0.2]<-""
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c11)
AS_labs<-c("W4_D11Z1mono1_Acc_No._M21452.1"="m1","W5_D11Z1mono2_Acc_No._M21452.1"="m2","W1_D11Z1mono3_Acc_No._M21452.1"="m3","W2_D11Z1mono4_Acc_No._M21452.1"="m4","W3_D11Z1mono5_Acc_No._M21452.1"="m5")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F)+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)


(c12<-read_csv("chr12_ASs_labels.csv")%>%unlist()%>%c(.))
for (j in c12){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}

df1<-dfJ1_D12Z3mono1_Acc_No._M28221.1%>%bind_rows(dfJ2_D12Z3mono2_Acc_No._M28221.1, dfJ1_D12Z3mono3_Acc_No._M28221.1, dfJ2_D12Z3mono4_Acc_No._M28221.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c12)
AS_labs<-c("J1_D12Z3mono1_Acc_No._M28221.1"="m1","J2_D12Z3mono2_Acc_No._M28221.1"="m2","J1_D12Z3mono3_Acc_No._M28221.1"="m3","J2_D12Z3mono4_Acc_No._M28221.1"="m4")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c13<-read_csv("chr13_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c13){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-`dfD1_D13Z1/D21Z1mono1_Acc_No._D29750.1`%>%bind_rows(`dfD2_D13Z1/D21Z1mono2_Acc_No._D29750.1`, `dfD2_D13Z1/D21Z1nmono3_Acc_No._D29750.1`, `dfD1_D13Z1/D21Z1mono4_Acc_No._D29750.1`, `dfD2_D13Z1/D21Z1mono5_Acc_No._D29750.1`, `dfD1_D13Z1/D21Z1mono6_Acc_No._D29750.1`, `dfD2_D13Z1/D21Z1mono7_Acc_No._D29750.1`, `dfD1_D13Z1/D21Z1mono8_Acc_No._D29750.1`, `dfD2_D13Z1/D21Z1mono9_Acc_No._D29750.1`, `dfD1_D13Z1/D21Z1mono10_Acc_No._D29750.1`,`dfD2_D13Z1/D21Z1mono11_Acc_No._D29750.1`)

df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c13)
AS_labs<-c("D1_D13Z1/D21Z1mono1_Acc_No._D29750.1"="m1","D2_D13Z1/D21Z1mono2_Acc_No._D29750.1"="m2","D2_D13Z1/D21Z1nmono3_Acc_No._D29750.1"="m3","D1_D13Z1/D21Z1mono4_Acc_No._D29750.1"="m4","D2_D13Z1/D21Z1mono5_Acc_No._D29750.1"="m5","D1_D13Z1/D21Z1mono6_Acc_No._D29750.1"="m6","D2_D13Z1/D21Z1mono7_Acc_No._D29750.1"="m7","D1_D13Z1/D21Z1mono8_Acc_No._D29750.1"="m8","D2_D13Z1/D21Z1mono9_Acc_No._D29750.1"="m9","D1_D13Z1/D21Z1mono10_Acc_No._D29750.1"="m10","D2_D13Z1/D21Z1mono11_Acc_No._D29750.1"="m11")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c14<-read_csv("chr14_ASs_labels.csv")%>%unlist()%>%c(.))

for (i in chromosomes){
  csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
  chr<-read_csv(csv)%>%unlist()%>%c(.)
  assign(paste("table",i,sep = ""), pairdist%>%filter(AS=="D1_D14Z1/D22Z1_Acc_No._M22273.1")%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
}
df1<-table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY)

df1
df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c14)
AS_labs<-c("D1_D14Z1/D22Z1_Acc_No._M22273.1" ="m1")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c15<-read_csv("chr15_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c15){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfD2_D15Z3mono1_Acc_No._GJ212045.1%>%bind_rows(`dfD1_D15Z3mono2_Acc_No._GJ212045.1`, `dfD2_D15Z3mono3_Acc_No._GJ212045.1`, `dfD1_D15Z3mono4_Acc_No._GJ212045.1`, `dfD1_D15Z3mono5_Acc_No._GJ212045.1`, `dfD2_D15Z3mono6_Acc_No._GJ212045.1`, `dfD1_D15Z3mono7_Acc_No._GJ212045.1`, `dfD2_D15Z3mono8_Acc_No._GJ212045.1`, `dfD1_D15Z3mono9_Acc_No._GJ212045.1`, dfD2_D15Z3mono10_Acc_No._GJ212045.1,`dfD1(D2)_D15Z3mono11_Acc_No._GJ212045.1`,dfD1_D15Z3mono12_Acc_No._GJ212045.1,dfD2_D15Z3mono13_Acc_No._GJ212045.1,dfD1_D15Z3mono14_Acc_No._GJ212045.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c15)
AS_labs<-c("D2_D15Z3mono1_Acc_No._GJ212045.1"="m1","D1_D15Z3mono2_Acc_No._GJ212045.1"="m2","D2_D15Z3mono3_Acc_No._GJ212045.1"="m3","D1_D15Z3mono4_Acc_No._GJ212045.1"="m4","D1_D15Z3mono5_Acc_No._GJ212045.1"="m5","D2_D15Z3mono6_Acc_No._GJ212045.1"="m6","D1_D15Z3mono7_Acc_No._GJ212045.1"="m7","D2_D15Z3mono8_Acc_No._GJ212045.1"="m8","D1_D15Z3mono9_Acc_No._GJ212045.1"="m9","D2_D15Z3mono10_Acc_No._GJ212045.1"="m10","D1(D2)_D15Z3mono11_Acc_No._GJ212045.1"="m11","D1_D15Z3mono12_Acc_No._GJ212045.1"="m12", "D2_D15Z3mono13_Acc_No._GJ212045.1"="m13","D1_D15Z3mono14_Acc_No._GJ212045.1"="m14")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c16<-read_csv("chr16_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c16){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-`dfJ2_D16Z2mono1_Acc_No._M58446.1` %>%bind_rows(`dfJ1_D16Z2mono2_Acc_No._M58446.1` , `dfJ2_D16Z2mono3_Acc_No._M58446.1` , `dfJ1_D16Z2mono4_Acc_No._M58446.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c16)
AS_labs<-c("J2_D16Z2mono1_Acc_No._M58446.1"="m1","J1_D16Z2mono2_Acc_No._M58446.1"="m2","J2_D16Z2mono3_Acc_No._M58446.1"="m3","J1_D16Z2mono4_Acc_No._M58446.1"="m4")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c17<-read_csv("chr17_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c17){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfW1_D17Z1mono1_Acc_No._M13882.1%>%bind_rows(`dfW2_D17Z1mono2_Acc_No._M13882.1`, `dfW3_D17Z1mono3_Acc_No._M13882.1`, `dfW4_D17Z1mono4_Acc_No._M13882.1`, `dfW5_D17Z1mono5_Acc_No._M13882.1`, `dfW1_D17Z1mono6_Acc_No._M13882.1`, `dfW2_D17Z1mono7_Acc_No._M13882.1`, `dfW3_D17Z1mono8_Acc_No._M13882.1`, `dfW4_D17Z1mono9_Acc_No._M13882.1`, dfW3_D17Z1mono10_Acc_No._M13882.1,`dfW4_D17Z1mono11_Acc_No._M13882.1`,dfW5_D17Z1mono12_Acc_No._M13882.1,dfW1_D17Z1mono13_Acc_No._M13882.1,dfW1_D17Z1mono14_Acc_No._M13882.1,dfW1_D17Z1mono15_Acc_No._M13882.1,dfW5_D17Z1mono16_Acc_No._M13882.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c17)
AS_labs<-c("W1_D17Z1mono1_Acc_No._M13882.1"="m1","W2_D17Z1mono2_Acc_No._M13882.1"="m2","W3_D17Z1mono3_Acc_No._M13882.1"="m3","W4_D17Z1mono4_Acc_No._M13882.1"="m4","W5_D17Z1mono5_Acc_No._M13882.1"="m5","W1_D17Z1mono6_Acc_No._M13882.1"="m6","W2_D17Z1mono7_Acc_No._M13882.1"="m7","W3_D17Z1mono8_Acc_No._M13882.1"="m8","W4_D17Z1mono9_Acc_No._M13882.1"="m9","W3_D17Z1mono10_Acc_No._M13882.1"="m10","W4_D17Z1mono11_Acc_No._M13882.1"="m11","W5_D17Z1mono12_Acc_No._M13882.1"="m12", "W1_D17Z1mono13_Acc_No._M13882.1"="m13","W1_D17Z1mono14_Acc_No._M13882.1"="m14","W1_D17Z1mono15_Acc_No._M13882.1"="m15","W5_D17Z1mono16_Acc_No._M13882.1"="m16","D1_D6Z1mono17_Acc_No_GJ211907.1"="m17", "D2_D6Z1mono18_Acc_No_GJ211907.1"="m18","D1_D6Z1mono19_Acc_No_GJ211907.1"="m19")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c17b<-read_csv("chr17b_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c17b){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfW4_D17Z1bmono1_Acc_No._GJ212053.1%>%bind_rows(`dfW5_D17Z1bmono2_Acc_No._GJ212053.1`, `dfW1_D17Z1bmono3_Acc_No._GJ212053.1`, `dfW2_D17Z1bmono4_Acc_No._GJ212053.1`, `dfW3_D17Z1bmono5_Acc_No._GJ212053.1`, `dfW2_D17Z1bmono6_Acc_No._GJ212053.1`, `dfW3_D17Z1bmono7_Acc_No._GJ212053.1`, `dfW4_D17Z1bmono8_Acc_No._GJ212053.1`, `dfW5_D17Z1bmono9_Acc_No._GJ212053.1`, dfW1_D17Z1bmono10_Acc_No._GJ212053.1,`dfW5_D17Z1bmono11_Acc_No._GJ212053.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c17b)
AS_labs<-c("W4_D17Z1bmono1_Acc_No._GJ212053.1"="m1","W5_D17Z1bmono2_Acc_No._GJ212053.1"="m2","W1_D17Z1bmono3_Acc_No._GJ212053.1"="m3","W2_D17Z1bmono4_Acc_No._GJ212053.1"="m4","W3_D17Z1bmono5_Acc_No._GJ212053.1"="m5","W2_D17Z1bmono6_Acc_No._GJ212053.1"="m6","W3_D17Z1bmono7_Acc_No._GJ212053.1"="m7","W4_D17Z1bmono8_Acc_No._GJ212053.1"="m8","W5_D17Z1bmono9_Acc_No._GJ212053.1"="m9","W1_D17Z1bmono10_Acc_No._GJ212053.1"="m10","W5_D17Z1bmono11_Acc_No._GJ212053.1"="m11")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c18<-read_csv("chr18_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c18){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfD1_D18Z1mono1_Acc_No._M65181.1%>%bind_rows(`dfD2_D18Z1mono2_Acc_No._M65181.1`, `dfD1_D18Z1mono3_Acc_No._M65181.1`, `dfD2_D18Z1mono4_Acc_No._M65181.1`, `dfD1_D18Z1mono5_Acc_No._M65181.1`, `dfD2_D18Z1mono6_Acc_No._M65181.1`, `dfD1_D18Z1mono7_Acc_No._M65181.1`, `dfD2_D18Z1mono8_Acc_No._M65181.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c18)
AS_labs<-c("D1_D18Z1mono1_Acc_No._M65181.1"="m1","D2_D18Z1mono2_Acc_No._M65181.1"="m2","D1_D18Z1mono3_Acc_No._M65181.1"="m3","D2_D18Z1mono4_Acc_No._M65181.1"="m4","D1_D18Z1mono5_Acc_No._M65181.1"="m5","D2_D18Z1mono6_Acc_No._M65181.1"="m6","D1_D18Z1mono7_Acc_No._M65181.1"="m7")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c19<-read_csv("chr19_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c19){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-`dfJ2_D5Z1/D19Z3mono1_Acc_no._M26920.1`%>%bind_rows(`dfJ1_D5Z1/D19Z3mono2_Acc_no._M26920.1`, `dfJ2_D5Z1/D19Z3mono3_Acc_no._M26920.1`, `dfJ1_D5Z1/D19Z3mono4_Acc_no._M26920.1`, `dfJ2_D5Z1/D19Z3mono5_Acc_no._M26920.1`, `dfJ1_D5Z1/D19Z3mono6_Acc_no._M26920.1`, `dfJ2_D5Z2/D19Z3mono1__Acc_No._M26919.1`, `dfJ1_D5Z2/D19Z3mono2_Acc_No._M26919.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c19)
AS_labs<-c("J2_D5Z1/D19Z3mono1_Acc_no._M26920.1"="m1","J1_D5Z1/D19Z3mono2_Acc_no._M26920.1"="m2","J2_D5Z1/D19Z3mono3_Acc_no._M26920.1"="m3","J1_D5Z1/D19Z3mono4_Acc_no._M26920.1"="m4","J2_D5Z1/D19Z3mono5_Acc_no._M26920.1"="m5","J1_D5Z1/D19Z3mono6_Acc_no._M26920.1"="m6","J2_D5Z2/D19Z3mono1__Acc_No._M26919.1"="m7","J1_D5Z2/D19Z3mono2_Acc_No._M26919.1"="m8")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(c20<-read_csv("chr20_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in c20){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-`dfD1_D20Z2mono1_Acc_No._X58269.1_X56450.1`%>%bind_rows(`dfD2_D20Z2mono2_Acc_No._X58269.1_X56450.1`, `dfD1_D20Z2mono3_Acc_No._X58269.1_X56450.1`, `dfD2_D20Z2mono4_Acc_No._X58269.1_X56450.1`, `dfD1_D20Z2mono5_Acc_No._X58269.1_X56450.1`)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=c20)
AS_labs<-c("D1_D20Z2mono1_Acc_No._X58269.1_X56450.1"="m1","D2_D20Z2mono2_Acc_No._X58269.1_X56450.1"="m2","D1_D20Z2mono3_Acc_No._X58269.1_X56450.1"="m3","D2_D20Z2mono4_Acc_No._X58269.1_X56450.1"="m4","D1_D20Z2mono5_Acc_No._X58269.1_X56450.1"="m5")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)

(cY<-read_csv("chrY_ASs_labels.csv")%>%unlist()%>%c(.))

for (j in cY){
  for (i in chromosomes){
    csv<-paste("chr",i,"_ASs_labels.csv", sep = "")
    chr<-read_csv(csv)%>%unlist()%>%c(.)
    assign(paste("table",i,sep = ""), pairdist%>%filter(AS==paste(j))%>%select(AS,chr)%>%gather(vs, distance, -AS)%>%mutate(chromosome=i)%>%mutate(mono=c(1:n())))
  }
  assign(paste("df",j,sep=""),table1%>%bind_rows(table2,table3,table4,table5,table6,table7,table8,table9,table10,table11,table12,table13,table14,table15,table16,table17,table17b,table18,table19,table20,table21,table22,tableX,tableY))
}
df1<-dfM1_DYZ3mono1_Acc_No._M29723.1%>%bind_rows(`dfM1_DYZ3mono2_Acc_No._M29723.1`, `dfM1_DYZ3mono3_Acc_No._M29723.1`, `dfM1_DYZ3mono4_Acc_No._M29723.1`, `dfM1_DYZ3mono5_Acc_No._M29723.1`, `dfM1_DYZ3mono6_Acc_No._M29723.1`, `dfM1_DYZ3mono7_Acc_No._M29723.1`, `dfM1_DYZ3mono8_Acc_No._M29723.1`, `dfM1_DYZ3mono9_Acc_No._M29723.1`, dfM1_DYZ3mono10_Acc_No._M29723.1,`dfM1_DYZ3mono11_Acc_No._M29723.1`,dfM1_DYZ3mono12_Acc_No._M29723.1,dfM1_DYZ3mono13_Acc_No._M29723.1,dfM1_DYZ3mono14_Acc_No._M29723.1)

df1<-df1%>%mutate(distance_per=distance/146)
df1$chr_f = factor(df1$chromosome, levels=chromosomes)
df1$AS_f = factor(df1$AS, levels=cY)
AS_labs<-c("M1_DYZ3mono1_Acc_No._M29723.1"="m1","M1_DYZ3mono2_Acc_No._M29723.1"="m2","M1_DYZ3mono3_Acc_No._M29723.1"="m3","M1_DYZ3mono4_Acc_No._M29723.1"="m4","M1_DYZ3mono5_Acc_No._M29723.1"="m5","M1_DYZ3mono6_Acc_No._M29723.1"="m6","M1_DYZ3mono7_Acc_No._M29723.1"="m7","M1_DYZ3mono8_Acc_No._M29723.1"="m8","M1_DYZ3mono9_Acc_No._M29723.1"="m9","M1_DYZ3mono10_Acc_No._M29723.1"="m10","M1_DYZ3mono11_Acc_No._M29723.1"="m11","M1_DYZ3mono12_Acc_No._M29723.1"="m12", "M1_DYZ3mono13_Acc_No._M29723.1"="m13","M1_DYZ3mono14_Acc_No._M29723.1"="m14")

gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+coord_fixed(10)+facet_grid(chr_f~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), chr_f=as_labeller(df1$chr_f)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.25))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+guides(fill="none")+theme(strip.text.y = element_text(angle=180))
print(gp1)





df1<-dfW3_DXZ1_monomer1_X02418.1%>%bind_rows(dfW4_DXZ1_monomer2_X02418.1, dfW5_DXZ1_monomer3_X02418.1, dfW1_DXZ1_monomer4_X02418.1, dfW2_DXZ1_monomer5_X02418.1,dfW3_DXZ1_monomer6_X02418.1, dfW3_DXZ1_monomer6_X02418.1,dfW4_DXZ1_monomer7_X02418.1, dfW5_DXZ1_monomer8_X02418.1,dfW1_DXZ1_monomer9_X02418.1, dfW2_DXZ1_monomer10_X02418.1,dfW3_DXZ1_monomer11_X02418.1,dfW4_DXZ1_monomer12_X02418.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))%>%filter(chromosome=="11")
df1$value[df1$distance_per>0.5]<-""
df1%>%select(value)
df1$AS_f = factor(df1$AS, levels=Xmono)
AS_labs<-c("W3_DXZ1_monomer1_X02418.1"="W3_m1","W4_DXZ1_monomer2_X02418.1"="W4_m2","W5_DXZ1_monomer3_X02418.1"="W5_m3","W1_DXZ1_monomer4_X02418.1"="W1_m4","W2_DXZ1_monomer5_X02418.1"="W2_m5","W3_DXZ1_monomer6_X02418.1"="W3_m6","W4_DXZ1_monomer7_X02418.1"="W4_m7","W5_DXZ1_monomer8_X02418.1"="W5_m8","W1_DXZ1_monomer9_X02418.1"="W1_m9","W2_DXZ1_monomer10_X02418.1"="W2_m10","W3_DXZ1_monomer11_X02418.1"="W3_m11","W4_DXZ1_monomer12_X02418.1"="W4_m12")
chr_lab<-c("1"="W4_m1", "2"="W5_m5", "3"="W1_m3", "4"="W2_m4", "5"="W3_m5")
gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F) +coord_fixed(10)+facet_grid(mono~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), mono=as_labeller(chr_lab)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.3))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)


df1<-dfW3_DXZ1_monomer1_X02418.1%>%bind_rows(dfW4_DXZ1_monomer2_X02418.1, dfW5_DXZ1_monomer3_X02418.1, dfW1_DXZ1_monomer4_X02418.1, dfW2_DXZ1_monomer5_X02418.1,dfW3_DXZ1_monomer6_X02418.1, dfW3_DXZ1_monomer6_X02418.1,dfW4_DXZ1_monomer7_X02418.1, dfW5_DXZ1_monomer8_X02418.1,dfW1_DXZ1_monomer9_X02418.1, dfW2_DXZ1_monomer10_X02418.1,dfW3_DXZ1_monomer11_X02418.1,dfW4_DXZ1_monomer12_X02418.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))%>%filter(chromosome=="17")
df1$value[df1$distance_per>0.5]<-""
df1%>%select(value)
df1$AS_f = factor(df1$AS, levels=Xmono)
AS_labs<-c("W3_DXZ1_monomer1_X02418.1"="W3_m1","W4_DXZ1_monomer2_X02418.1"="W4_m2","W5_DXZ1_monomer3_X02418.1"="W5_m3","W1_DXZ1_monomer4_X02418.1"="W1_m4","W2_DXZ1_monomer5_X02418.1"="W2_m5","W3_DXZ1_monomer6_X02418.1"="W3_m6","W4_DXZ1_monomer7_X02418.1"="W4_m7","W5_DXZ1_monomer8_X02418.1"="W5_m8","W1_DXZ1_monomer9_X02418.1"="W1_m9","W2_DXZ1_monomer10_X02418.1"="W2_m10","W3_DXZ1_monomer11_X02418.1"="W3_m11","W4_DXZ1_monomer12_X02418.1"="W4_m12")
chr_lab<-c("1"="W1_m1", "2"="W2_m5", "3"="W3_m3", "4"="W4_m4", "5"="W5_m5", "6"="W1_m6", "7"="W2_m7", "8"="W3_m8", "9"="W4_m9","10"="W3_m10", "11"="W4_m11", "12"="W5_m12", "13"="W1_m13", "14"="W1_m14", "15"="W1_m15", "16"="W5_m16")
gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F) +coord_fixed(10)+facet_grid(mono~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), mono=as_labeller(chr_lab)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.3))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)


df1<-dfW3_DXZ1_monomer1_X02418.1%>%bind_rows(dfW4_DXZ1_monomer2_X02418.1, dfW5_DXZ1_monomer3_X02418.1, dfW1_DXZ1_monomer4_X02418.1, dfW2_DXZ1_monomer5_X02418.1,dfW3_DXZ1_monomer6_X02418.1, dfW3_DXZ1_monomer6_X02418.1,dfW4_DXZ1_monomer7_X02418.1, dfW5_DXZ1_monomer8_X02418.1,dfW1_DXZ1_monomer9_X02418.1, dfW2_DXZ1_monomer10_X02418.1,dfW3_DXZ1_monomer11_X02418.1,dfW4_DXZ1_monomer12_X02418.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))%>%filter(chromosome=="17")
df1$value[df1$distance_per>0.5]<-""
df1%>%select(value)
df1$AS_f = factor(df1$AS, levels=Xmono)
AS_labs<-c("W3_DXZ1_monomer1_X02418.1"="W3_m1","W4_DXZ1_monomer2_X02418.1"="W4_m2","W5_DXZ1_monomer3_X02418.1"="W5_m3","W1_DXZ1_monomer4_X02418.1"="W1_m4","W2_DXZ1_monomer5_X02418.1"="W2_m5","W3_DXZ1_monomer6_X02418.1"="W3_m6","W4_DXZ1_monomer7_X02418.1"="W4_m7","W5_DXZ1_monomer8_X02418.1"="W5_m8","W1_DXZ1_monomer9_X02418.1"="W1_m9","W2_DXZ1_monomer10_X02418.1"="W2_m10","W3_DXZ1_monomer11_X02418.1"="W3_m11","W4_DXZ1_monomer12_X02418.1"="W4_m12")
chr_lab<-c("1"="W1_m1", "2"="W2_m5", "3"="W3_m3", "4"="W4_m4", "5"="W5_m5", "6"="W1_m6", "7"="W2_m7", "8"="W3_m8", "9"="W4_m9","10"="W3_m10", "11"="W4_m11", "12"="W5_m12", "13"="W1_m13", "14"="W1_m14", "15"="W1_m15", "16"="W5_m16")
gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F) +coord_fixed(10)+facet_grid(mono~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), mono=as_labeller(chr_lab)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.3))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)


df1<-dfW3_DXZ1_monomer1_X02418.1%>%bind_rows(dfW4_DXZ1_monomer2_X02418.1, dfW5_DXZ1_monomer3_X02418.1, dfW1_DXZ1_monomer4_X02418.1, dfW2_DXZ1_monomer5_X02418.1,dfW3_DXZ1_monomer6_X02418.1, dfW3_DXZ1_monomer6_X02418.1,dfW4_DXZ1_monomer7_X02418.1, dfW5_DXZ1_monomer8_X02418.1,dfW1_DXZ1_monomer9_X02418.1, dfW2_DXZ1_monomer10_X02418.1,dfW3_DXZ1_monomer11_X02418.1,dfW4_DXZ1_monomer12_X02418.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))%>%filter(chromosome=="X")
df1$value[df1$distance_per>0.5]<-""
df1%>%select(value)
df1$AS_f = factor(df1$AS, levels=Xmono)
AS_labs<-c("W3_DXZ1_monomer1_X02418.1"="W3_m1","W4_DXZ1_monomer2_X02418.1"="W4_m2","W5_DXZ1_monomer3_X02418.1"="W5_m3","W1_DXZ1_monomer4_X02418.1"="W1_m4","W2_DXZ1_monomer5_X02418.1"="W2_m5","W3_DXZ1_monomer6_X02418.1"="W3_m6","W4_DXZ1_monomer7_X02418.1"="W4_m7","W5_DXZ1_monomer8_X02418.1"="W5_m8","W1_DXZ1_monomer9_X02418.1"="W1_m9","W2_DXZ1_monomer10_X02418.1"="W2_m10","W3_DXZ1_monomer11_X02418.1"="W3_m11","W4_DXZ1_monomer12_X02418.1"="W4_m12")
chr_lab<-c("1"="W3_m1", "2"="W4_m5", "3"="W5_m3", "4"="W1_m4", "5"="W2_m5","6"="W3_m6", "7"="W4_m7", "8"="W5_m8", "9"="W1_m9", "10"="W2_m10","11"="W3_m11", "12"="W4_m12")
gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F) +coord_fixed(10)+facet_grid(mono~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), mono=as_labeller(chr_lab)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.3))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)


df1<-dfW3_DXZ1_monomer1_X02418.1%>%bind_rows(dfW4_DXZ1_monomer2_X02418.1, dfW5_DXZ1_monomer3_X02418.1, dfW1_DXZ1_monomer4_X02418.1, dfW2_DXZ1_monomer5_X02418.1,dfW3_DXZ1_monomer6_X02418.1, dfW3_DXZ1_monomer6_X02418.1,dfW4_DXZ1_monomer7_X02418.1, dfW5_DXZ1_monomer8_X02418.1,dfW1_DXZ1_monomer9_X02418.1, dfW2_DXZ1_monomer10_X02418.1,dfW3_DXZ1_monomer11_X02418.1,dfW4_DXZ1_monomer12_X02418.1)

df1
df1<-df1%>%mutate(distance_per=distance/146)%>%mutate(value=round(distance_per,digit=2))%>%filter(chromosome=="17b")
df1$value[df1$distance_per>0.5]<-""
df1%>%select(value)
df1$AS_f = factor(df1$AS, levels=Xmono)
AS_labs<-c("W3_DXZ1_monomer1_X02418.1"="W3_m1","W4_DXZ1_monomer2_X02418.1"="W4_m2","W5_DXZ1_monomer3_X02418.1"="W5_m3","W1_DXZ1_monomer4_X02418.1"="W1_m4","W2_DXZ1_monomer5_X02418.1"="W2_m5","W3_DXZ1_monomer6_X02418.1"="W3_m6","W4_DXZ1_monomer7_X02418.1"="W4_m7","W5_DXZ1_monomer8_X02418.1"="W5_m8","W1_DXZ1_monomer9_X02418.1"="W1_m9","W2_DXZ1_monomer10_X02418.1"="W2_m10","W3_DXZ1_monomer11_X02418.1"="W3_m11","W4_DXZ1_monomer12_X02418.1"="W4_m12")
chr_lab<-c("1"="W4_m1", "2"="W5_m5", "3"="W1_m3", "4"="W2_m4", "5"="W3_m5", "6"="W2_m6", "7"="W3_m7", "8"="W4_m8", "9"="W5_m9","10"="W1_m10", "11"="W5_m11", "12"="W5_m12")
gp1<-df1%>%ggplot(.,aes(y=reorder(vs,-mono), x=AS, fill=distance_per))+geom_tile(color="white")+ geom_text(aes(label=sprintf("%s", value),size=0.01), colour="black", show.legend = F) +coord_fixed(10)+facet_grid(mono~AS_f, scales="free", space="free",labeller=labeller(AS_f=as_labeller(AS_labs), mono=as_labeller(chr_lab)), switch="both")+scale_fill_gradient(low="red",high = "white",na.value="white",breaks=c(0,0.1,0.2),limits=c(0,0.3))
gp1<-gp1+theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank())+theme(strip.text.y = element_text(angle=180))
print(gp1)



chr11<-c("W4", "W5", "W1", "W2", "W3")

chr17<-c("W1", "W2", "W3", "W4", "W5", "W1", "W2", "W3", "W4","W3", "W4", "W5", "W1", "W1", "W1", "W5")

chrX<-c("1"="W3_m1", "2"="W4_m5", "3"="W5_m3", "4"="W1_m4", "5"="W2_m5","6"="W3_m6", "7"="W4_m7", "8"="W5_m8", "9"="W1_m9", "10"="W2_m10","11"="W3_m11", "12"="W4_m12")

chr17<-c("1"="W4_m1", "2"="W5_m5", "3"="W1_m3", "4"="W2_m4", "5"="W3_m5", "6"="W2_m6", "7"="W3_m7", "8"="W4_m8", "9"="W5_m9","10"="W1_m10", "11"="W5_m11", "12"="W5_m12")
