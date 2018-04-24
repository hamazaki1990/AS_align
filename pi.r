library(tidyverse)
getwd()
setwd("~/github/AS_align/")

(ab<-read_csv("AB_grid_profile.csv"))
pAB<-ab%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="AB")
pAB

(a<-read_csv("A_grid_profile.csv"))
pA<-a%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="A")
pA


(b<-read_csv("B_grid_profile.csv"))
pB<-b%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="B")
pB


df<-pAB%>%bind_rows(pA,pB)
df$g_f = factor(df$group, levels=c("AB","A","B"))
gp<-df%>%ggplot(.,aes(x=bp,y=pi))+geom_line()+facet_grid(g_f~.)
print(gp)


(j1<-read_csv("J1_grid_profile.csv"))
pJ1<-j1%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)
pJ1

gp<-pJ1%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(j2<-read_csv("J2_grid_profile.csv"))
pJ2<-j2%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)
pJ2

gp<-pJ2%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(d1<-read_csv("D1_grid_profile.csv"))
pD1<-d1%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)
pD1

gp<-pD1%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(d1<-read_csv("D1_grid_profile.csv"))
pD1<-d1%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)
pD1

gp<-pD1%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(d2<-read_csv("D2_grid_profile.csv"))
pD2<-d2%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)
pD2

gp<-pD2%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(w1<-read_csv("W1_grid_profile.csv"))
pW1<-w1%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w1")
pW1

gp<-pW1%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(w2<-read_csv("W2_grid_profile.csv"))
pW2<-w2%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w2")
pW2

gp<-pW2%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(w3<-read_csv("W3_grid_profile.csv"))
pW3<-w3%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w3")
pW3

gp<-pW3%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(w4<-read_csv("W4_grid_profile.csv"))
pW4<-w4%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi) %>%mutate(group="w4")
pW4

gp<-pW4%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(w5<-read_csv("W5_grid_profile.csv"))
pW5<-w5%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w5")
pW5

gp<-pW5%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)


(w123<-read_csv("W1W2W3_grid_profile.csv"))
pW123<-w123%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w123")
pW123

gp<-pW123%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)

(w45<-read_csv("W4W5_grid_profile.csv"))
pW45<-w45%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w45")
pW45

gp<-pW45%>%ggplot(.,aes(x=bp,y=pi))+geom_line()
print(gp)





