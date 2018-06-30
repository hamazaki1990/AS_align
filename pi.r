library(tidyverse)
getwd()
setwd("~/github/AS_align/")

(ab<-read_csv("AB_grid_profile.csv"))
pAB<-ab%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%mutate(conservation=1-pi)%>%select(bp,conservation)%>%mutate(group="AB")
pAB

(a<-read_csv("A_grid_profile.csv"))
pA<-a%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%mutate(conservation=1-pi)%>%select(bp,conservation)%>%mutate(group="A")
pA%>%filter(pi<0.1)


(b<-read_csv("B_grid_profile.csv"))
pB<-b%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%mutate(conservation=1-pi)%>%select(bp,conservation)%>%mutate(group="B")
pB%>%filter(pi<0.1)

(g_r<-tibble("xmin"=c(35,35,42,35,36,43.75,47),"xmax"=c(51,51,49,51,39,44.25,50),"ymin"=c(0,0,0,0,0,0,0),"ymax"=c(1,1,1,1,1,1,1),"group"=c("AB","A","A","B","B","B","B"), "fill"=c(gray,red,red,yellow,yellow,yellow,yellow)))
gp<-ggplot(g_r)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, alpha=0.1))+facet_wrap(~g_r$group)
print(gp)

df<-pAB
View(df)
df$g_f = factor(df$group, levels=c("AB","A","B"))
gp<-df%>%ggplot(.,aes(x=bp,y=conservation))+
  geom_segment(aes(x=35, y=-0.01,xend=51, yend=-0.01))+
  geom_line()+ggtitle("whole ASs")+theme_classic(base_size=11, base_family='')
print(gp)
ggsave("ABtype_pi.png",gp,width=5,height=3)

df<-pA
#gp<-df%>%ggplot(.,aes(x=bp,y=conservation))+geom_line()+geom_segment(aes(x=35, y=-0.01,xend=51, yend=-0.01))+geom_segment(aes(x=42, y=-0.005,xend=49, yend=-0.005))+ggtitle("Atype")+theme_classic(base_size=11, base_family='')
gp<-df%>%ggplot(.,aes(x=bp,y=conservation))+geom_line()+geom_rect(aes(xmin=35, ymin=-0,xmax=51, ymax=1), alpha=0.002, fill="#e41a1c")+geom_rect(aes(xmin=42, ymin=-0,xmax=49, ymax=1), alpha=0.002, fill="#e41a1c")+ggtitle("Atype")+theme_classic(base_size=11, base_family='')
print(gp)
ggsave("Atype_pi.png",gp,width=5,height=3)

df<-pB
gp<-df%>%ggplot(.,aes(x=bp,y=conservation))+geom_line()+geom_segment(aes(x=35, y=-0.01,xend=51, yend=-0.01))+geom_segment(aes(x=36, y=-0.005,xend=39, yend=-0.005))+geom_segment(aes(x=43.75, y=-0.005,xend=44.25, yend=-0.005))+geom_segment(aes(x=47, y=-0.005,xend=50, yend=-0.005))+ggtitle("Btype")+theme_classic(base_size=11, base_family='')
print(gp)
ggsave("Btype_pi.png",gp,width=5,height=3)


(df<-pAB%>%bind_rows(pA,pB,g_r))
df$g_f = factor(df$group, levels=c("AB","A","B"))
g_r$g_f = factor(g_r$group, levels=c("AB","A","B"))
gp<-ggplot(df)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill,alpha=0.1))+scale_fill_manual(values=c("#cccccc", "#4daf4a", "orange"))+geom_line(aes(x=bp,y=conservation))+facet_grid(g_f~.)+theme_classic(base_size=11, base_family='')+ theme(legend.position = 'none')
print(gp)
ggsave("pi.png",gp,width=5,height=5)

(g_r2<-tibble("xmin"=c(35,42,35,36,43.75,47),"xmax"=c(51,49,51,39,44.25,50),"ymin"=c(0.3,0.3,0.3,0.3,0.3,0.3),"ymax"=c(1,1,1,1,1,1),"group"=c("A","A","B","B","B","B"), "fill"=c("#4daf4a","#4daf4a","orange","orange","orange","orange")))

(df<-pA%>%bind_rows(pB,g_r2))
df$g_f = factor(df$group, levels=c("AB","A","B"))
g_r$g_f = factor(g_r$group, levels=c("AB","A","B"))
gp<-ggplot(df)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill,alpha=0.1))+scale_fill_manual(values=c( "#4daf4a", "#ff7f00"))+geom_line(aes(x=bp,y=conservation))+facet_grid(g_f~.)+theme_classic(base_size=11, base_family='')+ theme(legend.position = 'none')
print(gp)
ggsave("pi2.png",gp,width=5,height=4)




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

df<-pW1%>%bind_rows(pW2,pW3)
gp<-df%>%ggplot(.,aes(x=bp,y=pi, color=group))+geom_line()+facet_grid(group~.)
print(gp)


(w45<-read_csv("W4W5_grid_profile.csv"))
pW45<-w45%>%mutate(bp=c(1:167),ASs=`-`+A+C+G+T)%>%mutate(pi=1-((`-`)*(`-`-1)+A*(A-1)+C*(C-1)+G*(G-1)+T*(T-1))/(ASs*(ASs-1)))%>%select(bp,pi)%>%mutate(group="w45")
pW45

df<-pW4%>%bind_rows(pW5)
gp<-df%>%ggplot(.,aes(x=bp,y=pi, color=group))+geom_line()+facet_grid(group~.)
print(gp)





