library(ggplot2)
# Read in the files
prs <- pred_inf
hdl <- obj.bigsnp$fam$hdl
sex <- obj.bigsnp$fam$sex
CAD <- obj.bigsnp$fam$CAD
tg <- obj.bigsnp$fam$tg
ldl <- obj.bigsnp$fam$ldl
age <- obj.bigsnp$fam$age


dat2 <- data.frame(prs ,sex, CAD)
dat2["sex"][dat2["sex"] == 1] <- "M"
dat2["sex"][dat2["sex"] == 2] <- "F"
dat3 <- data.frame(dat2 , tg)





ggplot(dat2, aes(x=prs, y=cad, color=sex))+
  geom_point()+
  theme_classic()+
  labs(x="Polygenic Score", y="cad")




library(plyr)
mu <- ddply(dat2 , "sex", summarise, grp.mean=mean(prs))
head(mu)

ggplot(dat2 , aes(x=prs, color=sex)) +
  geom_density()
p<-ggplot(dat2, aes(x=prs, color=sex)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")
p

ggplot(dat2, aes(x=prs, y=tg, color=sex))+
  geom_point()+
  theme_classic()+
  labs(x="Polygenic Score", y="tg")



