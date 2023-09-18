rm(list=ls())
library(openxlsx)
library(dplyr)
sc1 <- read.xlsx("SC1.xlsx")
sc1 <- sc1 %>%
  mutate(sdra = as.factor(sdra),
         immunodep = as.factor(immunodep),
         centre = as.factor(centre),
         lactat = as.numeric(lactat),
         dcd.last.FU = as.factor(dcd.last.FU),
         bras_rando_1 = as.factor(bras_rando_1))

sc1_lactat <- sc1_lactat[which(sc1_lactat$etude == "SC"),]

sc1_lactat <- sc1%>%
  filter(lactat > 2)


sc1_lactat %>%
  count(bras_rando_1)
sc1_lactat %>%
  count(bras_rando_1,dcd.last.FU)

#33/53 no cooling and death
#23/51 cooling and death
pA = 23/51
pB  = 33/53
OR_unadjusted = (pA/(1-pA))/(pB/(1-pB))


sc1_lactat$centre <- droplevels(sc1_lactat$centre)
sc1_lactat$centre = factor(sc1_lactat$centre, levels = c("5","1","3","4","7","8","9")) #Set center 5 as the reference class 


res <- glm(data = sc1_lactat, family = "binomial", formula = dcd.last.FU ~ centre + age + immunodep + sdra + vp + lactat + PF + bras_rando_1 - 1)
coeff <-  res$coefficients
beta_adjusted <- coeff["bras_rando_11"] 
std <- summary(res)$coefficients["bras_rando_11",2]
save(coeff, file = "beta.rda")
OR_adjusted = exp(beta_adjusted) #Adjusted Odd ratio


pA=pA
pB=pB
kappa=1
alpha=0.05
beta=0.20
(OR=pA*(1-pB)/pB/(1-pA)) # 2
(nB=(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))*((qnorm(1-alpha/2)+qnorm(1-beta))/log(OR))^2)
ceiling(nB) # 156
z=log(OR)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))
(Power=pnorm(z-qnorm(1-alpha/2))+pnorm(-z-qnorm(1-alpha/2))) #270 patients necessary for a power of 80%, alpha = 0.05, unadjusted OR



pA=pA
pB=pB
kappa=1
alpha=0.05
beta=0.20
OR=OR_adjusted
(nB=(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))*((qnorm(1-alpha/2)+qnorm(1-beta))/log(OR))^2)
ceiling(nB) # 156
z=log(OR)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))
(Power=pnorm(z-qnorm(1-alpha/2))+pnorm(-z-qnorm(1-alpha/2))) #14548 patients necessary for a power of 80%, alpha = 0.05, adjusted OR


