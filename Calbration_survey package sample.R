library(survival)
library(survey)
library(dplyr)

## stratified on case status
str(nwtco)
dcchs <- twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel), subset=~I(in.subcohort | rel), data=nwtco)
svycoxph(Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12), design=dcchs)

## Using survival::cch 
subcoh <- nwtco$in.subcohort
selccoh <- with(nwtco, rel==1|subcoh==1)
coh.data <- nwtco[selccoh,]
coh.data$subcohort <- subcoh[selccoh]
cch(Surv(edrel, rel) ~ factor(stage) + factor(histol) + I(age/12), data =coh.data,
    subcoh = ~subcohort, id=~seqno, cohort.size=4028, method="LinYing")

## two-phase case-control
## Similar to Breslow & Chatterjee, Applied Statistics (1999) but with
## a slightly different version of the data set
nwtco$incc2 <- as.logical(with(nwtco, ifelse(rel | instit==2,1,rbinom(nrow(nwtco),1,.1))))
dccs2<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,instit)),
                data=nwtco, subset=~incc2)
dccs8<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,stage,instit)),
                data=nwtco, subset=~incc2)
summary(glm(rel~factor(stage)*factor(histol),data=nwtco,family=binomial()))
summary(svyglm(rel~factor(stage)*factor(histol),design=dccs2,family=quasibinomial()))
summary(svyglm(rel~factor(stage)*factor(histol),design=dccs8,family=quasibinomial()))

## Stratification on stage is really post-stratification, so we should use calibrate()
gccs8<-calibrate(dccs2, phase=2, formula=~interaction(rel,stage,instit))
summary(svyglm(rel~factor(stage)*factor(histol),design=gccs8,family=quasibinomial()))

## For this saturated model calibration is equivalent to estimating weights.
pccs8<-calibrate(dccs2, phase=2,formula=~interaction(rel,stage,instit), calfun="rrz")
summary(svyglm(rel~factor(stage)*factor(histol),design=pccs8,family=quasibinomial()))

## Since sampling is SRS at phase 1 and stratified RS at phase 2, we
## can use method="simple" to save memory.
dccs8_simple<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,stage,instit)),
                       data=nwtco, subset=~incc2,method="simple")
summary(svyglm(rel~factor(stage)*factor(histol),design=dccs8_simple,family=quasibinomial()))

## S3 method for class 'twophase'
calibrate(design=dcchs, phase=2, Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12), data=nwtco,
          calfun=c("linear","raking","logit","rrz"))
grake(mm,ww,calfun,eta=rep(0,NCOL(mm)),bounds,population,epsilon,
      verbose,maxit,variance=NULL)

#pubhealth.ku.dk
nwts <- read.csv("wilms2.csv")
head(nwts,5)
impmodel <- glm((histol-1)~instit+I(age>10)+I(stage==4)*study,
                data=nwts, subset=in.subsample, family=binomial)
nwts$imphist <- predict(impmodel, newdata=nwts, type="response")
nwts$imphist[nwts$in.subsample] <- nwts$histol[nwts$in.subsample]
ifmodel <- coxph(Surv(edrel,rel)~imphist*age+I(stage>2)*tumdiam, data=nwts)
inffun <- resid(ifmodel, "dfbeta")
colnames(inffun) <- paste("if",1:6,sep="")
nwts_if <- cbind(nwts, inffun)
if_design <- twophase(id = list(~1, ~1), subset = ~in.subsample,
                      strata = list(NULL, ~interaction(instit, rel)), data = nwts_if)
if_cal <- calibrate(if_design, phase=2, calfun="raking",
                    ~if1+if2+if3+if4+if5+if6+rel*instit)
m1 <- svycoxph(Surv(edrel,rel)~imphist*age+I(stage>2)*tumdiam, design=if_cal)
m2 <- svycoxph(Surv(edrel,rel)~histol*age+I(stage>2)*tumdiam, design=if_cal)
m3 <- coxph(Surv(edrel,rel)~imphist*age+I(stage>2)*tumdiam, data=nwts)
m4 <- coxph(Surv(edrel,rel)~histol*age+I(stage>2)*tumdiam, data=nwts)
