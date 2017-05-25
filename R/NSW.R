library(Matching)
data(lalonde)
dw<-lalonde

la_t<-read.table("nsw_treated.txt")
la_c<-read.table("nsw_control.txt")
colnames(la_t)<- c("treat", "age","educ","black","hisp", "married","nodegr","re75","re78")
colnames(la_c)<- c("treat", "age","educ","black","hisp", "married","nodegr","re75","re78")
lalonde<-rbind(la_t,la_c)

dw_t<-read.table("nswre74_treated.txt")
dw_c<-read.table("nswre74_control.txt")
colnames(dw_t)<- c("treat", "age","educ","black","hisp", "married","nodegr","re75","x","re78")
colnames(dw_c)<- c("treat", "age","educ","black","hisp", "married","nodegr","re75","x","re78")
dw_t<-dw_t[c("treat", "age","educ","black","hisp", "married","nodegr","re75","re78")]
dw_c<-dw_c[c("treat", "age","educ","black","hisp", "married","nodegr","re75","re78")]

dw<-rbind(dw,dw2)
dehejiawahha<-dw

cps<-read.table("cps_controls.txt")
colnames(cps)<- c("treat", "age","educ","black","hisp", "married","nodegr","re75","x","re78")
cps<-cps[c("treat", "age","educ","black","hisp", "married","nodegr","re75","re78")]

all<-rbind(la,dw)

all<- rbind(cps,dw_t)
all<-lalonde

fml<-as.formula("re78 ~ . -treat")
MH<- causalTree(fml,
                    all , treatment = all$treat,
                    split.Rule="user2", split.Honest=F, cv.option="CT", minsize = 50,
                    cv.alpha = 1, xval=0, cp=0,
                    propensity = rep(0.5, length(all$treat)))
FIT<- causalTree(fml,
                     all , treatment = all$treat,
                     split.Rule="fit", split.Honest=F, cv.option="CT", minsize = 50,
                     cv.alpha = 1, xval=0, cp=0,
                     propensity = rep(0.5, length(all$treat)))

CT<- causalTree(fml,
                    all , treatment = all$treat,
                    split.Rule="fit", split.Honest=T, cv.option="CT", minsize = 50,
                    cv.alpha = 1, xval=0, cp=0,
                    propensity = rep(0.5, length(all$treat)))

TOT<- causalTree(fml,
                all , treatment = all$treat,
                split.Rule="TOT", split.Honest=T, cv.option="CT", minsize = 50,
                cv.alpha = 1, xval=0, cp=0,
                propensity = rep(0.5, length(all$treat)))
TS <- causalTree(fml,
                all , treatment = all$treat,
                split.Rule="tstats", split.Honest=T, cv.option="CT", minsize = 50,
                cv.alpha = 1, xval=0, cp=0,
                propensity = rep(0.5, length(all$treat)))

MO <- causalTree(fml,
           all , treatment = all$treat,
           split.Rule="user5", split.Honest=T, cv.option="CT", minsize = 50,
           cv.alpha = 1, xval=0, cp=0,
           propensity = rep(0.5, length(all$treat)))

temp <- predict(brcaMH,newdata = lalonde)
temp2 <- predict(brcaMH, newdata = all, type = "vector")



