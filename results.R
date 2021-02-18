rm(list=ls())
data=read.csv("H:/covid-19/Design 2 Qi .csv")
dim(data)
require(survival)
require(survminer)
install.packages("survminer")
install.packages("dplyr")
install.packages("frailtyHL")
require(survival)
require(survminer)
require(dplyr)
require(frailtyHL)
cmuc_data=data[data$site==0,]
mgh_data=data[data$site==1,]
wcmc_data=data[data$site==2,]
data$three_cat[data$site==1]=data$three_cat[data$site==1]-1
data$outcome_1[data$outcome_1==3]=2
table(data$outcome_1[data$site==0])
table(data$recovered[data$site==0])
data$recovered[setdiff(which(data$recovered[data$site==0]==1),which(data$outcome_1[data$site==0]==2))]=0
data$outcome_1[data$outcome_1 %in% c(1,0)]=1-data$outcome_1[data$outcome_1 %in% c(1,0)]
##cox model
table(data$day_fiftyfive_roc_cat)
data$outcome_1=as.factor(data$outcome_1)
data$day_fiftyfive_roc_cat=as.factor(data$day_fiftyfive_roc_cat)
data$day_sev_roc_cat=as.factor(data$day_sev_roc_cat)
data$SEX=as.factor(data$SEX)
data$RACE_ETHNIC=as.factor(data$RACE_ETHNIC)
data$obese=as.factor(data$obese)
data$berlin_intubation=as.factor(data$berlin_intubation)
data$hgb_7=as.factor(data$hgb_7)
data$CVVH=as.factor(data$CVVH)
data$three_cat=as.factor(data$three_cat)
data$berlin_ever=as.factor(data$berlin_ever)
data$ever55=as.factor(data$ever55)
##model 1
survival_data=data[c("berlin_ever","three_cat","CVVH","hgb_7","berlin_intubation",
                     "obese","RACE_ETHNIC","SEX","recovered","outcome_1",
                     "day_fiftyfive_roc_cat","day_sev_roc_cat","drips","midaz4","days_paralytics",
                     "keta4","fent2","outcome_time","person_id","site")]
sfit1 <- survfit(Surv(time=outcome_time,event=recovered) ~ day_fiftyfive_roc_cat, data=survival_data)
ggsurvplot(sfit1,fun="event",conf.int=TRUE, pval=TRUE, risk.table=TRUE,
           legend.labs=c("ever55>0", "ever55>1","ever55>2"), legend.title="Hypoxemia",
           title="Kaplan-Meier Curve",
           risk.table.height=.15)
rfit1 <- coxph(Surv(time=outcome_time,event=recovered) ~ day_fiftyfive_roc_cat +frailty.gamma(site), data=data)
summary(rfit1)
rfit2 <- coxph(Surv(time=outcome_time,event=recovered) ~ day_fiftyfive_roc_cat+ age +frailty.gamma(site), data=data)
summary(rfit2)
rfit3 <- coxph(Surv(time=outcome_time,event=recovered) ~ day_fiftyfive_roc_cat+ SEX +frailty.gamma(site), data=data)
summary(rfit3)
rfit4 <- coxph(Surv(time=outcome_time,event=recovered) ~ day_fiftyfive_roc_cat+ RACE_ETHNIC +frailty.gamma(site), data=data)
rfit5 <- coxph(Surv(time=outcome_time,event=recovered) ~ day_fiftyfive_roc_cat+ berlin_intubation +frailty.gamma(site), data=data)
summary(rfit4)
survival_data=data[c("berlin_ever","three_cat","CVVH","hgb_7","berlin_intubation",
                     "obese","RACE_ETHNIC","SEX","recovered","outcome_1",
                     "day_fiftyfive_roc_cat","day_sev_roc_cat","drips","midaz4","days_paralytics",
                     "keta4","fent2","outcome_time","person_id","site","age")]

## model 1 and 2


data$surtime=data$outcome_time

data$status=varhandle::unfactor(data$outcome_1)

data$center=data$site

data$CHEMO=as.factor(data$ever55)

data_conti=data_surv=data

jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")

res1<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 33
# beta_h     se_beh   t_value    p_value
# -0.8843871 0.09244868 -9.566249 0.00000000
# -0.4452227 0.20403115 -2.182131 0.02909987
# alpha_h    rho_h
# [1,] 0.08533261 2.248363
# h0        hp     p.bvh
# [1,] -3835.808 -3835.839 -3842.408
# cAIC     rAIC
# [1,] 7679.479 7688.815
p.res1=res1$F.Est[1:(dim(res1$F.Est)[1]/2),]
p.res1=matrix(p.res1,nrow = 1)
rownames(p.res1)="ever55"
colnames(p.res1)=colnames(res1$F.Est)

print(p.res1)

#           beta_h     se_beh   t_value p_value
# ever55 -0.8843871 0.09244868 -9.566249       0

jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")

res2<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 54
# beta_h      se_beh     t_value      p_value
# -0.58751894 0.102681854  -5.7217407 1.054382e-08
# -0.50314288 0.060603285  -8.3022378 0.000000e+00
# -0.20753715 0.065556650  -3.1657681 1.546740e-03
# -0.15813608 0.051484179  -3.0715471 2.129526e-03
# 0.05344306 0.039253053   1.3615006 1.733555e-01
# -0.05147004 0.011853064  -4.3423402 1.409730e-05
# 0.77002662 0.251228553   3.0650442 2.176378e-03
# 0.74951281 0.230728292   3.2484652 1.160294e-03
# 0.89754717 0.237058362   3.7861865 1.529769e-04
# -0.33280617 0.101786829  -3.2696389 1.076849e-03
# -0.03423672 0.003390513 -10.0977992 0.000000e+00
# -0.28869223 0.091982285  -3.1385633 1.697782e-03
# 0.03685906 0.162767949   0.2264516 8.208502e-01
# 0.68455903 0.311156758   2.2000455 2.780367e-02
# 0.32623236 0.129347174   2.5221452 1.166415e-02
# 0.13904947 0.199009956   0.6987061 4.847357e-01
# -0.11660224 0.169579712  -0.6875955 4.917076e-01
# 0.17368488 0.146014184   1.1895069 2.342402e-01
# -0.30146505 0.216384142  -1.3931938 1.635612e-01
# -0.67217979 0.129965534  -5.1719850 2.316202e-07
# -0.63577924 0.155293421  -4.0940513 4.239004e-05
# 0.03666203 0.106416213   0.3445155 7.304587e-01
# 0.10348007 0.084621497   1.2228580 2.213833e-01
# 0.07484110 0.018744534   3.9926892 6.532816e-05
# 0.41934943 0.685353135   0.6118735 5.406214e-01
# 0.58515149 0.620831884   0.9425281 3.459223e-01
# 0.71313075 0.628550963   1.1345631 2.565584e-01
# 0.10934402 0.205438113   0.5322480 5.945543e-01
# 0.02103978 0.008579779   2.4522515 1.419654e-02
# 0.36899412 0.212457711   1.7367886 8.242450e-02
# 0.29465895 0.328858454   0.8960054 3.702499e-01
# 0.33048900 0.747978105   0.4418431 6.586027e-01
# 0.39567239 0.287741748   1.3750955 1.691018e-01
# 0.84956337 0.409827871   2.0729761 3.817451e-02
# 0.42995380 0.305410645   1.4077892 1.591935e-01
# 0.08886044 0.338367665   0.2626150 7.928473e-01
# alpha_h    rho_h
# [1,] 0.02309139 1.150885
# h0       hp     p.bvh
# [1,] -3529.755 -3527.62 -3578.205
# cAIC    rAIC
# [1,] 7134.553 7160.41

p.res2=res2$F.Est[1:(dim(res2$F.Est)[1]/2),]

p.res2=res2$F.Est[1:(dim(res2$F.Est)[1]/2),]
rownames(p.res2)= c("ever551","drips","midaz4", "fent2", "keta4",
                    "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                    "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res2)

#                      beta_h      se_beh     t_value      p_value
# ever551            -0.58751894 0.102681854  -5.7217407 1.054382e-08
# drips              -0.50314288 0.060603285  -8.3022378 0.000000e+00
# midaz4             -0.20753715 0.065556650  -3.1657681 1.546740e-03
# fent2              -0.15813608 0.051484179  -3.0715471 2.129526e-03
# keta4               0.05344306 0.039253053   1.3615006 1.733555e-01
# days_paralytics    -0.05147004 0.011853064  -4.3423402 1.409730e-05
# berlin_intubation1  0.77002662 0.251228553   3.0650442 2.176378e-03
# berlin_intubation2  0.74951281 0.230728292   3.2484652 1.160294e-03
# berlin_intubation3  0.89754717 0.237058362   3.7861865 1.529769e-04
# CVVH1              -0.33280617 0.101786829  -3.2696389 1.076849e-03
# age                -0.03423672 0.003390513 -10.0977992 0.000000e+00
# SEX1               -0.28869223 0.091982285  -3.1385633 1.697782e-03
# RACE_ETHNIC1        0.03685906 0.162767949   0.2264516 8.208502e-01
# RACE_ETHNIC2        0.68455903 0.311156758   2.2000455 2.780367e-02
# RACE_ETHNIC3        0.32623236 0.129347174   2.5221452 1.166415e-02
# RACE_ETHNIC4        0.13904947 0.199009956   0.6987061 4.847357e-01
# RACE_ETHNIC5       -0.11660224 0.169579712  -0.6875955 4.917076e-01
# RACE_ETHNIC6        0.17368488 0.146014184   1.1895069 2.342402e-01

##model 3 and 4

data$CHEMO=as.factor(data$day_fiftyfive_roc_cat)
data_conti=data_surv=data
jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")
res3<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 31
# beta_h    se_beh    t_value      p_value
# -0.3951719 0.1107883  -3.566909 0.0003612163
# -1.2747329 0.1131431 -11.266559 0.0000000000
# -0.3540133 0.2661624  -1.330065 0.1834967998
# -0.6290633 0.2221252  -2.832021 0.0046254847
# alpha_h    rho_h
# [1,] 0.1598773 1.672013
# h0        hp     p.bvh
# [1,] -3813.245 -3814.226 -3821.769
# cAIC     rAIC
# [1,] 7638.389 7647.537

p.res3=res3$F.Est[1:(dim(res3$F.Est)[1]/2),]

rownames(p.res3)=c("day_fiftyfive_roc_cat1","day_fiftyfive_roc_cat2")
print(p.res3)
#                      beta_h    se_beh    t_value      p_value
# day_fiftyfive_roc_cat1 -0.3951719 0.1107883  -3.566909 0.0003612163
# day_fiftyfive_roc_cat2 -1.2747329 0.1131431 -11.266559 0.0000000000


jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
res4<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 40
# beta_h      se_beh    t_value      p_value
# -0.31614300 0.120011635 -2.6342695 8.431852e-03
# -0.87240892 0.125588276 -6.9465793 3.742562e-12
# -0.50823789 0.061197122 -8.3049313 0.000000e+00
# -0.20169168 0.065371215 -3.0853286 2.033274e-03
# -0.14474848 0.052347540 -2.7651439 5.689769e-03
# 0.06266637 0.040531522  1.5461145 1.220769e-01
# -0.04462275 0.012033851 -3.7081024 2.088182e-04
# 0.79822485 0.251851511  3.1694265 1.527401e-03
# 0.73901790 0.231008644  3.1990920 1.378612e-03
# 0.93102611 0.237839862  3.9145083 9.058863e-05
# -0.33242705 0.101957710 -3.2604405 1.112393e-03
# -0.03275962 0.003386443 -9.6737548 0.000000e+00
# -0.31672374 0.092373333 -3.4287356 6.063999e-04
# 0.02240294 0.163414810  0.1370925 8.909577e-01
# 0.63187936 0.312567884  2.0215748 4.322030e-02
# 0.30257715 0.129689731  2.3330849 1.964369e-02
# 0.09398965 0.200575780  0.4685992 6.393562e-01
# -0.23761568 0.172975678 -1.3736942 1.695366e-01
# 0.14169862 0.146836809  0.9650075 3.345411e-01
# -0.20992574 0.282678747 -0.7426301 4.577057e-01
# -0.37312938 0.238164255 -1.5666893 1.171873e-01
# -0.67781347 0.131137971 -5.1687048 2.357218e-07
# -0.63641074 0.155039105 -4.1048401 4.045946e-05
# 0.04045426 0.106630561  0.3793871 7.044004e-01
# 0.11201596 0.084219868  1.3300420 1.835044e-01
# 0.07594847 0.018982853  4.0008985 6.310243e-05
# 0.42156893 0.685514282  0.6149674 5.385763e-01
# 0.57931667 0.620740888  0.9332665 3.506824e-01
# 0.71445557 0.628878772  1.1360784 2.559238e-01
# 0.11976215 0.205492089  0.5828066 5.600235e-01
# 0.02156143 0.008622916  2.5004805 1.240250e-02
# 0.37050810 0.212482178  1.7437138 8.120901e-02
# 0.31047773 0.328433366  0.9453294 3.444907e-01
# 0.32207698 0.747393852  0.4309334 6.665168e-01
# 0.38503088 0.288495801  1.3346152 1.820023e-01
# 0.85040304 0.409444302  2.0769688 3.780444e-02
# 0.40144345 0.308872575  1.2997057 1.937018e-01
# 0.08533646 0.338272132  0.2522716 8.008311e-01
# alpha_h     rho_h
# [1,] 0.05181345 0.8425671
# h0        hp     p.bvh
# [1,] -3520.979 -3520.156 -3571.518
# cAIC     rAIC
# [1,] 7121.399 7147.036

p.res4=res4$F.Est[1:(dim(res4$F.Est)[1]/2),]

rownames(p.res4)= c("day_fiftyfive_roc_cat1","day_fiftyfive_roc_cat2","drips","midaz4", "fent2", "keta4",
                    "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                    "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res4)
#                             beta_h      se_beh    t_value      p_value
# day_fiftyfive_roc_cat1 -0.31614300 0.120011635 -2.6342695 8.431852e-03
# day_fiftyfive_roc_cat2 -0.87240892 0.125588276 -6.9465793 3.742562e-12
# drips                  -0.50823789 0.061197122 -8.3049313 0.000000e+00
# midaz4                 -0.20169168 0.065371215 -3.0853286 2.033274e-03
# fent2                  -0.14474848 0.052347540 -2.7651439 5.689769e-03
# keta4                   0.06266637 0.040531522  1.5461145 1.220769e-01
# days_paralytics        -0.04462275 0.012033851 -3.7081024 2.088182e-04
# berlin_intubation1      0.79822485 0.251851511  3.1694265 1.527401e-03
# berlin_intubation2      0.73901790 0.231008644  3.1990920 1.378612e-03
# berlin_intubation3      0.93102611 0.237839862  3.9145083 9.058863e-05
# CVVH1                  -0.33242705 0.101957710 -3.2604405 1.112393e-03
# age                    -0.03275962 0.003386443 -9.6737548 0.000000e+00
# SEX1                   -0.31672374 0.092373333 -3.4287356 6.063999e-04
# RACE_ETHNIC1            0.02240294 0.163414810  0.1370925 8.909577e-01
# RACE_ETHNIC2            0.63187936 0.312567884  2.0215748 4.322030e-02
# RACE_ETHNIC3            0.30257715 0.129689731  2.3330849 1.964369e-02
# RACE_ETHNIC4            0.09398965 0.200575780  0.4685992 6.393562e-01
# RACE_ETHNIC5           -0.23761568 0.172975678 -1.3736942 1.695366e-01
# RACE_ETHNIC6            0.14169862 0.146836809  0.9650075 3.345411e-01

### model 5 and 6

data$CHEMO=data$day_fiftyfive_roc
data_conti=data_surv=data
jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO +(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")
res5<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 22
# beta_h     se_beh    t_value    p_value
# -0.17361193 0.01658933 -10.465279 0.00000000
# -0.04938983 0.02195893  -2.249191 0.02450033
# alpha_h  rho_h
# [1,] 0.2440543 1.3631
# h0       hp     p.bvh
# [1,] -3789.446 -3791.07 -3800.658
# cAIC     rAIC
# [1,] 7586.821 7605.316


p.res5=res5$F.Est[1:(dim(res5$F.Est)[1]/2),]
p.res5=matrix(p.res5,nrow = 1)
rownames(p.res5)="day_fiftyfive_roc"
colnames(p.res5)=colnames(res5$F.Est)
print(p.res5)

#                   beta_h     se_beh   t_value p_value
# day_fiftyfive_roc -0.1736119 0.01658933 -10.46528       0

jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
res6<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 30
# beta_h      se_beh    t_value      p_value
# -0.13067375 0.017978862 -7.2681878 3.643752e-13
# -0.51416164 0.063178095 -8.1382897 4.440892e-16
# -0.21942850 0.066710894 -3.2892453 1.004564e-03
# -0.10160533 0.055597868 -1.8275041 6.762401e-02
# 0.06667344 0.041017597  1.6254837 1.040595e-01
# -0.02351322 0.012949867 -1.8157113 6.941466e-02
# 0.83963768 0.254102376  3.3043283 9.520426e-04
# 0.76615850 0.231588630  3.3082734 9.387312e-04
# 0.87082660 0.236908298  3.6757961 2.371089e-04
# -0.33270642 0.103065806 -3.2280970 1.246167e-03
# -0.03233615 0.003388971 -9.5415835 0.000000e+00
# -0.37915710 0.093206539 -4.0679238 4.743389e-05
# -0.06484209 0.163396376 -0.3968392 6.914860e-01
# 0.56498478 0.311202909  1.8154868 6.944913e-02
# 0.25653372 0.128360369  1.9985430 4.565782e-02
# -0.01914359 0.197586389 -0.0968872 9.228160e-01
# -0.29838764 0.169520975 -1.7601813 7.837708e-02
# 0.11747751 0.145155480  0.8093219 4.183300e-01
# -0.01855695 0.025583561 -0.7253466 4.682394e-01
# -0.68142177 0.130853363 -5.2075220 1.913792e-07
# -0.65564010 0.155203604 -4.2243871 2.395919e-05
# 0.05456189 0.107296304  0.5085160 6.110915e-01
# 0.11560925 0.081981925  1.4101797 1.584866e-01
# 0.07845734 0.020779401  3.7757265 1.595420e-04
# 0.39075469 0.684854768  0.5705658 5.682940e-01
# 0.52068573 0.618439811  0.8419344 3.998247e-01
# 0.62686496 0.623167787  1.0059329 3.144479e-01
# 0.10101830 0.205428863  0.4917435 6.229007e-01
# 0.02033531 0.008588537  2.3677273 1.789772e-02
# 0.37234332 0.215138372  1.7307155 8.350251e-02
# 0.22979895 0.325434531  0.7061296 4.801076e-01
# 0.33673076 0.745441858  0.4517197 6.514709e-01
# 0.34803292 0.285127680  1.2206213 2.222294e-01
# 0.78746004 0.406599906  1.9366951 5.278263e-02
# 0.37232937 0.301484200  1.2349880 2.168350e-01
# 0.04736125 0.336384589  0.1407949 8.880319e-01
# alpha_h     rho_h
# [1,] 0.06941902 0.5565561
# h0        hp     p.bvh
# [1,] -3508.282 -3507.921 -3561.304
# cAIC     rAIC
# [1,] 7092.098 7126.608
p.res6=res6$F.Est[1:(dim(res6$F.Est)[1]/2),]
rownames(p.res6)= c("day_fiftyfive_roc","drips","midaz4", "fent2", "keta4",
                    "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                    "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res6)

#                    beta_h      se_beh    t_value      p_value
# day_fiftyfive_roc  -0.13067375 0.017978862 -7.2681878 3.643752e-13
# drips              -0.51416164 0.063178095 -8.1382897 4.440892e-16
# midaz4             -0.21942850 0.066710894 -3.2892453 1.004564e-03
# fent2              -0.10160533 0.055597868 -1.8275041 6.762401e-02
# keta4               0.06667344 0.041017597  1.6254837 1.040595e-01
# days_paralytics    -0.02351322 0.012949867 -1.8157113 6.941466e-02
# berlin_intubation1  0.83963768 0.254102376  3.3043283 9.520426e-04
# berlin_intubation2  0.76615850 0.231588630  3.3082734 9.387312e-04
# berlin_intubation3  0.87082660 0.236908298  3.6757961 2.371089e-04
# CVVH1              -0.33270642 0.103065806 -3.2280970 1.246167e-03
# age                -0.03233615 0.003388971 -9.5415835 0.000000e+00
# SEX1               -0.37915710 0.093206539 -4.0679238 4.743389e-05
# RACE_ETHNIC1       -0.06484209 0.163396376 -0.3968392 6.914860e-01
# RACE_ETHNIC2        0.56498478 0.311202909  1.8154868 6.944913e-02
# RACE_ETHNIC3        0.25653372 0.128360369  1.9985430 4.565782e-02
# RACE_ETHNIC4       -0.01914359 0.197586389 -0.0968872 9.228160e-01
# RACE_ETHNIC5       -0.29838764 0.169520975 -1.7601813 7.837708e-02
# RACE_ETHNIC6        0.11747751 0.145155480  0.8093219 4.183300e-01

##model 7 and 8

data$CHEMO=as.factor(data$ever70)
data_conti=data_surv=data
data$CHEMO
jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO +(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")
res7<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=200)

# [1] "iterations : "
# [1] 43
# beta_h    se_beh   t_value      p_value
# -1.117425 0.1581829 -7.064135 1.616263e-12
# -1.226608 0.3222225 -3.806710 1.408277e-04
# alpha_h    rho_h
# [1,] 0.07372755 2.414308
# h0       hp     p.bvh
# [1,] -3857.768 -3857.58 -3863.309
# cAIC     rAIC
# [1,] 7723.403 7730.619
p.res7=res7$F.Est[1:(dim(res7$F.Est)[1]/2),]
p.res7=res7$F.Est[1:(dim(res7$F.Est)[1]/2),]
p.res7=matrix(p.res7,nrow = 1)
rownames(p.res7)="ever70"
colnames(p.res7)=colnames(res7$F.Est)
print(p.res7)

#       beta_h    se_beh   t_value      p_value
# ever70 -1.117425 0.1581829 -7.064135 1.616263e-12

jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
res8<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=400)
# [1] "iterations : "
# [1] 400
# beta_h      se_beh      t_value      p_value
# 0.012890825 0.180782080   0.07130588 9.431543e-01
# -0.515294624 0.061583262  -8.36744606 0.000000e+00
# -0.222968224 0.064957700  -3.43251415 5.980127e-04
# -0.173723258 0.050507913  -3.43952553 5.827348e-04
# 0.060721735 0.035943843   1.68935011 9.115236e-02
# -0.053251143 0.010959927  -4.85871314 1.181512e-06
# 0.732193871 0.252202688   2.90319614 3.693752e-03
# 0.638946763 0.233950453   2.73112001 6.311948e-03
# 0.726679315 0.236840646   3.06822046 2.153377e-03
# -0.348141787 0.101221270  -3.43941336 5.829763e-04
# -0.036293048 0.003313978 -10.95150385 0.000000e+00
# -0.291767143 0.092293512  -3.16129635 1.570686e-03
# -0.096493516 0.159452618  -0.60515480 5.450761e-01
# 0.582655518 0.307677891   1.89371917 5.826230e-02
# 0.247788334 0.126869919   1.95308972 5.080897e-02
# -0.007145448 0.195720847  -0.03650837 9.708770e-01
# -0.263833592 0.169785107  -1.55392658 1.202019e-01
# 0.063306988 0.144595000   0.43782280 6.615147e-01
# -0.392724671 0.373180034  -1.05237321 2.926284e-01
# -0.667517859 0.131611956  -5.07186338 3.939392e-07
# -0.631958808 0.156224541  -4.04519549 5.227948e-05
# 0.042027307 0.107110791   0.39237230 6.947832e-01
# 0.092392074 0.086385924   1.06952695 2.848323e-01
# 0.074960688 0.018796397   3.98803487 6.662287e-05
# 0.584346683 0.706469852   0.82713605 4.081600e-01
# 0.770969603 0.655686035   1.17582129 2.396663e-01
# 0.896902379 0.664289187   1.35016857 1.769619e-01
# 0.089724039 0.206661194   0.43416007 6.641722e-01
# 0.020534807 0.008509425   2.41318387 1.581384e-02
# 0.399376591 0.212479605   1.87959965 6.016266e-02
# 0.227906664 0.329584647   0.69149660 4.892535e-01
# 0.250625067 0.755340697   0.33180400 7.400373e-01
# 0.363952422 0.286733231   1.26930674 2.043317e-01
# 0.784013843 0.407673721   1.92314050 5.446241e-02
# 0.394193347 0.301982243   1.30535274 1.917728e-01
# 0.057391562 0.337629084   0.16998406 8.650227e-01
# alpha_h    rho_h
# [1,] 1.719244e-06 118.0922
# h0        hp   p.bvh
# [1,] -3546.825 -3530.043 -3592.8
# cAIC   rAIC
# [1,] 7167.137 7189.6

p.res8=res8$F.Est[1:(dim(res8$F.Est)[1]/2),]
rownames(p.res8)= c("ever70","drips","midaz4", "fent2", "keta4",
                    "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                    "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res8)
#                      beta_h      se_beh      t_value      p_value
# ever70              0.012890825 0.180782080   0.07130588 9.431543e-01
# drips              -0.515294624 0.061583262  -8.36744606 0.000000e+00
# midaz4             -0.222968224 0.064957700  -3.43251415 5.980127e-04
# fent2              -0.173723258 0.050507913  -3.43952553 5.827348e-04
# keta4               0.060721735 0.035943843   1.68935011 9.115236e-02
# days_paralytics    -0.053251143 0.010959927  -4.85871314 1.181512e-06
# berlin_intubation1  0.732193871 0.252202688   2.90319614 3.693752e-03
# berlin_intubation2  0.638946763 0.233950453   2.73112001 6.311948e-03
# berlin_intubation3  0.726679315 0.236840646   3.06822046 2.153377e-03
# CVVH1              -0.348141787 0.101221270  -3.43941336 5.829763e-04
# age                -0.036293048 0.003313978 -10.95150385 0.000000e+00
# SEX1               -0.291767143 0.092293512  -3.16129635 1.570686e-03
# RACE_ETHNIC1       -0.096493516 0.159452618  -0.60515480 5.450761e-01
# RACE_ETHNIC2        0.582655518 0.307677891   1.89371917 5.826230e-02
# RACE_ETHNIC3        0.247788334 0.126869919   1.95308972 5.080897e-02
# RACE_ETHNIC4       -0.007145448 0.195720847  -0.03650837 9.708770e-01
# RACE_ETHNIC5       -0.263833592 0.169785107  -1.55392658 1.202019e-01
# RACE_ETHNIC6        0.063306988 0.144595000   0.43782280 6.615147e-01


##model 9 and 10

data$CHEMO=as.factor(data$day_sev_roc_cat)
data_conti=data_surv=data

jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO +(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")
res9<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 39
# beta_h    se_beh   t_value      p_value
# 0.4230019 0.2212009  1.912298 5.583804e-02
# -1.2214072 0.1589824 -7.682656 1.554312e-14
# -0.5747504 0.5413369 -1.061724 2.883610e-01
# -1.2967861 0.3241654 -4.000384 6.323984e-05
# alpha_h    rho_h
# [1,] 0.07639214 2.443784
# h0        hp     p.bvh
# [1,] -3826.426 -3826.293 -3832.682
# cAIC     rAIC
# [1,] 7664.724 7669.364

p.res9=res9$F.Est[1:(dim(res9$F.Est)[1]/2),]

rownames(p.res9)=c("day_sev_roc_cat1","day_sev_roc_cat2")

print(p.res9)

#                   beta_h    se_beh   t_value      p_value
# day_sev_roc_cat1  0.4230019 0.2212009  1.912298 5.583804e-02
# day_sev_roc_cat2 -1.2214072 0.1589824 -7.682656 1.554312e-14


jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")


res10<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=200)
# [1] "iterations : "
# [1] 143
# beta_h      se_beh       t_value      p_value
# 0.760908249 0.244224637   3.115608070 1.835661e-03
# -0.104316214 0.184260162  -0.566135475 5.713017e-01
# -0.495249627 0.062016991  -7.985708764 1.332268e-15
# -0.243508982 0.065268281  -3.730893136 1.908021e-04
# -0.136467917 0.050289569  -2.713642620 6.654792e-03
# 0.046945935 0.036890996   1.272558057 2.031749e-01
# -0.052932280 0.011178374  -4.735239511 2.187966e-06
# 0.713459159 0.253155632   2.818263039 4.828424e-03
# 0.631192870 0.234706602   2.689284686 7.160532e-03
# 0.736850441 0.238288867   3.092257097 1.986407e-03
# -0.369638087 0.101406176  -3.645123997 2.672629e-04
# -0.035208007 0.003358233 -10.484086787 0.000000e+00
# -0.233814719 0.093079142  -2.511999077 1.200494e-02
# -0.066394571 0.160270710  -0.414265155 6.786799e-01
# 0.632512111 0.308562662   2.049866004 4.037751e-02
# 0.203418168 0.128708622   1.580454868 1.140027e-01
# -0.032860036 0.196989266  -0.166811305 8.675185e-01
# -0.284549820 0.169850913  -1.675291673 9.387702e-02
# 0.026567332 0.145788386   0.182232157 8.554005e-01
# 0.001787614 0.574182738   0.003113318 9.975159e-01
# -0.471428454 0.377432347  -1.249040942 2.116501e-01
# -0.657007878 0.131701470  -4.988614625 6.081382e-07
# -0.643496599 0.156609701  -4.108919153 3.975152e-05
# 0.053093385 0.107106817   0.495704997 6.201026e-01
# 0.099677087 0.087002044   1.145686718 2.519248e-01
# 0.073032899 0.018934029   3.857229680 1.146794e-04
# 0.605310901 0.707593836   0.855449652 3.923023e-01
# 0.782831900 0.657102875   1.191338418 2.335208e-01
# 0.900435294 0.666284735   1.351427171 1.765586e-01
# 0.103262859 0.206518221   0.500018151 6.170623e-01
# 0.020898513 0.008528586   2.450407784 1.426945e-02
# 0.411473713 0.212252926   1.938600896 5.254995e-02
# 0.263698806 0.330360734   0.798214734 4.247459e-01
# 0.285032422 0.756654854   0.376700711 7.063960e-01
# 0.351728440 0.289276108   1.215891772 2.240262e-01
# 0.776562479 0.408677317   1.900184931 5.740886e-02
# 0.372701147 0.302871324   1.230559373 2.184877e-01
# 0.068359052 0.338091872   0.202190758 8.397676e-01
# alpha_h    rho_h
# [1,] 0.003015612 3.598284
# h0        hp     p.bvh
# [1,] -3536.73 -3531.341 -3583.821
# cAIC     rAIC
# [1,] 7151.699 7171.642


p.res10=res10$F.Est[1:(dim(res10$F.Est)[1]/2),]

rownames(p.res10)= c("day_sev_roc_cat1","day_sev_roc_cat2","drips","midaz4", "fent2", "keta4",
                    "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                    "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res10)

#                     beta_h      se_beh     t_value      p_value
# day_sev_roc_cat1    0.76090825 0.244224637   3.1156081 1.835661e-03
# day_sev_roc_cat2   -0.10431621 0.184260162  -0.5661355 5.713017e-01
# drips              -0.49524963 0.062016991  -7.9857088 1.332268e-15
# midaz4             -0.24350898 0.065268281  -3.7308931 1.908021e-04
# fent2              -0.13646792 0.050289569  -2.7136426 6.654792e-03
# keta4               0.04694593 0.036890996   1.2725581 2.031749e-01
# days_paralytics    -0.05293228 0.011178374  -4.7352395 2.187966e-06
# berlin_intubation1  0.71345916 0.253155632   2.8182630 4.828424e-03
# berlin_intubation2  0.63119287 0.234706602   2.6892847 7.160532e-03
# berlin_intubation3  0.73685044 0.238288867   3.0922571 1.986407e-03
# CVVH1              -0.36963809 0.101406176  -3.6451240 2.672629e-04
# age                -0.03520801 0.003358233 -10.4840868 0.000000e+00
# SEX1               -0.23381472 0.093079142  -2.5119991 1.200494e-02
# RACE_ETHNIC1       -0.06639457 0.160270710  -0.4142652 6.786799e-01
# RACE_ETHNIC2        0.63251211 0.308562662   2.0498660 4.037751e-02
# RACE_ETHNIC3        0.20341817 0.128708622   1.5804549 1.140027e-01
# RACE_ETHNIC4       -0.03286004 0.196989266  -0.1668113 8.675185e-01
# RACE_ETHNIC5       -0.28454982 0.169850913  -1.6752917 9.387702e-02
# RACE_ETHNIC6        0.02656733 0.145788386   0.1822322 8.554005e-01

##model 11 and 12
data$CHEMO=as.factor(data$twodays_70)
data_conti=data_surv=data

jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO +(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")

res11<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 39
# beta_h    se_beh   t_value      p_value
# -1.180921 0.1201982 -9.824781 0.000000e+00
# -1.194435 0.2552233 -4.679963 2.869273e-06
# alpha_h    rho_h
# [1,] 0.09633587 2.301741
# h0        hp     p.bvh
# [1,] -3835.651 -3835.871 -3841.914
# cAIC     rAIC
# [1,] 7679.193 7687.827

p.res11=res11$F.Est[1:(dim(res11$F.Est)[1]/2),]
p.res11=matrix(p.res11,nrow = 1)
rownames(p.res11)="twodays70"
colnames(p.res11)=colnames(res11$F.Est)


print(p.res11)

#             beta_h    se_beh   t_value p_value
# twodays70 -1.180921 0.1201982 -9.824781       0


jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
res12<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=200)
# [1] "iterations : "
# [1] 200
# beta_h      se_beh       t_value      p_value
# -0.2614230183 0.136982542  -1.908440414 5.633432e-02
# -0.4923642169 0.061576401  -7.995988805 1.332268e-15
# -0.2189593523 0.065027593  -3.367176039 7.594219e-04
# -0.1525995090 0.050892610  -2.998461072 2.713468e-03
# 0.0581052024 0.036792401   1.579271835 1.142737e-01
# -0.0558493017 0.011252700  -4.963191302 6.934422e-07
# 0.7570123655 0.251219585   3.013349318 2.583813e-03
# 0.6938534495 0.231046454   3.003090661 2.672528e-03
# 0.7904256414 0.235447950   3.357114131 7.876059e-04
# -0.3611196550 0.101764101  -3.548595752 3.872912e-04
# -0.0354666469 0.003350627 -10.585077089 0.000000e+00
# -0.2775252345 0.092651396  -2.995370253 2.741119e-03
# -0.0613764099 0.160199888  -0.383123925 7.016279e-01
# 0.6413347589 0.308603265   2.078185269 3.769230e-02
# 0.2504526923 0.127334660   1.966885475 4.919642e-02
# 0.0006737198 0.195861219   0.003439781 9.972555e-01
# -0.2372444081 0.168562545  -1.407456254 1.592921e-01
# 0.0838193085 0.144291030   0.580904500 5.613048e-01
# -0.6947353581 0.297005194  -2.339135384 1.932843e-02
# -0.6491197563 0.131546829  -4.934514656 8.035031e-07
# -0.6344187776 0.156461430  -4.054793420 5.017863e-05
# 0.0664860393 0.107296844   0.619645806 5.354910e-01
# 0.0982693287 0.088825630   1.106317277 2.685892e-01
# 0.0712734236 0.019165739   3.718793483 2.001766e-04
# 0.6608746198 0.702710061   0.940465573 3.469788e-01
# 0.8790754844 0.648323085   1.355921923 1.751240e-01
# 1.0052887071 0.657926627   1.527964769 1.265213e-01
# 0.1216057493 0.206897201   0.587759278 5.566939e-01
# 0.0223797828 0.008638701   2.590642232 9.579703e-03
# 0.4285872027 0.212726524   2.014733259 4.393261e-02
# 0.2959109670 0.330082990   0.896474450 3.699994e-01
# 0.3134820234 0.756389253   0.414445370 6.785480e-01
# 0.3474383356 0.287752194   1.207422020 2.272697e-01
# 0.7701030006 0.408017768   1.887425157 5.910316e-02
# 0.3684059765 0.302118926   1.219407142 2.226897e-01
# 0.1106324077 0.338233616   0.327088741 7.436008e-01
# alpha_h    rho_h
# [1,] 0.004942117 3.943058
# h0        hp     p.bvh
# [1,] -3542.472 -3537.964 -3590.016
# cAIC     rAIC
# [1,] 7159.743 7184.031

p.res12=res12$F.Est[1:(dim(res12$F.Est)[1]/2),]
rownames(p.res12)= c("twodays70","drips","midaz4", "fent2", "keta4",
                     "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                     "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res12)
#                     beta_h      se_beh       t_value      p_value
# twodays70          -0.2614230183 0.136982542  -1.908440414 5.633432e-02
# drips              -0.4923642169 0.061576401  -7.995988805 1.332268e-15
# midaz4             -0.2189593523 0.065027593  -3.367176039 7.594219e-04
# fent2              -0.1525995090 0.050892610  -2.998461072 2.713468e-03
# keta4               0.0581052024 0.036792401   1.579271835 1.142737e-01
# days_paralytics    -0.0558493017 0.011252700  -4.963191302 6.934422e-07
# berlin_intubation1  0.7570123655 0.251219585   3.013349318 2.583813e-03
# berlin_intubation2  0.6938534495 0.231046454   3.003090661 2.672528e-03
# berlin_intubation3  0.7904256414 0.235447950   3.357114131 7.876059e-04
# CVVH1              -0.3611196550 0.101764101  -3.548595752 3.872912e-04
# age                -0.0354666469 0.003350627 -10.585077089 0.000000e+00
# SEX1               -0.2775252345 0.092651396  -2.995370253 2.741119e-03
# RACE_ETHNIC1       -0.0613764099 0.160199888  -0.383123925 7.016279e-01
# RACE_ETHNIC2        0.6413347589 0.308603265   2.078185269 3.769230e-02
# RACE_ETHNIC3        0.2504526923 0.127334660   1.966885475 4.919642e-02
# RACE_ETHNIC4        0.0006737198 0.195861219   0.003439781 9.972555e-01
# RACE_ETHNIC5       -0.2372444081 0.168562545  -1.407456254 1.592921e-01
# RACE_ETHNIC6        0.0838193085 0.144291030   0.580904500 5.613048e-01

##model 13 and 14

data$CHEMO=data$day_sev_roc
data_conti=data_surv=data
jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO +(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")
res13<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 35
# beta_h      se_beh    t_value     p_value
# -0.08513331 0.005704998 -14.922583 0.00000e+00
# -0.08263303 0.012416868  -6.654901 2.83491e-11
# alpha_h    rho_h
# [1,] 0.5166756 1.369868
# h0        hp     p.bvh
# [1,] -3700.957 -3703.714 -3714.554
# cAIC     rAIC
# [1,] 7409.878 7433.107

p.res13=res13$F.Est[1:(dim(res13$F.Est)[1]/2),]
p.res13=res13$F.Est[1:(dim(res13$F.Est)[1]/2),]
p.res13=matrix(p.res13,nrow = 1)
rownames(p.res13)="day_sev_roc"
colnames(p.res13)=colnames(res11$F.Est)
print(p.res13)

#              beta_h      se_beh   t_value p_value
# day_sev_roc -0.08513331 0.005704998 -14.92258       0


jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
res14<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=200)

# [1] "iterations : "
# [1] 43
# beta_h      se_beh     t_value      p_value
# -0.0639750778 0.006862812 -9.32199196 0.000000e+00
# -0.3957073343 0.063774084 -6.20482976 5.475616e-10
# -0.1853912491 0.066854990 -2.77303533 5.553609e-03
# -0.1077573634 0.055288771 -1.94899183 5.129640e-02
# 0.0631885032 0.041381702  1.52696722 1.267692e-01
# 0.0004612005 0.013257431  0.03478807 9.722487e-01
# 0.8099802591 0.256521050  3.15755864 1.590962e-03
# 0.7449788582 0.234787672  3.17298967 1.508779e-03
# 0.8093356608 0.239233712  3.38303349 7.168989e-04
# -0.4634614288 0.104914314 -4.41752330 9.983831e-06
# -0.0294933639 0.003417551 -8.62997147 0.000000e+00
# -0.3079652028 0.092707276 -3.32190975 8.940361e-04
# 0.0551700679 0.164210969  0.33597066 7.368930e-01
# 0.6431304704 0.311884677  2.06207781 3.920033e-02
# 0.2384376703 0.128788933  1.85138323 6.411444e-02
# -0.0152782747 0.197034841 -0.07754098 9.381932e-01
# -0.3338127424 0.170697652 -1.95557900 5.051477e-02
# 0.1266978588 0.144739331  0.87535197 3.813824e-01
# -0.0700575353 0.015120917 -4.63315373 3.601369e-06
# -0.5998460637 0.137104385 -4.37510487 1.213741e-05
# -0.5768305481 0.156040824 -3.69666434 2.184509e-04
# 0.0744669993 0.112696843  0.66077272 5.087581e-01
# 0.1086193990 0.085124514  1.27600610 2.019534e-01
# 0.1235688556 0.022858012  5.40593189 6.447229e-08
# 0.4826585242 0.697519948  0.69196376 4.889601e-01
# 0.6296025101 0.632420074  0.99554479 3.194714e-01
# 0.7327087834 0.639410553  1.14591287 2.518312e-01
# 0.0187833412 0.211783812  0.08869111 9.293274e-01
# 0.0269905996 0.008712922  3.09776660 1.949849e-03
# 0.3473974028 0.214009833  1.62327777 1.045300e-01
# 0.3885905972 0.326859616  1.18886084 2.344944e-01
# 0.4243847058 0.748752447  0.56678907 5.708575e-01
# 0.3057230500 0.287407632  1.06372627 2.874527e-01
# 0.8114068975 0.407645594  1.99047140 4.653903e-02
# 0.3093383661 0.301428893  1.02623993 3.047785e-01
# 0.1866198530 0.337595495  0.55279130 5.804063e-01
# alpha_h     rho_h
# [1,] 0.175342 0.9442588
# h0        hp     p.bvh
# [1,] -3481.089 -3482.187 -3536.402
# cAIC     rAIC
# [1,] 7037.989 7076.804
p.res14=res14$F.Est[1:(dim(res14$F.Est)[1]/2),]
rownames(p.res14)= c("day_sev_roc","drips","midaz4", "fent2", "keta4",
                     "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                     "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
p.res14=res14$F.Est[1:(dim(res14$F.Est)[1]/2),]
rownames(p.res14)= c("day_sev_roc","drips","midaz4", "fent2", "keta4",
                     "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                     "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res14)

#                      beta_h      se_beh     t_value      p_value
# day_sev_roc        -0.0639750778 0.006862812 -9.32199196 0.000000e+00
# drips              -0.3957073343 0.063774084 -6.20482976 5.475616e-10
# midaz4             -0.1853912491 0.066854990 -2.77303533 5.553609e-03
# fent2              -0.1077573634 0.055288771 -1.94899183 5.129640e-02
# keta4               0.0631885032 0.041381702  1.52696722 1.267692e-01
# days_paralytics     0.0004612005 0.013257431  0.03478807 9.722487e-01
# berlin_intubation1  0.8099802591 0.256521050  3.15755864 1.590962e-03
# berlin_intubation2  0.7449788582 0.234787672  3.17298967 1.508779e-03
# berlin_intubation3  0.8093356608 0.239233712  3.38303349 7.168989e-04
# CVVH1              -0.4634614288 0.104914314 -4.41752330 9.983831e-06
# age                -0.0294933639 0.003417551 -8.62997147 0.000000e+00
# SEX1               -0.3079652028 0.092707276 -3.32190975 8.940361e-04
# RACE_ETHNIC1        0.0551700679 0.164210969  0.33597066 7.368930e-01
# RACE_ETHNIC2        0.6431304704 0.311884677  2.06207781 3.920033e-02
# RACE_ETHNIC3        0.2384376703 0.128788933  1.85138323 6.411444e-02
# RACE_ETHNIC4       -0.0152782747 0.197034841 -0.07754098 9.381932e-01
# RACE_ETHNIC5       -0.3338127424 0.170697652 -1.95557900 5.051477e-02
# RACE_ETHNIC6        0.1266978588 0.144739331  0.87535197 3.813824e-01

##model 15 and 16

data$CHEMO=data$three_cat
data_conti=data_surv=data
jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO +(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+(1|center),RandDist="Gamma")
res15<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=100)

# [1] "iterations : "
# [1] 35
# beta_h     se_beh  t_value      p_value
# 0.7976345 0.09625862 8.286369 2.220446e-16
# 1.5036985 0.16450800 9.140580 0.000000e+00
# 0.2992692 0.21750862 1.375896 1.688538e-01
# 1.4736608 0.34031599 4.330272 1.489254e-05
# alpha_h   rho_h
# [1,] 0.0895743 2.30331
# h0        hp     p.bvh
# [1,] -3823.068 -3823.174 -3830.786
# cAIC     rAIC
# [1,] 7658.008 7665.571
p.res15=res15$F.Est[1:(dim(res15$F.Est)[1]/2),]
p.res15=res15$F.Est[1:(dim(res15$F.Est)[1]/2),]
rownames(p.res15)=c("three_cat1","three_cat2")
print(p.res15)
#            beta_h     se_beh  t_value      p_value
# three_cat1 0.7976345 0.09625862 8.286369 2.220446e-16
# three_cat2 1.5036985 0.16450800 9.140580 0.000000e+00


jm1<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==2)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
jm2<-jointmodeling(Model="mean",RespDist="FM", Link="log",
                   LinPred=Surv(surtime,status==1)~CHEMO+drips +midaz4 +fent2 +keta4 +days_paralytics+berlin_intubation +CVVH +age+ SEX +RACE_ETHNIC+(1|center),RandDist="Gamma")
res16<-jmfit(jm1,jm2,data_conti,data_surv,Maxiter=200)

# [1] "iterations : "
# [1] 62
# beta_h      se_beh     t_value      p_value
# 0.61100298 0.104707040   5.8353572 5.367544e-09
# 0.40162699 0.195037025   2.0592346 3.947177e-02
# -0.52133113 0.062380175  -8.3573207 0.000000e+00
# -0.20756200 0.065637915  -3.1622272 1.565673e-03
# -0.16122702 0.051370954  -3.1384861 1.698230e-03
# 0.04788989 0.039065592   1.2258842 2.202423e-01
# -0.05058122 0.011822881  -4.2782480 1.883701e-05
# 0.72904979 0.253376084   2.8773426 4.010400e-03
# 0.69126651 0.235913672   2.9301672 3.387796e-03
# 0.84498452 0.241123260   3.5043675 4.576932e-04
# -0.32528221 0.101899157  -3.1921972 1.411949e-03
# -0.03453252 0.003395971 -10.1686732 0.000000e+00
# -0.28771091 0.091926248  -3.1298015 1.749245e-03
# 0.02157075 0.162611890   0.1326517 8.944688e-01
# 0.66874280 0.311261309   2.1484932 3.167460e-02
# 0.31427460 0.129538863   2.4261028 1.526195e-02
# 0.12527803 0.199149529   0.6290652 5.293064e-01
# -0.13677705 0.170570481  -0.8018800 4.226224e-01
# 0.15674052 0.146562076   1.0694480 2.848678e-01
# 0.26035354 0.227544372   1.1441880 2.525457e-01
# 0.53839278 0.386468406   1.3931094 1.635867e-01
# -0.65825555 0.131553773  -5.0036995 5.624039e-07
# -0.63315275 0.155501492  -4.0716828 4.667471e-05
# 0.03694806 0.106996389   0.3453207 7.298533e-01
# 0.11208921 0.085564920   1.3099903 1.901991e-01
# 0.07406326 0.018776607   3.9444431 7.998565e-05
# 0.53447392 0.705447243   0.7576384 4.486675e-01
# 0.72143140 0.652177811   1.1061882 2.686451e-01
# 0.84863716 0.660427315   1.2849819 1.987986e-01
# 0.10364752 0.206148184   0.5027816 6.151178e-01
# 0.02135141 0.008611429   2.4794271 1.315936e-02
# 0.38145187 0.213446883   1.7871044 7.392061e-02
# 0.31026725 0.330395945   0.9390771 3.476912e-01
# 0.26170062 0.755042539   0.3466038 7.288890e-01
# 0.39844100 0.287830283   1.3842914 1.662692e-01
# 0.83872452 0.410393654   2.0437073 4.098247e-02
# 0.42848221 0.305304494   1.4034586 1.604801e-01
# 0.09688384 0.339370331   0.2854812 7.752755e-01
# alpha_h    rho_h
# [1,] 0.01937234 1.325019
# h0       hp     p.bvh
# [1,] -3528.975 -3526.55 -3578.091
# cAIC     rAIC
# [1,] 7136.885 7160.181

p.res16=res16$F.Est[1:(dim(res16$F.Est)[1]/2),]
rownames(p.res16)= c("three_cat1","three_cat2","drips","midaz4", "fent2", "keta4",
                     "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                     "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
p.res16=res16$F.Est[1:(dim(res16$F.Est)[1]/2),]
rownames(p.res16)= c("three_cat1","three_cat2","drips","midaz4", "fent2", "keta4",
                     "days_paralytics", "berlin_intubation1", "berlin_intubation2", "berlin_intubation3" ,"CVVH1",
                     "age","SEX1","RACE_ETHNIC1","RACE_ETHNIC2","RACE_ETHNIC3","RACE_ETHNIC4","RACE_ETHNIC5", "RACE_ETHNIC6")
print(p.res16)

#                     beta_h      se_beh     t_value      p_value
# three_cat1          0.61100298 0.104707040   5.8353572 5.367544e-09
# three_cat2          0.40162699 0.195037025   2.0592346 3.947177e-02
# drips              -0.52133113 0.062380175  -8.3573207 0.000000e+00
# midaz4             -0.20756200 0.065637915  -3.1622272 1.565673e-03
# fent2              -0.16122702 0.051370954  -3.1384861 1.698230e-03
# keta4               0.04788989 0.039065592   1.2258842 2.202423e-01
# days_paralytics    -0.05058122 0.011822881  -4.2782480 1.883701e-05
# berlin_intubation1  0.72904979 0.253376084   2.8773426 4.010400e-03
# berlin_intubation2  0.69126651 0.235913672   2.9301672 3.387796e-03
# berlin_intubation3  0.84498452 0.241123260   3.5043675 4.576932e-04
# CVVH1              -0.32528221 0.101899157  -3.1921972 1.411949e-03
# age                -0.03453252 0.003395971 -10.1686732 0.000000e+00
# SEX1               -0.28771091 0.091926248  -3.1298015 1.749245e-03
# RACE_ETHNIC1        0.02157075 0.162611890   0.1326517 8.944688e-01
# RACE_ETHNIC2        0.66874280 0.311261309   2.1484932 3.167460e-02
# RACE_ETHNIC3        0.31427460 0.129538863   2.4261028 1.526195e-02
# RACE_ETHNIC4        0.12527803 0.199149529   0.6290652 5.293064e-01
# RACE_ETHNIC5       -0.13677705 0.170570481  -0.8018800 4.226224e-01
# RACE_ETHNIC6        0.15674052 0.146562076   1.0694480 2.848678e-01

