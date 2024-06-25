
############################ Analysis functions ################################


### get course effects
get_courses_effects<-function(N_courses=8,seed=100000,mu_beta=0.3,sd_beta=0.05){
  if (file.exists("./generated_data/courses_effects.RDS")){
    return(readRDS("./generated_data/courses_effects.RDS"))
  }else{
    library("tmvtnorm")
    set.seed(seed)
    for(x in 1:N_courses){
      beta<-rtmvnorm(n=1,mean = mu_beta,sigma = sd_beta,lower =0,upper = 1)
      if(x==1){#start of the matrix
        courses_effects<-c(beta)
      }else{
        courses_effects<-cbind(courses_effects,c(beta))
      }
    }
    colnames(courses_effects)<-c("c1","c2","c3","c4","c5","c6","c7","c8")
    row.names(courses_effects)<-c("beta")
    saveRDS(courses_effects,"./generated_data/courses_effects.RDS")
    return(courses_effects)
  }
}

### get information criterias
## Model nomenclature
# gm1-> simplest model without interactions
# gm2-> model with interactions
# pm-> extension of the generating models with the propensity scores
get_AIC_BIC<-function(init_seed=1,times,scenario=9){
  library(lme4)
  metrics<-data.frame(rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                      rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                      rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario))
  colnames(metrics)<-c("seed","AIC_gm1","BIC_gm1","AIC_gm2","BIC_gm2","AIC_pm","BIC_pm","scenario","decision")
  for(i in 1:times){
    general_seed=init_seed+(i-1)*10
    for(j in 1:scenario){
      seed=general_seed+j
      dataset<-readRDS(paste0("./generated_data/scenario_",j,"_seed_",seed,".RDS"))
      index=(i-1)*scenario+j
      print(index)
      if(j==1 | j==4 | j==7){#scenarios where the generating model is Model 1
        
        g_model<-lmer(formula=skill ~ c1 + c2 +c3 + c4 + c5 + c6 + c7 + c8 +
                       (1|student),data=dataset)
        p_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                       p1 + p2 + p3 + p4 + p5 + p6 + p7 +p8 + (1|student),data=dataset)
        anova_gm1_pm<-anova(g_model,p_model)
        metrics[index,2]<-as.numeric(anova_gm1_pm$AIC[1])
        metrics[index,3]<-as.numeric(anova_gm1_pm$BIC[1])
        metrics[index,4]<-NA
        metrics[index,5]<-NA
        metrics[index,6]<-as.numeric(anova_gm1_pm$AIC[2])
        metrics[index,7]<-as.numeric(anova_gm1_pm$BIC[2])
        if(metrics$AIC_gm1[index]<=metrics$AIC_pm[index] & metrics$BIC_gm1[index]<=metrics$BIC_pm[index]){
          metrics[index,9]<-"Generating Model"
        }else if(metrics$AIC_gm1[index]>metrics$AIC_pm[index] & metrics$BIC_gm1[index]>metrics$BIC_pm[index]){
          metrics[index,9]<-"Propensity Model"
        }else{
          metrics[index,9]<-"undecided"
        }
      }else{#scenarios where the generating model is Model 2
        g_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                       c1*c2 + c1*c3 + c2*c3 + c1*c2*c3 + c4*c5 + c6*c7 +c6*c8 + c7*c8 + c6*c7*c8 + (1|student),data=dataset)
        p_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        p1 + p2 + p3 + p4 + p5 + p6 + p7 +p8 + c1*c2 + c1*c3 + c2*c3 + c1*c2*c3 + c4*c5 + c6*c7 +c6*c8 + c7*c8 + c6*c7*c8 + (1|student),data=dataset)
        anova_gm2_pm<-anova(g_model,p_model)
        metrics[index,2]<-NA
        metrics[index,3]<-NA
        metrics[index,4]<-as.numeric(anova_gm2_pm$AIC[1])
        metrics[index,5]<-as.numeric(anova_gm2_pm$BIC[1])
        metrics[index,6]<-as.numeric(anova_gm2_pm$AIC[2])
        metrics[index,7]<-as.numeric(anova_gm2_pm$BIC[2])
        if(metrics$AIC_gm2[index]<=metrics$AIC_pm[index] & metrics$BIC_gm2[index]<=metrics$BIC_pm[index]){
          metrics[index,9]<-"Generating Model"
        }else if(metrics$AIC_gm2[index]>metrics$AIC_pm[index] & metrics$BIC_gm2[index]>metrics$BIC_pm[index]){
          metrics[index,9]<-"Propensity Model"
        }else{
          metrics[index,9]<-"undecided"
        }
      }
      metrics[index,1]<-seed
      metrics[index,8]<-j
    }
  }
  return(metrics)
}


### get course effect biases
## Model nomenclature
# gm1-> simplest model without interactions
# gm2-> model with interactions
# pm-> extension of the generating models with the propensity scores
get_biases<-function(init_seed=1,times,scenario=9){
  library(lme4)
  course_effects<-get_courses_effects()
  biases<-data.frame(rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                   rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario))
  colnames(biases)<-c("seed","scenario",
                      "ce1","bg1","bp1",
                      "ce2","bg2","bp2",
                      "ce3","bg3","bp3",
                      "ce4","bg4","bp4",
                      "ce5","bg5","bp5",
                      "ce6","bg6","bp6",
                      "ce7","bg7","bp7",
                      "ce8","bg8","bp8")
  for(i in 1:times){
    general_seed=init_seed+(i-1)*10
    for(j in 1:scenario){
      seed=general_seed+j
      dataset<-readRDS(paste0("./generated_data/scenario_",j,"_seed_",seed,".RDS"))
      index=(i-1)*scenario+j
      print(index)
      biases$seed[index]<-seed
      biases$scenario[index]<-j
      biases$ce1[index]<-course_effects[1]
      biases$ce2[index]<-course_effects[2]
      biases$ce3[index]<-course_effects[3]
      biases$ce4[index]<-course_effects[4]
      biases$ce5[index]<-course_effects[5]
      biases$ce6[index]<-course_effects[6]
      biases$ce7[index]<-course_effects[7]
      biases$ce8[index]<-course_effects[8]
      
      if(j==1 | j==4 | j==7){#scenarios where the generating model is Model 1
        g_model<-lmer(formula=skill ~ c1 + c2 +c3 + c4 + c5 + c6 + c7 + c8 +
                        (1|student),data=dataset)
        p_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        p1 + p2 + p3 + p4 + p5 + p6 + p7 +p8 + (1|student),data=dataset)
      }else{#scenarios where the generating model is Model 2
        g_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        c1*c2 + c1*c3 + c2*c3 + c1*c2*c3 + c4*c5 + c6*c7 +c6*c8 + c7*c8 + c6*c7*c8 + (1|student),data=dataset)
        p_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        p1 + p2 + p3 + p4 + p5 + p6 + p7 +p8 + c1*c2 + c1*c3 + c2*c3 + c1*c2*c3 + c4*c5 + c6*c7 +c6*c8 + c7*c8 + c6*c7*c8 + (1|student),data=dataset)
      }
      sum_g<-summary(g_model)
      sum_p<-summary(p_model)
      
      biases$bg1[index]<- sum_g$coefficients[2,1] - course_effects[1]
      biases$bg2[index]<- sum_g$coefficients[3,1] - course_effects[2]
      biases$bg3[index]<- sum_g$coefficients[4,1] - course_effects[3]
      biases$bg4[index]<- sum_g$coefficients[5,1] - course_effects[4]
      biases$bg5[index]<- sum_g$coefficients[6,1] - course_effects[5]
      biases$bg6[index]<- sum_g$coefficients[7,1] - course_effects[6]
      biases$bg7[index]<- sum_g$coefficients[8,1] - course_effects[7]
      biases$bg8[index]<- sum_g$coefficients[9,1] - course_effects[8]
      
      biases$bp1[index]<- sum_p$coefficients[2,1] + sum_p$coefficients[10,1]*mean(dataset$p1) - course_effects[1]
      biases$bp2[index]<- sum_p$coefficients[3,1] + sum_p$coefficients[11,1]*mean(dataset$p2) - course_effects[2]
      biases$bp3[index]<- sum_p$coefficients[4,1] + sum_p$coefficients[12,1]*mean(dataset$p3) - course_effects[3]
      biases$bp4[index]<- sum_p$coefficients[5,1] + sum_p$coefficients[13,1]*mean(dataset$p4) - course_effects[4]
      biases$bp5[index]<- sum_p$coefficients[6,1] + sum_p$coefficients[14,1]*mean(dataset$p5) - course_effects[5]
      biases$bp6[index]<- sum_p$coefficients[7,1] + sum_p$coefficients[15,1]*mean(dataset$p6) - course_effects[6]
      biases$bp7[index]<- sum_p$coefficients[8,1] + sum_p$coefficients[16,1]*mean(dataset$p7) - course_effects[7]
      biases$bp8[index]<- sum_p$coefficients[9,1] + sum_p$coefficients[17,1]*mean(dataset$p8) - course_effects[8]
    }
  }
  return(biases)
}


get_biases<-function(init_seed=1,times,scenario=9){
  library(lme4)
  course_effects<-get_courses_effects()
  biases<-data.frame(rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario),
                     rep(0,times*scenario),rep(0,times*scenario),rep(0,times*scenario))
  colnames(biases)<-c("seed","scenario",
                      "ce1","bg1","bp1",
                      "ce2","bg2","bp2",
                      "ce3","bg3","bp3",
                      "ce4","bg4","bp4",
                      "ce5","bg5","bp5",
                      "ce6","bg6","bp6",
                      "ce7","bg7","bp7",
                      "ce8","bg8","bp8")
  for(i in 1:times){
    general_seed=init_seed+(i-1)*10
    for(j in 7:9){
      seed=general_seed+j
      dataset<-readRDS(paste0("./generated_data/scenario_",j,"_seed_",seed,".RDS"))
      index=(i-1)*scenario+j
      print(index)
      biases$seed[index]<-seed
      biases$scenario[index]<-j
      biases$ce1[index]<-course_effects[1]
      biases$ce2[index]<-course_effects[2]
      biases$ce3[index]<-course_effects[3]
      biases$ce4[index]<-course_effects[4]
      biases$ce5[index]<-course_effects[5]
      biases$ce6[index]<-course_effects[6]
      biases$ce7[index]<-course_effects[7]
      biases$ce8[index]<-course_effects[8]
      
      if(j==7){#scenarios where the generating model is Model 1
        g_model<-lmer(formula=skill ~ c1 + c2 +c3 + c4 + c5 + c6 + c7 + c8 +
                        (1|student),data=dataset)
        p_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        p1 + p2 + p3 + p4 + p5 + p6 + p7 +p8 + (1|student),data=dataset)
      }else if(j==8 | j==9){#scenarios where the generating model is Model 2
        g_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        c1*c2 + c1*c3 + c2*c3 + c1*c2*c3 + c4*c5 + c6*c7 +c6*c8 + c7*c8 + c6*c7*c8 + (1|student),data=dataset)
        p_model<-lmer(formula=skill ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + 
                        p1 + p2 + p3 + p4 + p5 + p6 + p7 +p8 + c1*c2 + c1*c3 + c2*c3 + c1*c2*c3 + c4*c5 + c6*c7 +c6*c8 + c7*c8 + c6*c7*c8 + (1|student),data=dataset)
      }
      sum_g<-summary(g_model)
      sum_p<-summary(p_model)
      
      biases$bg1[index]<- sum_g$coefficients[2,1] - course_effects[1]
      biases$bg2[index]<- sum_g$coefficients[3,1] - course_effects[2]
      biases$bg3[index]<- sum_g$coefficients[4,1] - course_effects[3]
      biases$bg4[index]<- sum_g$coefficients[5,1] - course_effects[4]
      biases$bg5[index]<- sum_g$coefficients[6,1] - course_effects[5]
      biases$bg6[index]<- sum_g$coefficients[7,1] - course_effects[6]
      biases$bg7[index]<- sum_g$coefficients[8,1] - course_effects[7]
      biases$bg8[index]<- sum_g$coefficients[9,1] - course_effects[8]
      
      biases$bp1[index]<- sum_p$coefficients[2,1] + sum_p$coefficients[10,1]*mean(dataset$p1) - course_effects[1]
      biases$bp2[index]<- sum_p$coefficients[3,1] + sum_p$coefficients[11,1]*mean(dataset$p2) - course_effects[2]
      biases$bp3[index]<- sum_p$coefficients[4,1] + sum_p$coefficients[12,1]*mean(dataset$p3) - course_effects[3]
      biases$bp4[index]<- sum_p$coefficients[5,1] + sum_p$coefficients[13,1]*mean(dataset$p4) - course_effects[4]
      biases$bp5[index]<- sum_p$coefficients[6,1] + sum_p$coefficients[14,1]*mean(dataset$p5) - course_effects[5]
      biases$bp6[index]<- sum_p$coefficients[7,1] + sum_p$coefficients[15,1]*mean(dataset$p6) - course_effects[6]
      biases$bp7[index]<- sum_p$coefficients[8,1] + sum_p$coefficients[16,1]*mean(dataset$p7) - course_effects[7]
      biases$bp8[index]<- sum_p$coefficients[9,1] + sum_p$coefficients[17,1]*mean(dataset$p8) - course_effects[8]
    }
  }
  return(biases)
}


### Bias Analysis ###
library(dplyr)
biases<-get_biases(init_seed=20001,times=2000)
#c1
biases[abs(biases$bg1)<abs(biases$bp1),27]<-"No prop"
biases[abs(biases$bg1)>=abs(biases$bp1),27]<-"prop"
#c2
biases[abs(biases$bg2)<abs(biases$bp2),28]<-"No prop"
biases[abs(biases$bg2)>=abs(biases$bp2),28]<-"prop"
#c3
biases[abs(biases$bg3)<abs(biases$bp3),29]<-"No prop"
biases[abs(biases$bg3)>=abs(biases$bp3),29]<-"prop"
#c4
biases[abs(biases$bg4)<abs(biases$bp4),30]<-"No prop"
biases[abs(biases$bg4)>=abs(biases$bp4),30]<-"prop"
#c5
biases[abs(biases$bg5)<abs(biases$bp5),31]<-"No prop"
biases[abs(biases$bg5)>=abs(biases$bp5),31]<-"prop"
#c6
biases[abs(biases$bg6)<abs(biases$bp6),32]<-"No prop"
biases[abs(biases$bg6)>=abs(biases$bp6),32]<-"prop"
#c7
biases[abs(biases$bg7)<abs(biases$bp7),33]<-"No prop"
biases[abs(biases$bg7)>=abs(biases$bp7),33]<-"prop"
#c8
biases[abs(biases$bg8)<abs(biases$bp8),34]<-"No prop"
biases[abs(biases$bg8)>=abs(biases$bp8),34]<-"prop"

names(biases)[27:34]<-c("dc1","dc2","dc3","dc4","dc5","dc6","dc7","dc8")



c1_freq<-biases %>% group_by(scenario) %>% count(dc1)
c1_freq$n<-c1_freq$n/2000
c2_freq<-biases %>% group_by(scenario) %>% count(dc2)
c2_freq$n<-c2_freq$n/2000
c3_freq<-biases %>% group_by(scenario) %>% count(dc3)
c3_freq$n<-c3_freq$n/2000
c4_freq<-biases %>% group_by(scenario) %>% count(dc4)
c4_freq$n<-c4_freq$n/2000
c5_freq<-biases %>% group_by(scenario) %>% count(dc5)
c5_freq$n<-c5_freq$n/2000
c6_freq<-biases %>% group_by(scenario) %>% count(dc6)
c6_freq$n<-c6_freq$n/2000
c7_freq<-biases %>% group_by(scenario) %>% count(dc7)
c7_freq$n<-c7_freq$n/2000
c8_freq<-biases %>% group_by(scenario) %>% count(dc8)
c8_freq$n<-c8_freq$n/2000

G1_freq<-c1_freq$n+c2_freq$n+c3_freq$n
G1_freq<-G1_freq/3
print(G1_freq)

G2_freq<-c4_freq$n+c5_freq$n
G2_freq<-G2_freq/2
print(G2_freq)

G3_freq<-c6_freq$n+c7_freq$n+c8_freq$n
G3_freq<-G3_freq/3
print(G3_freq)



### Parsimonuosity Analysis ###

library(dplyr)
metrics<-readRDS("metrics.RDS")
metrics[,9]<-rep("undecided",length(metrics))
metrics<-cbind(metrics,rep("undecided",length(metrics)))
colnames(metrics)[9]<-"decision_AIC"
colnames(metrics)[10]<-"decision_BIC"
metrics[!is.na(metrics$AIC_gm1) & metrics$AIC_gm1<=metrics$AIC_pm,9]<-"Generating Model"
metrics[!is.na(metrics$AIC_gm1) & metrics$AIC_gm1>metrics$AIC_pm,9]<-"Propensity Model"
metrics[!is.na(metrics$AIC_gm2) & metrics$AIC_gm2<=metrics$AIC_pm,9]<-"Generating Model"
metrics[!is.na(metrics$AIC_gm2) & metrics$AIC_gm2>metrics$AIC_pm,9]<-"Propensity Model"

metrics[!is.na(metrics$BIC_gm1) & metrics$BIC_gm1<=metrics$BIC_pm,10]<-"Generating Model"
metrics[!is.na(metrics$BIC_gm1) & metrics$BIC_gm1>metrics$BIC_pm,10]<-"Propensity Model"
metrics[!is.na(metrics$BIC_gm2) & metrics$BIC_gm2<=metrics$BIC_pm,10]<-"Generating Model"
metrics[!is.na(metrics$BIC_gm2) & metrics$BIC_gm2>metrics$BIC_pm,10]<-"Propensity Model"


metrics$decision_AIC[metrics$decision_AIC=="Generating Model"]="No prop"
metrics$decision_AIC[metrics$decision_AIC=="Propensity Model"]="Prop"
lattice::histogram(~ as.factor(decision_AIC)| factor(scenario,labels = 
                                                       c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Scenario 5","Scenario 6","Scenario 7","Scenario 8","Scenario 9")),index.cond=list(c(7,8,9,4,5,6,1,2,3)), data = metrics, type="percent",main="Model preference by scenario",
                   xlab="Model Preference",ylab="Percentage")

metrics %>% group_by(scenario,decision_AIC) %>% summarize(frequency = n()/100,.groups = "keep")



### Figure of interest for paper ###


par(mfrow=c(1,3))
N_student=800
val=-0.1
pnorm(val,mean = 0,sd=0.2,lower.tail = TRUE)
set.seed(1000)
den<-density(x = rnorm(n=N_student,mean=0,sd=0.2))
plot(den,xlab="interest",main="")
polygon(c(den$x[den$x <= val], val),
        c(den$y[den$x <= val], 0),
        col = "lightblue",
        border = 1)
curve(den)
polygon(c(den$x[den$x >= val], val),
        c(den$y[den$x >= val], 0),
        col = "lightgreen",
        border = 1)
legend(x=-0.5,y=1.9,legend=c("probG1=30.8%","probG3=69.2%"),col=c("lightblue","lightgreen","white"),lty=1:1,cex=1)

val=0
plot(den,xlab="interest",main="")
polygon(c(den$x[den$x <= val], val),
        c(den$y[den$x <= val], 0),
        col = "lightblue",
        border = 1)
curve(den)
polygon(c(den$x[den$x >= val], val),
        c(den$y[den$x >= val], 0),
        col = "lightgreen",
        border = 1)
legend(x=-0.5,y=1.9,legend=c("probG1=50%","probG3=50%"),col=c("lightblue","lightgreen","white"),lty=1:1,cex=1)

val=0.2
pnorm(val,mean = 0,sd=0.2,lower.tail = TRUE)
plot(den,xlab="interest",main="")
polygon(c(den$x[den$x <= val], val),
        c(den$y[den$x <= val], 0),
        col = "lightblue",
        border = 1)
curve(den)
polygon(c(den$x[den$x >= val], val),
        c(den$y[den$x >= val], 0),
        col = "lightgreen",
        border = 1)
legend(x=-0.5,y=1.9,legend=c("probG1=84.1%","probG3=15.9%"),col=c("lightblue","lightgreen","white"),lty=1:1,cex=1)

##### End of figure of interest ######

