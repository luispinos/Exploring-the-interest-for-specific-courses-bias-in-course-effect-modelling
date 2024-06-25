
######################## Data Generation functions #############################

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

### get interaction effects
get_interactions<-function(){
  if (file.exists(paste0("./generated_data/interaction_effects.RDS"))){
    return(readRDS(paste0("./generated_data/interaction_effects.RDS")))
  }else{
    interactions<-data.frame(rbind(c(0,0.025,0.05),c(0,0.05,0.1)))
    colnames(interactions)<-c("no_effect","small_effect","large_effect")
    rownames(interactions)<-c("double_int","triple_int")
    saveRDS(interactions,paste0("./generated_data/interaction_effects.RDS"))
    return(interactions)
  }
}

### get confounders
get_confounders<-function(N_student=800,seed,mu_interest=0,sd_interest=0.2){
  if (file.exists(paste0("./generated_data/student_confounders_seed_",seed,".RDS"))){
    return(readRDS(paste0("./generated_data/student_confounders_seed_",seed,".RDS")))
  }else{
    library("MASS")
    set.seed(seed)
    student_confounders<-data.frame(student=seq(1,N_student),stage=rep(0,times=N_student),interest=rep(0,times=N_student))
    student_confounders<-rbind(student_confounders,data.frame(student=seq(1,N_student),stage=rep(1,times=N_student),interest=rnorm(n = N_student, mean=mu_interest,sd=sd_interest)))
    student_confounders<-rbind(student_confounders,data.frame(student=seq(1,N_student),stage=rep(2,times=N_student),interest=rnorm(n = N_student, mean=mu_interest,sd=sd_interest)))
    saveRDS(student_confounders,paste0("./generated_data/student_confounders_seed_",seed,".RDS"))
    return(student_confounders)
  }
}

### fixed effect of interest over skill
get_interest_effect<-function(confounder_value,scenario){
  if(scenario==4 | scenario==5 | scenario==6){#scenarios with small effect of interest
    return(1*confounder_value)
  }else{#scenarios 7,8 and 9 (large effect of interest)
    return(3*confounder_value)
  }
}

### get all combinations for stage 1
get_combinations_stage1<-function(N_courses=8,N_courses_stage1=3){
  if (file.exists(paste0("./generated_data/combinations_stage1.RDS"))){
    return(readRDS(paste0("./generated_data/combinations_stage1.RDS")))
  }else{
    all_combinations<-expand.grid(rep(list(0:1),N_courses))
    all_combinations<-cbind(all_combinations,rowSums(all_combinations),rowSums(all_combinations[,1:3]),rowSums(all_combinations[,6:8]))
    combinations_stage1<-all_combinations[all_combinations$`rowSums(all_combinations)`==N_courses_stage1,]
    colnames(combinations_stage1)<-c("c1","c2","c3","c4","c5","c6","c7","c8","sumT","sumG1","sumG3")
    saveRDS(combinations_stage1,paste0("./generated_data/combinations_stage1.RDS"))
    return(combinations_stage1)
  }
}

### get all combinations for stage 2
get_combinations_stage2<-function(N_courses=8,N_courses_stage2=5){
  if (file.exists(paste0("./generated_data/combinations_stage2.RDS"))){
    return(readRDS(paste0("./generated_data/combinations_stage2.RDS")))
  }else{
    all_combinations<-expand.grid(rep(list(0:1),N_courses))
    all_combinations<-cbind(all_combinations,rowSums(all_combinations),rowSums(all_combinations[,1:3]),rowSums(all_combinations[,6:8]))
    combinations_stage2<-all_combinations[all_combinations$`rowSums(all_combinations)`==N_courses_stage2,]
    colnames(combinations_stage2)<-c("c1","c2","c3","c4","c5","c6","c7","c8","sumT","sumG1","sumG3")
    saveRDS(combinations_stage2,paste0("./generated_data/combinations_stage2.RDS"))
    return(combinations_stage2)
  }
}

### get empty student skills data frame
get_empty_student_skills<-function(N_student=800){
  if (file.exists(paste0("./generated_data/empty_student_skills.RDS"))){
    return(readRDS(paste0("./generated_data/empty_student_skills.RDS")))
  }else{
    student_skills<-data.frame(student=seq(1,N_student),stage=rep(0,times=N_student),skill=rep(0,times=N_student))
    student_skills<-rbind(student_skills,data.frame(student=seq(1,N_student),stage=rep(1,times=N_student),skill=rep(0,times=N_student)))
    student_skills<-rbind(student_skills,data.frame(student=seq(1,N_student),stage=rep(2,times=N_student),skill=rep(0,times=N_student)))
    saveRDS(student_skills,paste0("./generated_data/empty_student_skills.RDS"))
    return(student_skills)
  }
}

### build datasets
build_datasets<-function(N_student=800,seed,N_courses=8,N_courses_stage1=3,N_courses_stage2=5,
                         sd_re_student=0.25,sd_re_course=0.03,sd_residual=0.2,times,scenario,beta_0=-1.5){
  courses_effects<-get_courses_effects()
  interactions<-get_interactions()
  stage<-c(rep(0,N_student),rep(1,N_student),rep(2,N_student))
  combinations_stage1<-get_combinations_stage1()
  combinations_stage2<-get_combinations_stage2()
  student_skills<-get_empty_student_skills()
  for(i in 1:times){
    building_seed=seed+(i-1)*10
    student_confounders<-get_confounders(seed=building_seed)
    if(scenario==1 | scenario==2 | scenario==3){
      build_datasets_no_interest(seed=building_seed+scenario,student_skills = student_skills,student_confounders = student_confounders,
                       combinations_stage1 = combinations_stage1,combinations_stage2 = combinations_stage2,beta_0 = beta_0,
                       courses_effects = courses_effects,stage=stage,scenario=scenario)
    }else if(scenario==4 | scenario==5 | scenario==6 | scenario==7 | scenario==8 | scenario==9){
      build_datasets_interest(seed=building_seed+scenario,student_skills = student_skills,student_confounders = student_confounders,
                       combinations_stage1 = combinations_stage1,combinations_stage2 = combinations_stage2,beta_0 = beta_0,
                       courses_effects = courses_effects,stage=stage,scenario=scenario)
    }else{
      print("Not valid scenario")
    }
  }
}

### build datasets with no interest
build_datasets_no_interest<-function(N_student=800,seed,N_courses=8,sd_re_student=0.25,sd_re_course=0.03,sd_residual=0.2,re_course=0,
                       student_skills,student_confounders,combinations_stage1,combinations_stage2,beta_0=-1.5,
                       courses_effects,stage,scenario){
  set.seed(seed)
  courses_stage0<-rep(0,N_courses)
  courses_stage1<-NULL
  courses_stage2<-NULL
  for(i in 1:(N_student-1)){
    courses_stage0<-rbind(courses_stage0,rep(0,N_courses))  
  }
  courses_stage0<-as.data.frame(courses_stage0)
  
  r_effect_st<-rnorm(n = N_student,mean=0,sd=sd_re_student)
  residual<-rnorm(n=nrow(student_skills),mean=0,sd=sd_residual)
  if(re_course==1){
    r_effect_c<-rnorm(n = N_courses,mean=0,sd=0.1)
  }else{
    r_effect_c<-rep(0,N_courses)
  }
  for(i in 1:nrow(student_skills)){
    if(i==N_student+1){#first student stage 1
      chosen_sequences=matrix(rep(0,N_student*2),ncol=2)
    }
    if(i>=1 & i<=N_student){#stage 0
      index=i
      combinations=matrix(rep(0,N_courses),ncol=N_courses)
    }
    if(i>=(N_student+1) & i<=(N_student*2)){#stage 1
      index=i-N_student
      combinations=combinations_stage1
    }
    if(i>=(1+N_student*2) & i<=nrow(student_skills)){#stage 2
      index=i-N_student*2
      combinations<-NULL
      current_student=i-N_student*2
      courses_stage1_current_student=combinations_stage1[chosen_sequences[current_student,1],]
      for(j in 1:nrow(combinations_stage2)){
        if(sum(courses_stage1_current_student[1:N_courses]==combinations_stage2[j,1:N_courses])==6){
          if(is.null(combinations)){
            combinations<-combinations_stage2[j,]
          }else{
            combinations<-rbind(combinations,combinations_stage2[j,])
          }
        }
      }
    }
    course_sequence=combinations[sample(nrow(combinations),size=1),]
    student_skills$skill[i]=as.numeric(beta_0 + 
                                         course_sequence[1]*(courses_effects[1,1] + r_effect_c[1]) +
                                         course_sequence[2]*(courses_effects[1,2] + r_effect_c[2]) +
                                         course_sequence[3]*(courses_effects[1,3] + r_effect_c[3]) +
                                         course_sequence[4]*(courses_effects[1,4] + r_effect_c[4]) +
                                         course_sequence[5]*(courses_effects[1,5] + r_effect_c[5]) +
                                         course_sequence[6]*(courses_effects[1,6] + r_effect_c[6]) +
                                         course_sequence[7]*(courses_effects[1,7] + r_effect_c[7]) +
                                         course_sequence[8]*(courses_effects[1,8] + r_effect_c[8]) +
                                         r_effect_st[index] + residual[i])
    if(scenario==2 | scenario==3){#interactions
      interactions<-get_interactions()
      student_skills$skill[i]=as.numeric(student_skills$skill[i]+
        course_sequence[1]*course_sequence[2]*(interactions[1,scenario]) +
        course_sequence[1]*course_sequence[3]*(interactions[1,scenario]) +
        course_sequence[2]*course_sequence[3]*(interactions[1,scenario]) +
        course_sequence[4]*course_sequence[5]*(interactions[1,scenario]) +
        course_sequence[6]*course_sequence[7]*(interactions[1,scenario]) +
        course_sequence[6]*course_sequence[8]*(interactions[1,scenario]) +
        course_sequence[7]*course_sequence[8]*(interactions[1,scenario]) +
        course_sequence[1]*course_sequence[2]*course_sequence[3]*(interactions[2,scenario]) +
        course_sequence[6]*course_sequence[7]*course_sequence[8]*(interactions[2,scenario]))
    }
    if(i>=(N_student+1) & i<=(N_student*2)){#stage 1
      chosen_sequences[index,1]=which(rownames(course_sequence)==rownames(combinations_stage1))
    }
    if(i>=(1+N_student*2) & i<=nrow(student_skills)){#stage 2
      chosen_sequences[index,2]=which(rownames(course_sequence)==rownames(combinations_stage2))
    }
  }
  courses_stage1<-combinations_stage1[chosen_sequences[,1],1:N_courses]
  courses_stage2<-combinations_stage2[chosen_sequences[,2],1:N_courses]
  
  colnames(courses_stage0)<-c("c1","c2","c3","c4","c5","c6","c7","c8")
  colnames(courses_stage1)<-names(courses_stage0)
  colnames(courses_stage2)<-names(courses_stage0)
  courses_stages<-data.frame(rbind(courses_stage0,courses_stage1,courses_stage2),stage)
  row.names(courses_stages)<-seq(1,nrow(courses_stages))
  
  dataset<-data.frame(student_skills$student,student_skills$skill,courses_stages)
  colnames(dataset)[1:2]<-c("student","skill")
  dataset<-add_propensity_scores(dataset=dataset)
  saveRDS(dataset,paste0("./generated_data/scenario_",scenario,"_seed_",seed,".RDS"))
}

### build datasets with interest
build_datasets_interest<-function(N_student=800,seed,N_courses=8,sd_re_student=0.25,sd_re_course=0.03,sd_residual=0.2,re_course=0,
                         student_skills,student_confounders,combinations_stage1,combinations_stage2,beta_0=-1.5,
                         courses_effects,stage,scenario){
  set.seed(seed)
  courses_stage0<-rep(0,N_courses)
  courses_stage1<-NULL
  courses_stage2<-NULL
  for(i in 1:(N_student-1)){
    courses_stage0<-rbind(courses_stage0,rep(0,N_courses))  
  }
  courses_stage0<-as.data.frame(courses_stage0)
  
  r_effect_st<-rnorm(n = N_student,mean=0,sd=sd_re_student)
  residual<-rnorm(n=nrow(student_skills),mean=0,sd=sd_residual)
  if(re_course==1){
    r_effect_c<-rnorm(n = N_courses,mean=0,sd=0.1)
  }else{
    r_effect_c<-rep(0,N_courses)
  }
  for(i in 1:nrow(student_skills)){
    if(i==N_student+1){#first student stage 1
      chosen_sequences=matrix(rep(0,N_student*2),ncol=2)
    }
    if(i>=1 & i<=N_student){#stage 0
      index=i
      combinations=matrix(rep(0,N_courses),ncol=N_courses)
      course_sequence=combinations[sample(nrow(combinations),size=1),]
    }
    if(i>=(N_student+1) & i<=(N_student*2)){#stage 1
      index=i-N_student
      combinations=combinations_stage1
      combinations$rowNumber<-seq.int(nrow(combinations))#rowNumber on combinations
      index_G1<-combinations$rowNumber[combinations$sumG1>=combinations$sumG3]#rows which G1 is chosen more or equal than G3
      index_G3<-combinations$rowNumber[combinations$sumG1<combinations$sumG3]#rows which G1 is chosen less than G3
      pG1<-pnorm(student_confounders$interest[i],mean = mean(student_confounders$interest),sd=sd(student_confounders$interest),lower.tail = TRUE)
      if(length(index_G1)==0 | length(index_G3==0)){
        course_sequence=combinations[sample(nrow(combinations),size=1),]
      }else{
        prob_combinations<-rep(0,nrow(combinations))
        prob_combinations[index_G1]=pG1/length(index_G1)
        prob_combinations[index_G3]=(1-pG1)/length(index_G3)
        course_sequence=combinations[sample(nrow(combinations), size=1,prob=prob_combinations),]
      }
    }
    if(i>=(1+N_student*2) & i<=nrow(student_confounders)){#stage 2
      index=i-N_student*2
      combinations<-NULL
      current_student=i-N_student*2
      courses_stage1_current_student=combinations_stage1[chosen_sequences[current_student,1],]
      for(j in 1:nrow(combinations_stage2)){
        if(sum(courses_stage1_current_student[1:N_courses]==combinations_stage2[j,1:N_courses])==6){
          if(is.null(combinations)){
            combinations<-combinations_stage2[j,]
          }else{
            combinations<-rbind(combinations,combinations_stage2[j,])
          }
        }
      }
      combinations$rowNumber<-seq.int(nrow(combinations))#rowNumber on combinations
      index_G1<-combinations$rowNumber[combinations$sumG1>=combinations$sumG3]#rows which G1 is chosen more or equal than G3
      index_G3<-combinations$rowNumber[combinations$sumG1<combinations$sumG3]#rows which G1 is chosen less than G3
      pG1<-pnorm(student_confounders$interest[i],mean = mean(student_confounders$interest),sd=sd(student_confounders$interest),lower.tail = TRUE)
      if(length(index_G1)==0 | length(index_G3==0)){
        course_sequence=combinations[sample(nrow(combinations),size=1),]
      }else{
        prob_combinations<-rep(0,nrow(combinations))
        prob_combinations[index_G1]=pG1/length(index_G1)
        prob_combinations[index_G3]=(1-pG1)/length(index_G3)
        course_sequence=combinations[sample(nrow(combinations), size=1,prob=prob_combinations),]
      }
    }
    student_skills$skill[i]=as.numeric(beta_0 + 
                                         course_sequence[1]*(courses_effects[1,1] + r_effect_c[1]) +
                                         course_sequence[2]*(courses_effects[1,2] + r_effect_c[2]) +
                                         course_sequence[3]*(courses_effects[1,3] + r_effect_c[3]) +
                                         course_sequence[4]*(courses_effects[1,4] + r_effect_c[4]) +
                                         course_sequence[5]*(courses_effects[1,5] + r_effect_c[5]) +
                                         course_sequence[6]*(courses_effects[1,6] + r_effect_c[6]) +
                                         course_sequence[7]*(courses_effects[1,7] + r_effect_c[7]) +
                                         course_sequence[8]*(courses_effects[1,8] + r_effect_c[8]) +
                                         get_interest_effect(confounder_value = as.numeric(student_confounders$interest[i]),
                                         scenario=scenario)
                                         + r_effect_st[index] + residual[i])
    if(scenario==5 | scenario==6 | scenario==8 | scenario==9){#interactions
      interactions<-get_interactions()
      if(scenario==5 | scenario==8){
        small_large=0
      }else{
        small_large=1
      }
      student_skills$skill[i]=as.numeric(student_skills$skill[i]+
                                           course_sequence[1]*course_sequence[2]*(interactions[1,2+small_large]) +
                                           course_sequence[1]*course_sequence[3]*(interactions[1,2+small_large]) +
                                           course_sequence[2]*course_sequence[3]*(interactions[1,2+small_large]) +
                                           course_sequence[4]*course_sequence[5]*(interactions[1,2+small_large]) +
                                           course_sequence[6]*course_sequence[7]*(interactions[1,2+small_large]) +
                                           course_sequence[6]*course_sequence[8]*(interactions[1,2+small_large]) +
                                           course_sequence[7]*course_sequence[8]*(interactions[1,2+small_large]) +
                                           course_sequence[1]*course_sequence[2]*course_sequence[3]*(interactions[2,2+small_large]) +
                                           course_sequence[6]*course_sequence[7]*course_sequence[8]*(interactions[2,2+small_large]))
    }
    if(i>=(N_student+1) & i<=(N_student*2)){#stage 1
      chosen_sequences[index,1]=which(rownames(course_sequence)==rownames(combinations_stage1))
    }
    if(i>=(1+N_student*2) & i<=nrow(student_skills)){#stage 2
      chosen_sequences[index,2]=which(rownames(course_sequence)==rownames(combinations_stage2))
    }
  }
  
  courses_stage1<-combinations_stage1[chosen_sequences[,1],1:N_courses]
  courses_stage2<-combinations_stage2[chosen_sequences[,2],1:N_courses]
  
  colnames(courses_stage0)<-c("c1","c2","c3","c4","c5","c6","c7","c8")
  colnames(courses_stage1)<-names(courses_stage0)
  colnames(courses_stage2)<-names(courses_stage0)
  courses_stages<-data.frame(rbind(courses_stage0,courses_stage1,courses_stage2),stage)
  row.names(courses_stages)<-seq(1,nrow(courses_stages))
  
  dataset<-data.frame(student_skills$student,student_skills$skill,courses_stages)
  colnames(dataset)[1:2]<-c("student","skill")
  dataset<-add_propensity_scores(dataset=dataset)
  saveRDS(dataset,paste0("./generated_data/scenario_",scenario,"_seed_",seed,".RDS"))
}

### add propensity scores
add_propensity_scores<-function(N_student=800,dataset){
  total<-nrow(dataset)
  library(MatchIt)
  temp_dataset<-dataset[N_student+1:total,]
  #propensity score course 1
  m.c1 = matchit(c1 ~ c2 + c3 + c4 + c5 + c6 + c7 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")
  p1<-rbind(rep(0,times=N_student),m.c1$distance)
  dataset$p1<-p1
  
  #propensity score course 2
  m.c2 = matchit(c2 ~ c1 + c3 + c4 + c5 + c6 + c7 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")
  
  p2<-rbind(rep(0,times=N_student),m.c2$distance)
  dataset$p2<-p2
  
  #propensity score course 3
  m.c3 = matchit(c3 ~ c1 + c2 + c4 + c5 + c6 + c7 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")

  p3<-rbind(rep(0,times=N_student),m.c3$distance)
  dataset$p3<-p3
  
  #propensity score course 4
  m.c4 = matchit(c4 ~ c1 + c2 + c3 + c5 + c6 + c7 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")
  
  p4<-rbind(rep(0,times=N_student),m.c4$distance)
  dataset$p4<-p4
  
  #propensity score course 5
  m.c5 = matchit(c5 ~ c1 + c2 + c3 + c4 + c6 + c7 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")
  
  p5<-rbind(rep(0,times=N_student),m.c5$distance)
  dataset$p5<-p5
  
  #propensity score course 6
  m.c6 = matchit(c6 ~ c1 + c2 + c3 + c4 + c5 + c7 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")

  p6<-rbind(rep(0,times=N_student),m.c6$distance)
  dataset$p6<-p6
  
  #propensity score course 7
  m.c7 = matchit(c7 ~ c1 + c2 + c3 + c4 + c5 + c6 + c8,
                 data = temp_dataset, method = "nearest",distance="glm")

  p7<-rbind(rep(0,times=N_student),m.c7$distance)
  dataset$p7<-p7
  
  #propensity score course 8
  m.c8 = matchit(c8 ~ c1 + c2 + c3 + c4 + c5 + c6 + c7,
                 data = temp_dataset, method = "nearest",distance="glm")
  
  p8<-rbind(rep(0,times=N_student),m.c8$distance)
  dataset$p8<-p8
  
  #dataset$p1[1:N_student]=0
  dataset$p1[dataset$c1==0]=0
  #dataset$p2[1:N_student]=0
  dataset$p2[dataset$c2==0]=0
  #dataset$p3[1:N_student]=0
  dataset$p3[dataset$c3==0]=0
  #dataset$p4[1:N_student]=0
  dataset$p4[dataset$c4==0]=0
  #dataset$p5[1:N_student]=0
  dataset$p5[dataset$c5==0]=0
  #dataset$p6[1:N_student]=0
  dataset$p6[dataset$c6==0]=0
  #dataset$p7[1:N_student]=0
  dataset$p7[dataset$c7==0]=0
  #dataset$p8[1:N_student]=0
  dataset$p8[dataset$c8==0]=0
  
  return(dataset)
  
}

