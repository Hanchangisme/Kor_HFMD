### set directory

files<-list.files('DATASET/')
data_directory=paste0(getwd(),'/DATASET/')
result_folder='result'
if (file.exists(result_folder)) {
  print("The folder already exists")
  print("This folder includes annual counterfactual effects")
} else {
  dir.create(result_folder)
  print('The folder created')
  print("This folder includes annual counterfactual effects")
}

###  regression for years
### put years you want 
years<-list(2014,2015,2016,2017,2018,2019)

for (year in years){
  x<-sprintf('dataset_%s_epidemic_period.csv',year)
  x<-paste0(data_directory,x)
  df<-read.csv(x)
  
  y<-sprintf('dataset_%s_epidemic_weekly.csv',year)
  y<-paste0(data_directory,y)
  df2<-read.csv(y)
  
  z<-sprintf('dataset_%s_weekly.csv',year)
  z<-paste0(data_directory,z)
  df3<-read.csv(z)
  
  ### regression model
  if (year!=2015){
    model1 <- lm(log(Rt) ~ Cu_hfmd+schoolvac, data =df)
  } else {
    model1 <- lm(log(Rt) ~ Cu_hfmd+schoolvac+intervention, data =df)
  }
  a<-summary(model1)
  
  if (year!=2015){
    model2 <- lm(log(Rt) ~ Cu_hfmd+Rainyseason, data =df)
  } else {
    model2 <- lm(log(Rt) ~ Cu_hfmd+Rainyseason+intervention, data =df)
  }
  b<-summary(model2)
  
  if (year!=2015){
    model3 <- lm(log(Rt) ~ Cu_hfmd+schoolvac+Rainyseason,data =df)
    cat("currnet year: ",year,'\n')
  } else {
    model3 <-lm(log(Rt) ~ Cu_hfmd+schoolvac+Rainyseason+intervention,data=df)
    cat("currnet year: ",year,'\n')
  }
  c<-summary(model3)
  
  models<-list(a,b,c)
  
  file<-sprintf('reg_result_%s.txt',year)
  for (model in models){
    cat(year,"\n",toString(model[1]),"\n",file=file,append=TRUE)
    for (no in c(4,8,9,10,11)){
      write.table(model[no],file=file,append=TRUE)
      cat("\n",file=file,append=TRUE)
    }
  }
  
  est<-as.data.frame(c[4])[3,1]
  est_ci<-as.data.frame(confint(model3))[3,1]
  est_CI<-as.data.frame(confint(model3))[3,2]
  cat(
    "ESTIMATE school vac coeff\n",
    est,
    "\n 95CI under \n",
    est_ci,
    "\n 95CI up \n",
    est_CI,
    "\n",
    file=file,append=TRUE
  )
  reg_coef_mat_phsm=c(0.0,est,est_ci,est_CI)
  phsm=df2$schoolvac
  
  ### rain coefficient
  rain_est<-as.data.frame(c[4])[4,1]
  rain_est_ci<-as.data.frame(confint(model3))[4,1]
  rain_est_CI<-as.data.frame(confint(model3))[4,2]
  cat(
    "ESTIMATE rain \n",
    rain_est,
    "\n 95CI under \n",
    rain_est_ci,
    "\n 95CI up \n",
    rain_est_CI,
    "\n",
    file=file,append=TRUE
  )
  reg_coef_mat_rain=c(0.0, rain_est, rain_est_ci, rain_est_CI)
  rain=df2$Rainyseason
  
  if (year==2015){
    est2<-as.data.frame(c[4])[5,1]
    est2_ci<-as.data.frame(confint(model3))[5,1]
    est2_CI<-as.data.frame(confint(model3))[5,2]
    
    cat(
      "ESTIMATE intervention \n",
      est2,
      "\n 95CI under \n",
      est2_ci,
      "\n 95CI up \n",
      est2_CI,
      "\n",
      file=file,append=TRUE
    )
    
    reg_coef_mat_int=c(0.0,est2,est2_ci,est2_CI)
    intervention=df2$intervention
  } else{
  
  }
  
  "S0"<-1-(df2$cases_total[1]/df3$Cu_hfmd[length(df3$Cu_hfmd)])
  "E0"<-0.00015
  "I0"<-df2$cases_total[1]/df3$Cu_hfmd[length(df3$Cu_hfmd)]
  "R0"<-1-(S0+E0+I0)
  #sigma: Rate of exposed individuals become infected
  "sigma"=0.57
  #gamma: recovery rate
  "gamma"=0.50
  
  "R_0"=2.55
  
  beta0=(S0*R_0/S0)*gamma
  beta=beta0
  
  np=length(reg_coef_mat_phsm)
  n=length(phsm)
  #Setting for model output 
  S = c()
  E = c()
  I = c()
  R = c()
  t = c()
  N = c()
  Inc  = c()
  Rt = c()
  betaT=c()
  #SEIR model and save output
  for(j in 1:np){
    dt=1
    S[1]<-S0
    E[1]=E0
    I[1]=I0
    R[1]=R0
    t[1]=0
    Inc[1]=0
    N[1]=S[1]+E[1]+I[1]+R[1] 
    betaT[1]=beta
    Rt[1]=beta/gamma
    for(i in 2:n){
      #betaT[i]=beta*(exp(reg_coef_mat_hdsc[j] *hday[i]))
      
      #betaT[i]= beta + reg_coef_mat_Hum[j]*Hum[i] - reg_coef_mat_hdsc[j]*hday[i]
      
      counter_rain=(exp(reg_coef_mat_rain[2]*rain[i]))
      counter_phsm=(exp(reg_coef_mat_phsm[j]*phsm[i]))
      if(year!=2015){
        betaT[i]= beta*counter_rain*counter_phsm
      } else{
        counter_inter=(exp(reg_coef_mat_int[2]*intervention[i]))
        betaT[i]= beta*counter_rain*counter_phsm*counter_inter
      }
    }
    for(i in 2:n-1){
      dS=-betaT[i]*I[i]*S[i]
      dE=betaT[i]*I[i]*S[i]-sigma*E[i]
      dI=sigma*E[i]-gamma*I[i]
      dR=gamma*I[i]
      S[i+1]=S[i]+dt*dS
      E[i+1]=E[i]+dt*dE
      I[i+1]=I[i]+dt*dI
      R[i+1]=R[i]+dt*dR
      t[i+1]=t[i]+dt
    }
    for(i in 2:n){
      Inc[i]=betaT[i]*I[i]*S[i]
      Rt[i]=betaT[i]*S[i]/gamma
    }
    
    matrix=data.frame(time=t,s=S,e=E,i=I,r=R,inc=Inc,rt=Rt,real_inc=Inc*df3$Cu_hfmd[length(df3$Cu_hfmd)])
    form=sprintf('/result/SEIR_Rt_%s_%s.csv',year,j)
    form=paste0(getwd(),form)
    write.csv(matrix[2:n,],file=form)
  }
  file1=sprintf("result/SEIR_Rt_%s_1.csv",year)
  file2=sprintf("result/SEIR_Rt_%s_2.csv",year)
  file3=sprintf("result/SEIR_Rt_%s_3.csv",year)
  file4=sprintf("result/SEIR_Rt_%s_4.csv",year)
  
  SEIR_without_phsm_rain <-read.csv(file1)
  SEIR_with_phsm_rain <-read.csv(file2)
  SEIR_with_phsm_rainLb <-read.csv(file3)
  SEIR_with_phsm_rainUb <-read.csv(file4)
  
  without_phsm<-SEIR_without_phsm_rain$real_inc
  with_phsm_est<-SEIR_with_phsm_rain$real_inc
  with_phsm_est_ci<-SEIR_with_phsm_rainLb$real_inc
  with_phsm_est_CI<-SEIR_with_phsm_rainUb$real_inc
  
  
  final=sprintf("Counterfactual_Result_%s.csv",year)
  if(year!=2015){
    result<-data.frame(df2$date[2:n],
                       df2$schoolvac[2:n],
                       df2$Rainyseason[2:n],
                       without_phsm,
                       with_phsm_est,
                       with_phsm_est_ci,
                       with_phsm_est_CI
    )
  } else{
    result<-data.frame(df2$date[2:n],
                       df2$schoolvac[2:n],
                       df2$Rainyseason[2:n],
                       df2$intervention[2:n],
                       without_phsm,
                       with_phsm_est,
                       with_phsm_est_ci,
                       with_phsm_est_CI
    )
  }
  
  write.csv(result,file=final)
}