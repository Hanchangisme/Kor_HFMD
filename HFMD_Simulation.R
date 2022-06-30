### set directory

files<-list.files('DATASET/')
data_directory=paste0(getwd(),'/DATASET/')
result_folder='simulation'
if (file.exists(result_folder)) {
  print("The folder already exists")
} else {
  dir.create(result_folder)
  print('The folder created')
}

### put years you want 
years<-list(2022)

for (year in years){

  y<-sprintf('dataset_%s_epidemic_weekly.csv',year)
  y<-paste0(data_directory,y)
  df2<-read.csv(y)
  
  z<-sprintf('dataset_%s_weekly.csv',year)
  z<-paste0(data_directory,z)
  df3<-read.csv(z)
  
  "S0"<-0.99999999
  "E0"<-0.001
  "I0"<-3000/700000
  "R0"<-1-(S0+E0+I0)
  #sigma: Rate of exposed individuals become infected
  "sigma"=0.57
  #gamma: recovery rate
  "gamma"=0.15
  
  "R_0"=2.5
  
  beta0=(S0*R_0/S0)*gamma
  beta=beta0
  
  reg_coef_mat_phsm=c(0.0,-0.0457,-0.0299,-0.0634)
  np=length(reg_coef_mat_phsm)
  
  phsm=df2$schoolvac
  n=length(phsm)
  
  ### rain coefficient
  
  reg_coef_mat_rain=c(0.0, 0.0377)
  rain=df2$Rainyseason
  
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
      betaT[i]= beta*(exp(reg_coef_mat_rain[2]*(rain[i])))*(exp(reg_coef_mat_phsm[j]*phsm[i]))
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
    
    matrix=data.frame(time=t,s=S,e=E,i=I,r=R,inc=Inc,rt=Rt)
    form=sprintf('/simulation/SEIR_Rt_%s_%s.csv',year,j)
    form=paste0(getwd(),form)
    write.csv(matrix[2:n,],file=form)
  }
  file1=sprintf("simulation/SEIR_Rt_%s_1.csv",year)
  file2=sprintf("simulation/SEIR_Rt_%s_2.csv",year)
  file3=sprintf("simulation/SEIR_Rt_%s_3.csv",year)
  file4=sprintf("simulation/SEIR_Rt_%s_4.csv",year)

  
  SEIR_without_phsm_rain <-read.csv(file1)
  SEIR_with_phsm_rain <-read.csv(file2)
  SEIR_with_phsm_rainLb <-read.csv(file3)
  SEIR_with_phsm_rainUb <-read.csv(file4)

  
  without_phsm<-SEIR_without_phsm_rain$inc
  with_phsm_est<-SEIR_with_phsm_rain$inc
  with_phsm_est_ci<-SEIR_with_phsm_rainLb$inc
  with_phsm_est_CI<-SEIR_with_phsm_rainUb$inc

  without_phsm=without_phsm*700000
  with_phsm_est=with_phsm_est*700000
  with_phsm_est_ci=with_phsm_est_ci*700000
  with_phsm_est_CI=with_phsm_est_CI*700000
  
  

  final=sprintf("Simulation_%s.csv",year)
  result<-data.frame(df2$date[2:n],
                     df2$schoolvac[2:n],
                     df2$Rainyseason[2:n],
                     without_phsm,
                     with_phsm_est,
                     with_phsm_est_ci,
                     with_phsm_est_CI
                     )
  
  write.csv(result,file=final)
}