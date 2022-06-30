library('dplyr')
### set directory

files<-list.files('DATASET/')
data_directory=paste0(getwd(),'/DATASET/')
result_folder='sensitivity'
if (file.exists(result_folder)) {
    print("The folder already exists")
    print("This folder includes annual sensitivity")
} else {
    dir.create(result_folder)
    print('The folder created')
    print("This folder includes annual sensitivity")
}

### automatical regression by years
### put years you want 
years<-list(2014,2015,2016,2017,2018,2019)
s0_list<-c(0.999,0.899,0.799,0.699)
R_0_list<-c(2.5,3.0,3.5,4.0,4.5)

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
        model3 <- lm(log(Rt) ~ Cu_hfmd+schoolvac+Rainyseason,data =df)
        cat("currnet year: ",year,'\n')
    } else {
        model3 <-lm(log(Rt) ~ Cu_hfmd+schoolvac+Rainyseason+intervention,data=df)
        cat("currnet year: ",year,'\n')
    }
    c<-summary(model3)
    
    est<-as.data.frame(c[4])[3,1]
    est_ci<-as.data.frame(confint(model3))[3,1]
    est_CI<-as.data.frame(confint(model3))[3,2]
    reg_coef_mat_phsm=c(0.0,est,est_ci,est_CI)

    phsm=df2$schoolvac


    ### rain coefficient
    rain_est<-as.data.frame(c[4])[4,1]
    rain_est_ci<-as.data.frame(confint(model3))[4,1]
    rain_est_CI<-as.data.frame(confint(model3))[4,2]

    reg_coef_mat_rain=c(0.0, rain_est, rain_est_ci, rain_est_CI)
    rain=df2$Rainyseason

    if (year==2015){
        est2<-as.data.frame(c[4])[5,1]
        est2_ci<-as.data.frame(confint(model3))[5,1]
        est2_CI<-as.data.frame(confint(model3))[5,2]
        
        reg_coef_mat_int=c(0.0,est2,est2_ci,est2_CI)
        intervention=df2$intervention
        
    } else{
        
    }
    np=length(reg_coef_mat_phsm)
    n=length(phsm)
    
    for(ss in s0_list){
        "S0"<-ss
        "E0"<-0.00015
        "I0"<-df2$cases_total[1]/df3$Cu_hfmd[length(df3$Cu_hfmd)]
        "R0"<-1-(S0+E0+I0)
        #sigma: Rate of exposed individuals become infected
        "sigma"=0.57
        #gamma: recovery rate
        "gamma"=0.50
        
        for (R_0 in R_0_list){
            "R_0"<-R_0
            beta0=(S0*R_0/S0)*gamma
            beta=beta0
            
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
                form=sprintf('/sensitivity/SEIR_Rt_%s_%s_%s_%s.csv',year,j,ss,R_0)
                form=paste0(getwd(),form)
                write.csv(matrix[2:n,],file=form)
            }
            file1=sprintf("sensitivity/SEIR_Rt_%s_1_%s_%s.csv",year,ss,R_0)
            file2=sprintf("sensitivity/SEIR_Rt_%s_2_%s_%s.csv",year,ss,R_0)
            file3=sprintf("sensitivity/SEIR_Rt_%s_3_%s_%s.csv",year,ss,R_0)
            file4=sprintf("sensitivity/SEIR_Rt_%s_4_%s_%s.csv",year,ss,R_0)
            
            SEIR_without_phsm_rain <-read.csv(file1)
            SEIR_with_phsm_rain <-read.csv(file2)
            SEIR_with_phsm_rainLb <-read.csv(file3)
            SEIR_with_phsm_rainUb <-read.csv(file4)
            
            without_phsm<-SEIR_without_phsm_rain$real_inc
            with_phsm_est<-SEIR_with_phsm_rain$real_inc
            with_phsm_est_ci<-SEIR_with_phsm_rainLb$real_inc
            with_phsm_est_CI<-SEIR_with_phsm_rainUb$real_inc
            
            
            final=sprintf("sensitivity/Result_%s_%s_%s.csv",year,ss,R_0)
            
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
    }
}

for (year in years){
    #define dataframe 
    li<-list()
    for(ss in s0_list){
        for (R_0 in R_0_list){
            final=sprintf("sensitivity/Result_%s_%s_%s.csv",year,ss,R_0)
            df<-read.csv(final)
            
            df2<-data.frame(df$without_phsm,
                            df$with_phsm_est,
                            df$with_phsm_est_ci,
                            df$with_phsm_est_CI)
            colnames(df2)[1]<-paste("without_phsm",paste(ss,R_0,sep="_"),sep="_")
            colnames(df2)[2]<-paste("with_phsm_est",paste(ss,R_0,sep="_"),sep="_")
            colnames(df2)[3]<-paste("with_phsm_est_ci",paste(ss,R_0,sep="_"),sep="_")  
            colnames(df2)[4]<-paste("with_phsm_est_CI",paste(ss,R_0,sep="_"),sep="_")
            li<-append(li,df2)
        }
    }
    x=bind_rows(li, .id = "column_label")
    fin=sprintf("Sensitivity_%s.csv",year)
    write.csv(x,file=fin)
}



