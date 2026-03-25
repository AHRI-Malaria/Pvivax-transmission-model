library(tidyverse)
library(readr)
library(ggplot2)
library(ggprism)


post <- read_delim("data/bayesian_posterior_samples_all_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

dt <- read_csv("data/pvivax_transmission_data.csv")


agec = mean(dt$age)



post <- post %>%
  mutate(across(everything(), ~as.numeric(gsub(",", ".", .))))

logistic <- function(x) {
  1 / (1 + exp(-x))
}


mu_P=mean(post$mu_P)         
gamma_P_age  = mean(post$gamma_P_age)
mu_G=mean(post$mu_G)
gamma_G_P=mean(post$gamma_G_P)     
gamma_G_age = mean(post$gamma_G_age)
alpha = mean(post$alpha)         
beta1 = mean(post$beta_study.1.)
beta2 = mean(post$beta_study.2.)
beta_G = mean(post$beta_G)      
beta_P = mean(post$beta_P)      
beta_age = mean(post$beta_age)
beta_fever = mean(post$beta_fever)


#the total effects of age cant be looked at without fever but can make gams and parasites disappear from the model entirely

age_data = expand.grid(age=sort(unique(dt$age)),
                       study=c(1,2,3),
                       fever=c(0,1)) %>%
  mutate(study1 = as.numeric(study==1),study2 = as.numeric(study==2))

age_data_list = list()

for(i in 1:nrow(post)){

  
  age_data_list[[i]]=age_data
  
  mu_P=post$mu_P[i]         
  gamma_P_age  = post$gamma_P_age[i]
  mu_G=post$mu_G[i]
  gamma_G_P=post$gamma_G_P[i]     
  gamma_G_age = post$gamma_G_age[i]
  alpha = post$alpha[i]         
  beta1 = post$beta_study.1.[i]
  beta2 = post$beta_study.2.[i]
  beta_G = post$beta_G[i]      
  beta_P = post$beta_P[i]      
  beta_age = post$beta_age[i]
  beta_fever = post$beta_fever[i]

age_lo = function(age,fever,study1,study2){
  logistic(alpha+beta_P*mu_P+beta_G*mu_G+beta_age*age+(beta_P*gamma_P_age+beta_G*gamma_G_age+beta_G*gamma_G_P*gamma_P_age)*(age-agec) + beta_fever*fever+beta1*study1+beta2*study2)
}

age_data_list[[i]]$pred0 = age_lo(age_data_list[[i]]$age, age_data_list[[i]]$fever, age_data_list[[i]]$study1, age_data_list[[i]]$study2)

age_data_list[[i]]$pred1 = age_lo(age_data_list[[i]]$age+5, age_data_list[[i]]$fever, age_data_list[[i]]$study1, age_data_list[[i]]$study2)

age_data_list[[i]]$diff = age_data_list[[i]]$pred1-age_data_list[[i]]$pred0

age_data_list[[i]]$fever = factor(age_data_list[[i]]$fever, labels=c("No fever","Fever"))
age_data_list[[i]]$study = factor(age_data_list[[i]]$study, labels=c("Study 1","Study 2","Study 3"))
print(i)
}

age_all <- bind_rows(age_data_list, .id = "draw")
age_all$draw <- as.numeric(age_all$draw)


age_summary <- age_all %>%
  group_by(study, fever) %>%
  summarise(
    med  = median(diff),
    l95  = quantile(diff, 0.025),
    u95  = quantile(diff, 0.975),
    .groups = "drop"
  )

rm(age_data_list, age_all, age_data)


write.csv2(age_summary,"total_age.csv",row.names = F)

ggplot(data=age_summary, aes(x=study, y=med, fill=fever))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=l95, ymax=u95),
                position=position_dodge(width=0.8),
                width=0.2)+  xlab("Age")+
  theme_prism()+
  scale_fill_viridis_d(begin=0.2)+
  ylab("Difference in proportion infected mosquitos")+
  xlab("Age")+
  ggtitle("Total effects of a 5 year higher age")+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position = "bottom", axis.title.x = element_blank())

ggsave("age_total_diff.tiff", device="tiff", width=16, height=12, units="cm", dpi=1200)

ggsave("age_total_diff.png", device="png", width=16, height=12, units="cm", dpi=1200)

# the effects of asexual parasite density on the proportion infected mosquitos through gam density

p_data = expand.grid(age=sort(unique(dt$age)),
                     pd = seq(0,6,by=1),
                     study=c(1,2,3),
                     fever=c(0,1)) %>%
  mutate(study1 = as.numeric(study==1),study2 = as.numeric(study==2))


# ggplot(data= p_data, aes(x=age, y=pd,fill=diff))+
#   geom_tile()+
#   facet_grid(~study+fever )+
#   theme_prism()+
#   ylab("Asexual parasite density (log10)")+
#   xlab("Age")+
#   scale_fill_viridis_c(labels=scales::percent)+
#   ggtitle("Effects of asexual parasite density")
# 


p_data_list=list()

for (i in 1:nrow(post)){

  p_data_list[[i]] = p_data
  
  
  mu_P=post$mu_P[i]         
  gamma_P_age  = post$gamma_P_age[i]
  mu_G=post$mu_G[i]
  gamma_G_P=post$gamma_G_P[i]     
  gamma_G_age = post$gamma_G_age[i]
  alpha = post$alpha[i]         
  beta1 = post$beta_study.1.[i]
  beta2 = post$beta_study.2.[i]
  beta_G = post$beta_G[i]      
  beta_P = post$beta_P[i]      
  beta_age = post$beta_age[i]
  beta_fever = post$beta_fever[i]
  
  
par_lo = function(age,pd,fever,study1,study2){
  logistic(alpha+ beta_age*age + beta_P*mu_P + beta_G*(mu_G + gamma_G_age*(age-agec) + gamma_G_P*(pd-mu_P)) + beta_fever*fever+beta1*study1+beta2*study2)
}



p_data_list[[i]]$pred0 = par_lo(p_data_list[[i]]$age, p_data_list[[i]]$pd, p_data_list[[i]]$fever, p_data_list[[i]]$study1, p_data_list[[i]]$study2)

p_data_list[[i]]$pred1 = par_lo(p_data_list[[i]]$age, p_data_list[[i]]$pd+1, p_data_list[[i]]$fever, p_data_list[[i]]$study1, p_data_list[[i]]$study2)

p_data_list[[i]]$diff = p_data_list[[i]]$pred1-p_data_list[[i]]$pred0

p_data_list[[i]]$fever = factor(p_data_list[[i]]$fever, labels=c("No fever","Fever"))
p_data_list[[i]]$study = factor(p_data_list[[i]]$study, labels=c("Study 1","Study 2","Study 3"))

print(i)

}



p_all <- bind_rows(p_data_list, .id = "draw")

p_all = p_all %>% 
  mutate(age = factor(ifelse(age<5,1,ifelse(age>15,3,2)),labels=c("<5yrs","5-15yrs",">15yrs")),
         par=cut(pd, c(-1,2,3,4,5,6))
         )

p_summary <- p_all %>%
  group_by(age, par, fever, study) %>%
  summarise(
    med  = median(diff),
    l95  = quantile(diff, 0.025),
    u95  = quantile(diff, 0.975),
    .groups = "drop"
  ) %>%
  mutate(par = factor(par, labels=c("\U2264 100","100-1k","1k-10k","10k-100k",">100k")))

rm(p_data_list,p_all,p_data)


write.csv2(p_summary,"mediated_parasitedens.csv",row.names = F)


ggplot(data=p_summary, aes(x=par, y=med, fill=age))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=l95, ymax=u95),
                position=position_dodge(width=0.8),
                width=0.2)+  xlab("Age")+
  theme_prism()+
  scale_fill_brewer()+
  ylab("Difference in proportion infected mosquitos")+
  xlab("Asexual parasite density (current)")+
  ggtitle("Effects of a 10X higher asexual parasite density (mediated through gametocyte density)")+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position = "top")+
  facet_grid(vars(study),vars(fever))

ggsave("par_mediated_diff.tiff", device="tiff", width=27, height=18, units="cm", dpi=1200)
ggsave("par_mediated_diff.png", device="png", width=27, height=18, units="cm", dpi=1200)

# direct effects of age, gam density and fever and asexual parasite density (for age and asexual parasite density these are )

## gam density direct effect

g_data = expand.grid(age=sort(unique(dt$age)),
                     pd = seq(0,6,by=1),
                     gd = seq(0,6,by=1),
                     study=c(1,2,3),
                     fever=c(0,1)) %>%
  mutate(study1 = as.numeric(study==1),study2 = as.numeric(study==2))


g_data_list = list()

for (i in 1:nrow(post)){
  
  g_data_list[[i]] = g_data

  mu_P=post$mu_P[i]         
  gamma_P_age  = post$gamma_P_age[i]
  mu_G=post$mu_G[i]
  gamma_G_P=post$gamma_G_P[i]     
  gamma_G_age = post$gamma_G_age[i]
  alpha = post$alpha[i]         
  beta1 = post$beta_study.1.[i]
  beta2 = post$beta_study.2.[i]
  beta_G = post$beta_G[i]      
  beta_P = post$beta_P[i]      
  beta_age = post$beta_age[i]
  beta_fever = post$beta_fever[i]
  
  full_lo = function(age,pd,gd,fever,study1,study2){
    logistic(alpha+ beta_age*age + beta_P*pd + beta_G*gd + beta_fever*fever+beta1*study1+beta2*study2)
  }
  
g_data_list[[i]]$diff = full_lo(g_data_list[[i]]$age, g_data_list[[i]]$pd, g_data_list[[i]]$gd+1, g_data_list[[i]]$fever, g_data_list[[i]]$study1, g_data_list[[i]]$study2) -full_lo(g_data_list[[i]]$age, g_data_list[[i]]$pd, g_data_list[[i]]$gd, g_data_list[[i]]$fever, g_data_list[[i]]$study1, g_data_list[[i]]$study2)

g_data_list[[i]]$fever = factor(g_data_list[[i]]$fever, labels=c("No fever","Fever"))
g_data_list[[i]]$study = factor(g_data_list[[i]]$study, labels=c("Study 1","Study 2","Study 3"))

print(i)


}

g_all <- bind_rows(g_data_list, .id = "draw")

g_all = g_all %>% 
  mutate(age = factor(ifelse(age<5,1,ifelse(age>15,3,2)),labels=c("<5yrs","5-15yrs",">15yrs")),
         par=cut(pd, c(-1,2,3,4,5,6)),
         gam=cut(gd, c(-1,2,3,4,5,6))
  )

g_summary <- g_all %>%
  group_by(age, gam, par, fever, study) %>%
  summarise(
    med  = median(diff),
    l95  = quantile(diff, 0.025),
    u95  = quantile(diff, 0.975),
    .groups = "drop"
  ) %>%
  mutate(par = factor(par, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k")),gam = factor(gam, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k"))
         )

rm(g_data_list,g_all,g_data)



ggplot(data=g_summary, aes(x=gam, y=med, fill=par))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=l95, ymax=u95),
                position=position_dodge(width=0.8),
                width=0.2)+  xlab("Age")+
  theme_prism()+
  scale_fill_viridis_d(name = "Asexual\nparasite\ndensity"  ,begin=0.3, end=0.9)+
  ylab("Difference in proportion infected mosquitos")+
  xlab("Gametocyte density (current)")+
  ggtitle("Effects of a 10X higher gametocyte density")+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position = "top")+
  facet_grid(vars(study),vars(fever,age))+
  theme(legend.title = element_text())

ggsave("gam_direct_diff.tiff", device="tiff", width=60, height=30, units="cm", dpi=1200)

ggsave("gam_direct_diff.png", device="png", width=60, height=30, units="cm", dpi=1200)


## direct effect of fever


f_data = expand.grid(age=sort(unique(dt$age)),
                     pd = seq(0,6,by=1),
                     gd = seq(0,6,by=1),
                     study=c(1,2,3),
                     fever=c(0)) %>%
  mutate(study1 = as.numeric(study==1),study2 = as.numeric(study==2))


f_data_list = list()

for (i in 1:nrow(post)){
  
  f_data_list[[i]] = f_data
  
  mu_P=post$mu_P[i]         
  gamma_P_age  = post$gamma_P_age[i]
  mu_G=post$mu_G[i]
  gamma_G_P=post$gamma_G_P[i]     
  gamma_G_age = post$gamma_G_age[i]
  alpha = post$alpha[i]         
  beta1 = post$beta_study.1.[i]
  beta2 = post$beta_study.2.[i]
  beta_G = post$beta_G[i]      
  beta_P = post$beta_P[i]      
  beta_age = post$beta_age[i]
  beta_fever = post$beta_fever[i]
  
  full_lo = function(age,pd,gd,fever,study1,study2){
    logistic(alpha+ beta_age*age + beta_P*pd + beta_G*gd + beta_fever*fever+beta1*study1+beta2*study2)
  }
  
  f_data_list[[i]]$diff = full_lo(f_data_list[[i]]$age, f_data_list[[i]]$pd, f_data_list[[i]]$gd, f_data_list[[i]]$fever+1, f_data_list[[i]]$study1, f_data_list[[i]]$study2) -full_lo(f_data_list[[i]]$age, f_data_list[[i]]$pd, f_data_list[[i]]$gd, f_data_list[[i]]$fever, f_data_list[[i]]$study1, f_data_list[[i]]$study2)
  

  f_data_list[[i]]$study = factor(f_data_list[[i]]$study, labels=c("Study 1","Study 2","Study 3"))
  
  print(i)
  
  
}

f_all <- bind_rows(f_data_list, .id = "draw")

f_all = f_all %>% 
  mutate(age = factor(ifelse(age<5,1,ifelse(age>15,3,2)),labels=c("<5yrs","5-15yrs",">15yrs")),
         par=cut(pd, c(-1,2,3,4,5,6)),
         gam=cut(gd, c(-1,2,3,4,5,6))
  )

f_summary <- f_all %>%
  group_by(age, gam, par, fever, study) %>%
  summarise(
    med  = median(diff),
    l95  = quantile(diff, 0.025),
    u95  = quantile(diff, 0.975),
    .groups = "drop"
  ) %>%
  mutate(par = factor(par, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k")),gam = factor(gam, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k"))
  )

rm(f_data_list,f_all,f_data)



ggplot(data=f_summary, aes(x=gam, y=med, fill=par))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=l95, ymax=u95),
                position=position_dodge(width=0.8),
                width=0.2)+  xlab("Age")+
  theme_prism()+
  scale_fill_viridis_d(name = "Asexual\nparasite\ndensity"  ,begin=0.3, end=0.9)+
  ylab("Difference in proportion infected mosquitos")+
  xlab("Gametocyte density (current)")+
  ggtitle("Effects of differences in fever status")+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position = "top")+
  facet_grid(vars(study),vars(age))+
  theme(legend.title = element_text())

ggsave("fever_direct_diff.tiff", device="tiff", width=60, height=30, units="cm", dpi=1200)

ggsave("fever_direct_diff.png", device="png", width=60, height=30, units="cm", dpi=1200)


## residual effect of asexual parasiete densiyt

p2_data = expand.grid(age=sort(unique(dt$age)),
                     pd = seq(0,6,by=1),
                     gd = seq(0,6,by=1),
                     study=c(1,2,3),
                     fever=c(0,1)) %>%
  mutate(study1 = as.numeric(study==1),study2 = as.numeric(study==2))


p2_data_list = list()

for (i in 1:nrow(post)){
  
  p2_data_list[[i]] = p2_data
  
  mu_P=post$mu_P[i]         
  gamma_P_age  = post$gamma_P_age[i]
  mu_G=post$mu_G[i]
  gamma_G_P=post$gamma_G_P[i]     
  gamma_G_age = post$gamma_G_age[i]
  alpha = post$alpha[i]         
  beta1 = post$beta_study.1.[i]
  beta2 = post$beta_study.2.[i]
  beta_G = post$beta_G[i]      
  beta_P = post$beta_P[i]      
  beta_age = post$beta_age[i]
  beta_fever = post$beta_fever[i]
  
  full_lo = function(age,pd,gd,fever,study1,study2){
    logistic(alpha+ beta_age*age + beta_P*pd + beta_G*gd + beta_fever*fever+beta1*study1+beta2*study2)
  }
  
  p2_data_list[[i]]$diff = full_lo(p2_data_list[[i]]$age, p2_data_list[[i]]$pd+1, p2_data_list[[i]]$gd, p2_data_list[[i]]$fever, p2_data_list[[i]]$study1, p2_data_list[[i]]$study2) -full_lo(p2_data_list[[i]]$age, p2_data_list[[i]]$pd, p2_data_list[[i]]$gd, p2_data_list[[i]]$fever, p2_data_list[[i]]$study1, p2_data_list[[i]]$study2)
  
  p2_data_list[[i]]$fever = factor(p2_data_list[[i]]$fever, labels=c("No fever","Fever"))
  p2_data_list[[i]]$study = factor(p2_data_list[[i]]$study, labels=c("Study 1","Study 2","Study 3"))
  
  print(i)
  
  
}

p2_all <- bind_rows(p2_data_list, .id = "draw")

p2_all = p2_all %>% 
  mutate(age = factor(ifelse(age<5,1,ifelse(age>15,3,2)),labels=c("<5yrs","5-15yrs",">15yrs")),
         par=cut(pd, c(-1,2,3,4,5,6)),
         gam=cut(gd, c(-1,2,3,4,5,6))
  )

p2_summary <- p2_all %>%
  group_by(age, gam, par, fever, study) %>%
  summarise(
    med  = median(diff),
    l95  = quantile(diff, 0.025),
    u95  = quantile(diff, 0.975),
    .groups = "drop"
  ) %>%
  mutate(par = factor(par, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k")),gam = factor(gam, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k"))
  )

rm(p2_data_list,p2_all,p2_data)



ggplot(data=p2_summary, aes(x=gam, y=med, fill=par))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=l95, ymax=u95),
                position=position_dodge(width=0.8),
                width=0.2)+  xlab("Age")+
  theme_prism()+
  scale_fill_viridis_d(name = "Current asexual\nparasite\ndensity"  ,begin=0.3, end=0.9)+
  ylab("Difference in proportion infected mosquitos")+
  xlab("Gametocyte density (current)")+
  ggtitle("Residual ffects of a 10X higher asexual parasite density")+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position = "top")+
  facet_grid(vars(study),vars(fever,age))+
  theme(legend.title = element_text())

ggsave("par_residual_diff.tiff", device="tiff", width=60, height=30, units="cm", dpi=1200)

ggsave("par_residual_diff.png", device="png", width=60, height=30, units="cm", dpi=1200)




## residual effect of age

age2_data = expand.grid(age=sort(unique(dt$age)),
                      pd = seq(0,6,by=1),
                      gd = seq(0,6,by=1),
                      study=c(1,2,3),
                      fever=c(0,1)) %>%
  mutate(study1 = as.numeric(study==1),study2 = as.numeric(study==2))


age2_data_list = list()

for (i in 1:nrow(post)){
  
  age2_data_list[[i]] = age2_data
  
  mu_P=post$mu_P[i]         
  gamma_P_age  = post$gamma_P_age[i]
  mu_G=post$mu_G[i]
  gamma_G_P=post$gamma_G_P[i]     
  gamma_G_age = post$gamma_G_age[i]
  alpha = post$alpha[i]         
  beta1 = post$beta_study.1.[i]
  beta2 = post$beta_study.2.[i]
  beta_G = post$beta_G[i]      
  beta_P = post$beta_P[i]      
  beta_age = post$beta_age[i]
  beta_fever = post$beta_fever[i]
  
  full_lo = function(age,pd,gd,fever,study1,study2){
    logistic(alpha+ beta_age*age + beta_P*pd + beta_G*gd + beta_fever*fever+beta1*study1+beta2*study2)
  }
  
  age2_data_list[[i]]$diff = full_lo(age2_data_list[[i]]$age+5, age2_data_list[[i]]$pd, age2_data_list[[i]]$gd, age2_data_list[[i]]$fever, age2_data_list[[i]]$study1, age2_data_list[[i]]$study2) -full_lo(age2_data_list[[i]]$age, age2_data_list[[i]]$pd, age2_data_list[[i]]$gd, age2_data_list[[i]]$fever, age2_data_list[[i]]$study1, age2_data_list[[i]]$study2)
  
  age2_data_list[[i]]$fever = factor(age2_data_list[[i]]$fever, labels=c("No fever","Fever"))
  age2_data_list[[i]]$study = factor(age2_data_list[[i]]$study, labels=c("Study 1","Study 2","Study 3"))
  
  print(i)
  
  
}

age2_all <- bind_rows(age2_data_list, .id = "draw")

age2_all = age2_all %>% 
  mutate(age = factor(ifelse(age<5,1,ifelse(age>15,3,2)),labels=c("<5yrs","5-15yrs",">15yrs")),
         par=cut(pd, c(-1,2,3,4,5,6)),
         gam=cut(gd, c(-1,2,3,4,5,6))
  )

age2_summary <- age2_all %>%
  group_by(age, gam, par, fever, study) %>%
  summarise(
    med  = median(diff),
    l95  = quantile(diff, 0.025),
    u95  = quantile(diff, 0.975),
    .groups = "drop"
  ) %>%
  mutate(par = factor(par, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k")),gam = factor(gam, labels=c("≤100","100-1k","1k-10k","10k-100k",">100k"))
  )

rm(age2_data_list,age2_all,age2_data)



ggplot(data=age2_summary, aes(x=gam, y=med, fill=par))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_errorbar(aes(ymin=l95, ymax=u95),
                position=position_dodge(width=0.8),
                width=0.2)+  xlab("Age")+
  theme_prism()+
  scale_fill_viridis_d(name = "Asexual\nparasite\ndensity"  ,begin=0.3, end=0.9)+
  ylab("Difference in proportion infected mosquitos")+
  xlab("Gametocyte density (current)")+
  ggtitle("Residual effects of a 5 year higher age")+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position = "top")+
  facet_grid(vars(study),vars(fever,age))+
  theme(legend.title = element_text())

ggsave("age_residual_diff.tiff", device="tiff", width=60, height=30, units="cm", dpi=1200)

ggsave("age_residual_diff.png", device="png", width=60, height=30, units="cm", dpi=1200)

