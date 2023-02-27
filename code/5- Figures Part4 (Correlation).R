
rm(list=setdiff(ls(), c("df_MI", "MI_info", "Number_barcode", "df_scatter", "Bubble_plot", "Stacked_histo", "Select_and_order_Organs")))

## --- Correlation MDA ----

df_cor = filter(df_scatter, Cells == "MDA-231")[,c(2:6,11)]

## by mouse
cor_list_temp = list()
for (i in unique(df_cor$Mouse)) { # Create list of matrices of correlation for each MI
  cor_list_temp[[i]] = as.data.frame(cor(log(filter(df_cor, Mouse == i)[1:4]+1), use="pairwise.complete.obs"))
}

correlation_df = list()
for (k in colnames(cor_list_temp[["358"]])) { # Create list of matrices of correlation for each organs by MI
  correlation_df[[k]] = as.data.frame(cbind(t(do.call(cbind,lapply(cor_list_temp, function(x) select(x,all_of(k))))), Mouse = names(cor_list_temp)))
}


p1 = left_join(correlation_df[["Tum"]], MI_info[,c(2,6)], by="Mouse") %>% 
  pivot_longer(!c("Mouse", "MI")) %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFP","IMFPc","SC","ID","IV")),
         name = factor(name, levels = c("CTC", "Liver", "Lung", "Tum"))) %>%
  #filter(MI != "IV") %>%
  filter(name != "Tum") %>%
  group_by(MI, name) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), 
                   sd = sd(value, na.rm = TRUE)) %>%
  ggplot(aes(x=MI, y=name, fill=mean))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=paste0(round(mean, digits = 2), "\n+/-(", round(sd, digits = 2), ")")), size=3, color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+ labs(x="", y="", subtitle = "Correlation with Tumor")+
  theme(panel.background = element_blank())


p2 = left_join(correlation_df[["Lung"]], MI_info[,c(2,6)], by="Mouse") %>% 
  pivot_longer(!c("Mouse", "MI")) %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFP","IMFPc","SC","ID","IV")),
         name = factor(name, levels = c("CTC", "Liver", "Lung", "Tum"))) %>%
  filter(name != c("Tum")) %>%
  group_by(MI, name) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), 
                   sd = sd(value, na.rm = TRUE)) %>%
  ggplot(aes(x=MI, y=name, fill=mean))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=paste0(round(mean, digits = 2), "\n+/-(", round(sd, digits = 2), ")")), size=3 , color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+ labs(x="", y="", subtitle = "Correlation with Lung")+
  theme(panel.background = element_blank())

p3 = left_join(correlation_df[["Liver"]], MI_info[,c(2,6)], by="Mouse") %>% 
  pivot_longer(!c("Mouse", "MI")) %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFP","IMFPc","SC","ID","IV")),
         name = factor(name, levels = c("CTC", "Lung", "Liver", "Tum"))) %>%
  filter(name %in% c("Lung", "CTC", "Liver")) %>%
  group_by(MI, name) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), 
                   sd = sd(value, na.rm = TRUE)) %>%
  ggplot(aes(x=MI, y=name, fill=mean))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=paste0(round(mean, digits = 2), "\n+/-(", round(sd, digits = 2), ")")), size=3, color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+ labs(x="", y="", subtitle = "Correlation with Liver")+
  theme(panel.background = element_blank())

plots = ggarrange(p1, p2, p3, ncol = 1, labels = "Log transformed (All)")

plots


## Correlation in lung with Top 10 barcodes in liver:
set.seed(44)

## by mouse with top 95% OR top 10 barcodes:
cor_list_temp = list()
for (i in unique(df_cor$Mouse)) { # Create list of matrices of correlation for each MI
  df_cor_top95 = df_cor %>% filter(Mouse == i) %>%
    arrange(desc(Liver)) %>%
    slice(1:10)
  #dplyr::mutate(lag_cumul_frq = lag(cumsum(Liver), default = 0)) %>%
  #filter(lag_cumul_frq<=950000)
  cor_list_temp[[i]] = as.data.frame(cor(log(df_cor_top95[1:4]+1), use="pairwise.complete.obs"))
}

correlation_df = list()
for (k in colnames(cor_list_temp[["358"]])) { # Create list of matrices of correlation for each organs by MI
  correlation_df[[k]] = as.data.frame(cbind(t(do.call(cbind,lapply(cor_list_temp, function(x) select(x,all_of(k))))), Mouse = names(cor_list_temp)))
}

left_join(correlation_df[["Liver"]], MI_info[,c(2,6)], by="Mouse") %>% 
  pivot_longer(!c("Mouse", "MI"))%>%
  #filter(MI =="IV") %>%
  mutate(value = as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFP","IMFPc","SC","ID","IV")))%>%
  filter(name =="Lung")%>%
  ggplot()+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey80")+
  geom_point(aes(y=value, x=MI, colour=Mouse))+
  geom_boxplot(aes(y=value, x=MI))+
  scale_y_continuous(limits = c(-1,1))+
  stat_summary(aes(y=value, x=MI, colour=Mouse))+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  labs(x="", y="Pearson correlation", subtitle = "Top 10 barcodes in liver correlation with lung")+
  theme_bw()+
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) # 5x6 Landscape
  
  group_by(name) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), sd= sd(value, na.rm = TRUE))


## single dataset (df_scatter):
cor_list_temp = list()
for (i in unique(df_cor$MI)) { # Create list of matrices of correlation for each MI
  cor_list_temp[[i]] = as.data.frame(cor(log(filter(df_cor, MI == i)[1:4]+1), use="pairwise.complete.obs"))
}
correlation_df = list()
for (k in colnames(cor_list_temp[["IMFP"]])) { # Create list of matrices of correlation for each organs by MI
  correlation_df[[k]] = as.data.frame(cbind(t(do.call(cbind,lapply(cor_list_temp, function(x) select(x,all_of(k))))), MI = names(cor_list_temp)))
}

# MDA Tumour correlation:
p1 = melt(correlation_df[["Tum"]], id.vars = "MI") %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFPc","IMFP","SC","ID","IV")),
         variable = factor(variable, levels = c("CTC", "Liver", "Lung", "Tum"))) %>%
  #filter(MI != "IV") %>%
  filter(variable != "Tum") %>%
  ggplot(aes(x=MI, y=variable, fill=value))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=round(value, digits = 2)),size=3, color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+ labs(x="", y="", subtitle = "Correlation with Tumor")+
  theme(panel.background = element_blank())

# MDA Lung correlation:
p2 = melt(correlation_df[["Lung"]], id.vars = "MI") %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFPc","IMFP","SC","ID","IV")),
         variable = factor(variable, levels = c("CTC", "Liver", "Lung", "Tum"))) %>%
  filter(variable != "Tum") %>%
  ggplot(aes(x=MI, y=variable, fill=value))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=round(value, digits = 2)),size=3, color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+labs(x="", y="", subtitle = "Correlation with Lung")+
  theme(panel.background = element_blank())
  
# MDA Liver correlation:
p3 = melt(correlation_df[["Liver"]], id.vars = "MI") %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFPc","IMFP","SC","ID","IV")),
         variable = factor(variable, levels = c("CTC", "Lung", "Liver"))) %>%
  filter(variable != "Tum") %>%
  ggplot(aes(x=MI, y=variable, fill=value))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=round(value, digits = 2)),size=3, color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+labs(x="", y="", subtitle = "Correlation with Liver")+
  theme(panel.background = element_blank())

plots_log_all = ggarrange(p1, p2, p3, ncol = 1, labels = "Log transformed (All)")

ggarrange(plots_no_log_all, plots_log_all, plots_no_log, plots_log, nrow = 1)


## --- Correlation PDX ----

df_cor = filter(df_scatter, Cells == "PDX-T412")[,c(2:6,11)]

## by mouse
cor_list_temp = list()
for (i in unique(df_cor$Mouse)) { # Create list of matrices of correlation for each MI
  cor_list_temp[[i]] = as.data.frame(cor(log(filter(df_cor, Mouse == i)[1:4]+1), use="pairwise.complete.obs"))
}

correlation_df = list()
for (k in colnames(cor_list_temp[["592"]])) { # Create list of matrices of correlation for each organs by MI
  correlation_df[[k]] = as.data.frame(cbind(t(do.call(cbind,lapply(cor_list_temp, function(x) select(x,all_of(k))))), Mouse = names(cor_list_temp)))
}


p1 = left_join(correlation_df[["Tum"]], MI_info[,c(2,6)], by="Mouse") %>% 
  pivot_longer(!c("Mouse", "MI")) %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFP","IMFPc","SC","ID","IV")),
         name = factor(name, levels = c("CTC", "Lung", "Tum"))) %>%
  #filter(MI != "IV") %>%
  filter(name != "Tum") %>%
  group_by(MI, name) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), 
                   sd = sd(value, na.rm = TRUE)) %>%
  ggplot(aes(x=MI, y=name, fill=mean))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=paste0(round(mean, digits = 2), "\n+/-(", round(sd, digits = 2), ")")), size=3, color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+ labs(x="", y="", subtitle = "Correlation with Tumor")+
  theme(panel.background = element_blank())


p2 = left_join(correlation_df[["Lung"]], MI_info[,c(2,6)], by="Mouse") %>% 
  pivot_longer(!c("Mouse", "MI")) %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFP","IMFPc","SC","ID","IV")),
         name = factor(name, levels = c("CTC","Lung", "Tum"))) %>%
  filter(name != c("Tum")) %>%
  group_by(MI, name) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), 
                   sd = sd(value, na.rm = TRUE)) %>%
  ggplot(aes(x=MI, y=name, fill=mean))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=paste0(round(mean, digits = 2), "\n+/-(", round(sd, digits = 2), ")")), size=3 , color = "white")+
  scale_fill_viridis()+ #limits= c(0,1)
  coord_fixed()+ labs(x="", y="", subtitle = "Correlation with Lung")+
  theme(panel.background = element_blank())

plots = ggarrange(p1, p2, ncol = 1, labels = "Log transformed (All)")

plots







# PDX Tumours correlation:
melt(correlation_df[["Tum"]], id.vars = "MI") %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFPc","IMFP","SC","ID","IV")),
         variable = factor(variable, levels = c("CTC", "Liver", "Lung", "Tum"))) %>%
  filter(MI != "IV" & variable != "Liver") %>%
  ggplot(aes(x=MI, y=variable, fill=value))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=round(value, digits = 2)), color = "white")+
  scale_fill_viridis(limits= c(0,1))+
  coord_fixed()+ labs(x="", y="")+
  theme(panel.background = element_blank())

# PDX Lung correlation:
melt(correlation_df[["Lung"]], id.vars = "MI") %>%
  mutate(value =  as.numeric(as.character(value)),
         MI = factor(MI, levels = c("IMFPc","IMFP","SC","ID","IV")),
         variable = factor(variable, levels = c("CTC", "Tum", "Lung", "Liver"))) %>%
  filter(!variable %in% c("Liver")) %>%
  ggplot(aes(x=MI, y=variable, fill=value))+
  geom_tile(color = "black", lwd=0.5)+
  geom_text(aes(label=round(value, digits = 2)), color = "white")+
  scale_fill_viridis(limits= c(0,1))+
  coord_fixed()+labs(x="", y="")+
  theme(panel.background = element_blank())
