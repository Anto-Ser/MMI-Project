
rm(list=setdiff(ls(), c("df_MI", "MI_info", "Bubble_plot", "Stacked_histo", "Select_and_order_Organs")))

## --- Number of barcode all organs csv export  ----

Number_barcode = data.frame(nam=names(df_MI), bc =colSums(df_MI>0), shannon= diversity(t(df_MI)))
Number_barcode = cbind(Number_barcode, str_split(Number_barcode$nam, "_", simplify = T))
names(Number_barcode) = c("Full_name","bc","Shannon", "Exp","Mouse","Org")

Number_barcode = Number_barcode %>%
  select(Full_name, Mouse,Org,bc, Shannon) %>%
  filter(Org %in% c("Tum","CTC","Liver","Lung","Liver1","Lung1","Liver2","SNLiver","SNLung","Blood", "Ov","Heart","PE")) %>%
  as.data.frame()

Number_barcode= left_join(Number_barcode, MI_info, by="Mouse")

Number_barcode$MI = factor(Number_barcode$MI, levels = c("IMFP","IMFPc","SC","ID","IV"))

Number_barcode[is.na(Number_barcode)] <- 0

#write.csv(Number_barcode, file="Nb_barcode_ALL_MI.csv")

## --- Figure 1 barplot with number of barcode per model: ----

Number_barcode = Number_barcode %>% group_by(Org, Cells, MI) %>% dplyr::mutate(mean_bc = mean(bc),
                                                                       sd_bc = sd(bc),
                                                                       mean_shannon = mean(Shannon),
                                                                       sd_shannon = sd(Shannon))

my_comparisons = list(c("IMFPc", "ID"), c("IMFP", "ID"), c("SC", "ID"), c("IMFP", "SC"), c("IMFP", "IMFPc"), c("IMFPc", "SC"),
                      c("IV", "IMFP"), c("IV", "IMFPc"), c("IV", "SC"), c("IV", "ID"))
my_comparisons = list(c("IMFPc", "ID"), c("IMFP", "ID"), c("SC", "ID"), c("IMFP", "SC"), c("IMFP", "IMFPc"), c("IMFPc", "SC"))


## --- Figure 1 barplot with Fold change number of barcode per model: ----

Organ = "Tum"
Model = "PDX-T412"

Mean_bc_MDA_IMFP = Number_barcode %>% filter(Org == Organ & Cells == Model) %>%
  group_by(MI) %>%
  dplyr::summarise(mean_bc = mean(bc)) %>%
  filter(MI=="IMFP") %>% pull(mean_bc)

Number_barcode_fig = Number_barcode %>%
  filter(Org == Organ & Cells == Model) %>%
  dplyr::mutate(norm_bc = bc/Mean_bc_MDA_IMFP)

Number_barcode_fig = Number_barcode_fig %>% group_by(Org, Cells, MI) %>%
  dplyr::mutate(fold = log2(norm_bc), mean_fold = mean(fold), sd_fold = sd(fold))

Number_barcode_fig %>% group_by(Org, Cells, MI) %>%
  dplyr::summarise(log2 = mean(fold), fold = mean(norm_bc))


# Log2 fold change:
Number_barcode_fig %>% filter(Org == Organ) %>%
  ggplot(aes(x=MI,y=fold))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_fold-sd_fold, ymax=mean_fold+sd_fold), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Log2 fold change", subtitle = paste(Model,Organ), caption = "t.test")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 5.3, method = "anova")+
  theme_classic2()+
  facet_wrap(~Cells)+ # 5x4 landscapre
  scale_y_continuous(limits = c(-3,10), breaks = c(-3,-2,-1,0,1,2,3))
  #scale_x_discrete(limits = rev)+
  #coord_flip()


# Raw number:
Number_barcode_fig %>% filter(Org == Organ) %>%
  ggplot(aes(x=MI,y=bc))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_bc-sd_bc, ymax=mean_bc+sd_bc), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Log2 fold change", subtitle = paste(Model,Organ), caption = "wilcox.test")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 5.3, method = "anova")+
  theme_classic2()+
  facet_wrap(~Cells) # 5x4 landscapre
  #scale_y_continuous(limits = c(-3,10), breaks = c(-3,-2,-1,0,1,2,3))


## ANOVA and Tukey follow test: ----
model <- aov(bc~MI, data=Number_barcode_fig)
summary(model)
TukeyHSD(model, conf.level=.95)



## --- PDX Organs number of barcode (including 0 (IV and ID mice)) ----

MI = c("IMFPc", "IMFP", "SC", "ID", "IV")
com = expand.grid(MI,MI)
com = as.data.frame(unique(t(apply(com, 1, sort))))
com = as.matrix(com[!com$V1 == com$V2,])
my_comparisons = split(com, seq(nrow(com)))


Organ = "CTC"
PDX_Ref = Number_barcode %>% filter(Cells == "PDX-T412" & Org == "Tum") %>% select(Exp, Mouse, Cells, MI, Org)
PDX_bc = Number_barcode %>% filter(Cells == "PDX-T412" & Org == Organ) %>% select(Exp, Mouse, Cells, MI, Org, bc, Shannon)
info_pdx_iv = MI_info %>% filter(Cells == "PDX-T412", Project == "MI", MI == "IV") %>% select(Exp, Mouse, Cells, MI)
PDX_bc = full_join(PDX_bc, info_pdx_iv)
PDX_bc = full_join(PDX_bc, PDX_Ref[,1:4])
PDX_bc$Org = Organ
PDX_bc[is.na(PDX_bc)] <- 0

PDX_bc = PDX_bc %>% group_by(Org, Cells, MI) %>% dplyr::mutate(mean_bc = mean(bc),
                                                                               sd_bc = sd(bc),
                                                                               mean_shannon = mean(Shannon),
                                                                               sd_shannon = sd(Shannon))
PDX_bc$MI = factor(PDX_bc$MI, levels = c("IMFP","IMFPc","SC","ID","IV"))

ggplot(PDX_bc,aes(x=MI,y=bc))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_bc-sd_bc, ymax=mean_bc+sd_bc), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Number of barcodes", subtitle = paste("PDX", Organ, "number of barcodes"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 10, method = "anova")+
  theme_classic()+
  facet_wrap(~Org, scales = "free") # 4x5

ggplot(PDX_bc,aes(x=MI,y=Shannon))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_shannon-sd_shannon, ymax=mean_shannon+sd_shannon), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Number of barcodes", subtitle = paste("PDX", Organ, "number of barcodes"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 4, method = "anova")+
  theme_classic()+
  facet_wrap(~Org, scales = "free") #4x10


model <- aov(bc~MI, data=PDX_bc)
summary(model)
TukeyHSD(model, conf.level=.95)

## --- MDA Organs number of barcode  ----

Number_barcode %>% filter(Cells == "MDA-231" & Org %in% c("CTC", "Lung", "Liver")) %>%
  ggplot(aes(x=MI,y=Shannon))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_shannon-sd_shannon, ymax=mean_shannon+sd_shannon), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Number of barcodes", subtitle = paste(Model, Organ, "number of barcodes"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 6, method = "anova")+
  theme_classic()+
  facet_wrap(~Org) #4x10



Number_barcode %>% filter(Cells == "MDA-231" & Org %in% c("Tum", "Lung")) %>%
  ggplot(aes(x=MI,y=bc))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_bc-sd_bc, ymax=mean_bc+sd_bc), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Number of barcodes", subtitle = paste(Model, Organ, "number of barcodes"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 6, method = "anova")+
  theme_classic()+
  facet_wrap(~Org) #4x10


## Tumour + IV lung comparison :

Number_barcode %>% filter(Cells == "MDA-231" & Org %in% c("Lung") & MI == "IV") %>%
  full_join(Number_barcode %>% filter(Cells == "MDA-231" & Org %in% c("Tum")))%>%
  dplyr::mutate(Org_2 = "Tum2") %>%
  ggplot(aes(x=MI,y=Shannon))+
  geom_bar(aes(fill=MI), stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=mean_shannon-sd_shannon, ymax=mean_shannon+sd_shannon), width=0.1)+
  geom_quasirandom(aes(shape=Exp), size=2, width = 0.1)+
  labs(x="", y="Number of barcodes", subtitle = paste(Model, Organ, "number of barcodes"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = FALSE)+
  stat_compare_means(label.y = 6, method = "anova")+
  theme_classic()+
  facet_wrap(~Org_2) #4x5

lung_tum = Number_barcode %>% filter(Cells == "MDA-231" & Org %in% c("Lung") & MI == "IV") %>%
  full_join(Number_barcode %>% filter(Cells == "MDA-231" & Org %in% c("Tum")))%>%
  dplyr::mutate(Org_2 = "Tum2")

model <- aov(Shannon~MI, data=lung_tum)
summary(model)
TukeyHSD(model, conf.level=.95)
