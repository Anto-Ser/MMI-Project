## --- Stacked Histograms and bubble plots----

rm(list=setdiff(ls(), c("df_MI", "MI_info", "Number_barcode", "Bubble_plot", "Stacked_histo", "Select_and_order_Organs")))

## --- MDA-231  ----
MDA.color = read.csv("./data/Colour_MDA_Histo.csv")

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Tum", info_dataframe = Number_barcode)

MDA.tum.histo = Stacked_histo(MDA.tum, angle.x = 90, my_col = MDA.color$x)

MDA.tum.histo$Stack.plot # Exported as PDF 5x12
#MDA.tum.histo$Colour

MDA.tum.exp1= MDA.tum[,grep("Exp38_", names(MDA.tum))]
MDA.tum.exp2= MDA.tum[,grep("Exp1009_", names(MDA.tum))]
MDA.tum.exp3= MDA.tum[,grep("Exp1011_", names(MDA.tum))]

MDA.tum.histo1 = Stacked_histo(MDA.tum.exp1, angle.x = 90, my_col = MDA.color$x)
MDA.tum.histo2 = Stacked_histo(MDA.tum.exp2, angle.x = 90, my_col = MDA.color$x)
MDA.tum.histo3 = Stacked_histo(MDA.tum.exp3, angle.x = 90, my_col = MDA.color$x)

MDA.tum.histo1$Stack.plot #Exported as PDF 5x8
MDA.tum.histo2$Stack.plot #Exported as PDF 5x11
MDA.tum.histo3$Stack.plot #Exported as PDF 5x8


## --- Figure 1 c) Bubble plot for MDA tumours ----

Col.MDA.bbp = read.csv("./data/Y.Position_Colours_MDA_bubbles.csv")$Colours
Y.MDA.bbp = read.csv("./data/Y.Position_Colours_MDA_bubbles.csv")$Y.Position

MDA.tum = Select_and_order_Organs(df_MI, "MDA-231", "Tum", info_dataframe = Number_barcode)

MDA.tum = MDA.tum %>% select(c("Exp38_232_Tum","Exp38_234_Tum","Exp38_235_Tum",
                               "Exp38_236_Tum","Exp1009_362_Tum","Exp1009_363_Tum",
                               "Exp1009_364_Tum","Exp1009_365_Tum","Exp1011_830_Tum",
                               "Exp1011_831_Tum","Exp1011_834_Tum","Exp1011_835_Tum",
                               "Exp1011_836_Tum","Exp1011_837_Tum"),everything())
names(MDA.tum)[11:24]
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x20

## --- Supp Figure 2 d) Bubble plot MDA tumours by experiments ----
MDA.tum.exp1= MDA.tum[,grep("Exp38_", names(MDA.tum))]
MDA.tum.exp2= MDA.tum[,grep("Exp1009_", names(MDA.tum))]
MDA.tum.exp3= MDA.tum[,grep("Exp1011_", names(MDA.tum))]

MDA.tum.bbp_1 = Bubble_plot(MDA.tum.exp1, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp_2 = Bubble_plot(MDA.tum.exp2, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp_3 = Bubble_plot(MDA.tum.exp3, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)

MDA.tum.bbp_1$BBP #5x8
MDA.tum.bbp_2$BBP #5x11
MDA.tum.bbp_3$BBP#5x8

## --- Supp Figure 1 i) Stack histogram PDX-412 tumours ----

PDX.tum = Select_and_order_Organs(df_MI, "PDX-T412", "Tum", info_dataframe = Number_barcode)
PDX.color = read.csv("./data/Colour_PDX_Histo.csv")

names(PDX.tum)[11:20]

PDX.tum = PDX.tum %>% select(c("Exp72_585_Tum","Exp72_586_Tum","Exp72_587_Tum","Exp72_588_Tum",
                               "Exp72_589_Tum","Exp1012_847_Tum","Exp1012_848_Tum",
                               "Exp1012_849_Tum","Exp1012_850_Tum","Exp1012_851_Tum"),everything())

PDX.tum.histo = Stacked_histo(PDX.tum, angle.x = 90, my_col = PDX.color$x)

PDX.tum.histo$Stack.plot #5x12

## --- Supp Figure 7 Bubble plot MDA-231 organs and PDX-412 organs ----


MDA.tum = Select_and_order_Organs(df_MI, "MDA-231", "Lung", info_dataframe = Number_barcode)
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x20

MDA.tum = Select_and_order_Organs(df_MI, "MDA-231", "Liver", info_dataframe = Number_barcode)
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x20

MDA.tum = Select_and_order_Organs(df.MI, "MDA-231", "CTC", info_dataframe = Number_barcode)
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x17

#write.csv(data.frame(Y.Position = Y.MDA.bbp, Colours = Col.MDA.bbp ), file="Y.Position_Colours_MDA_bubbles.csv")

## --- Bubble plot for all organs grouped for PDX ----

Col.MDA.bbp = read.csv("Y.Position_Colours_PDX_bubbles.csv")$Colours
Y.MDA.bbp = read.csv("Y.Position_Colours_PDX_bubbles.csv")$Y.Position

MDA.tum = Select_and_order_Organs(df_MI, Model = "PDX-T412", Organ = "Lung", info_dataframe = Number_barcode)

MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x15

MDA.tum = Select_and_order_Organs(df_MI, "PDX-T412", "CTC", info_dataframe = Number_barcode)
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x8

#write.csv(data.frame(Y.Position = Y.MDA.bbp, Colours = Col.MDA.bbp ), file="Y.Position_Colours_PDX_bubbles.csv")

## --- Specific examples stacked histograms ----
Examples_MDA =c(359,362,371, 377) # Figure 3
Examples_PDX =c(592,588,582, 859) # Figure 3

Examples_MDA_Organs =c(356,365,371, 377,366) # Figure 4
Examples_PDX_Organs =c(843,851,582, 857) # Not included

pdf(file= "AAAAAAA.pdf", height = 5, width = 5)
for (i in Examples_MDA) {

  Ms_Tum_exp = df_MI[,grep(paste("_", i, "_", sep = ""), names(df_MI)), drop=F]

  colnames(Ms_Tum_exp) <- str_split(names(Ms_Tum_exp), "_", simplify = T)[,3]
  Ms_Tum_exp = Ms_Tum_exp[,c("Tum", "Lung", "Liver", "CTC")]

  p2 = Stacked_histo(Ms_Tum_exp, my_col = MDA.color$x)

  print(p2$Stack.plot + labs(title = paste("Mouse",i)))#4x4(MDA) #4x3(PDX)
}
dev.off()
