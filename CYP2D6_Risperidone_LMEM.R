#Script for Making the CYP2D6 variants dataset and constructing a Linear Mixed Effects Model for Risperidone
#Mahmoud Zidan         17/01/2023      Center for Psychopharmacology, Diakonhjemmet Sykhus

#Import the Risperidone dataset
library('readxl')
Risp= read_excel("Excel_file_with_risperidone_patient_data.xslx")

#Import the CYP2D6 variants dataset
CYP2D6= read_excel("Excel_file_with_risperidone_patients_imputed_genotypes_for_the_NDUFA6-DT_SNPs.xslx")

#Clean the datasets before merging them
Risp= Risp[Risp$GeneticID %in% CYP2D6$Genetic_ID,]
CYP2D6= CYP2D6[, c('Genetic_ID','rs28842514', 'rs116907259', 'rs9623522', 'rs2839709', 'rs71329131', 'rs118017929', 'rs111988926')]
rep= as.data.frame(table(Risp$GeneticID))
rep= rep[match(rep$Var1, CYP2D6$Genetic_ID),]
CYP2D6$Observations= rep$Freq
CYP2D6= data.frame(CYP2D6[rep(seq_len(dim(CYP2D6)[1]), CYP2D6$Observations), c(1:8), drop = FALSE], row.names=NULL)
CYP2D6= CYP2D6[match(Risp$GeneticID, CYP2D6$Genetic_ID),]
Risp= cbind(Risp, CYP2D6)
Risp= Risp[,-26]
Risp1= Risp[!Risp$CYP2D6_ph1== 'Unknown dup',]
save(Risp, file = "Dataframe_with_risperidone_patient_data_and_SNPs_genotypes.rda")

#Make a dataframe for NM only 
NM= Risp[Risp$CYP2D6=='*1/*1',]
LnMR= NM$lnMR2
SNP1= as.factor(NM$rs28842514)
SNP2= as.factor(NM$rs116907259)
SNP3= as.factor(NM$rs9623522)
SNP4= as.factor(NM$rs2839709)
SNP5= as.factor(NM$rs71329131)
SNP6= as.factor(NM$rs118017929)
SNP7= as.factor(NM$rs111988926)
TAD= NM$TAD_impute
ID= NM$GeneticID

#Linear Mixed Model Analysis for the 7 SNPs
library(nlme)
library(lme4)
SNP1[SNP1=='TG']= 'TT'
SNP2[SNP2=='CT']= 'CC'
SNP3[SNP3=='TC']= 'TT'
SNP4[SNP4=='AG']= 'GG'
SNP5[SNP5=='CT']= 'CC'
SNP6[SNP6=='GC']= 'GG'
SNP7[SNP7=='CG']= 'CC'

LMM= lme(LnMR ~ SNP1 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP2 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP3 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP4 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP5 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP6 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP7 + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

##################
#Plotting the results
NM$rs28842514= paste0('SNP1', NM$rs28842514)
NM$rs116907259= paste0('SNP2', NM$rs116907259)
NM$rs9623522= paste0('SNP3', NM$rs9623522)
NM$rs2839709= paste0('SNP4', NM$rs2839709)
NM$rs71329131= paste0('SNP5', NM$rs71329131)
NM$rs118017929= paste0('SNP6', NM$rs118017929)
NM$rs111988926= paste0('SNP7', NM$rs111988926)
NMPlot=data.frame(matrix(NA, ncol = 2, nrow = 5950))
colnames(NMPlot)= c('lnMR2', 'SNPs')
NMPlot$lnMR2=rep(NM$lnMR2)
NMPlot$SNPs[1:850]= NM$rs28842514
NMPlot$SNPs[851:1700]= NM$rs116907259
NMPlot$SNPs[1701:2550]= NM$rs9623522
NMPlot$SNPs[2551:3400]= NM$rs2839709
NMPlot$SNPs[3401:4250]= NM$rs71329131
NMPlot$SNPs[4251:5100]= NM$rs118017929
NMPlot$SNPs[5101:5950]= NM$rs111988926
NMPlot$SNPs[NMPlot$SNPs=='SNP1TG']= 'SNP1TT'
NMPlot$SNPs[NMPlot$SNPs== 'SNP2CT']= 'SNP2CC'
NMPlot$SNPs[NMPlot$SNPs== 'SNP6GC']= 'SNP6GG'

ggplot(NMPlot, aes(x= SNPs, 
                   y=lnMR2, fill=SNPs ))+
  geom_boxplot() +
  labs(x= 'Genotype', y='Ln (MR)', title = '')+
  theme_bw() + theme(panel.border = element_blank())+
  theme(axis.line = element_line(color = 'black')) + theme(legend.position = 'none')+
  scale_fill_manual(values=c("orange3", "orange2", "deepskyblue", 'deepskyblue3','indianred', 'indianred2', 
                                       'purple3', 'purple1','springgreen2','springgreen3',
                                      'pink','pink2','gold4','gold3'))+
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=1, alpha= 0.2)+
  geom_signif(comparisons = list(c("SNP1GG", "SNP1TT")), annotations="*",textsize = 8,
              map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("SNP3CC", "SNP3TT")), annotations="*",textsize = 8,
              map_signif_level=TRUE)

################################################################
#Linear Mixed Models for haplotypes in NM only
Haplotype= as.factor(NM$Haplotype_All_Phenotypes_1_2_3)

LMM= lme(LnMR ~ Haplotype + TAD, random = ~1|ID)
summary(LMM)

#Linear Mixed Models for haplotypes in UM only 
UM= Risp[Risp$CYP2D6_ph1=='UM',]

LnMR= UM$lnMR2 
TAD= UM$TAD_impute
Age= UM$Age
ID= UM$GeneticID
Haplotype= UM$Haplotype_All_Phenotypes_1_2_3

library(nlme)
library(lme4)

LMM= lme(LnMR ~ Haplotype + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

#Linear Mixed Models for haplotypes in IM only 
IM= Risp[Risp$CYP2D6_ph1=='IM',]

LnMR= IM$lnMR2
TAD= IM$TAD_impute
Age= IM$Age
ID= IM$GeneticID
Haplotype= IM$Haplotype_All_Phenotypes_1_2_3

library(nlme)
library(lme4)

LMM= lme(LnMR ~ Haplotype + TAD, random = ~1|ID, method = 'ML')
summary(LMM)

#linear Mixed Models for haplotypes in PM only 
PM= Risp[Risp$CYP2D6_ph1=='PM',]

LnMR= PM$lnMR2
TAD= PM$TAD_impute
Age= PM$Age
ID= PM$GeneticID
Haplotype= PM$Haplotype_All_Phenotypes_1_2_3

LMM= lme(LnMR ~ Haplotype + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

#Plot the results
ggplot(Risp, aes(x= Haplotype_All_Phenotypes_1_2_3, 
                 y=lnMR2, fill=Haplotype_All_Phenotypes_1_2_3 ))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x= 'Diplotype', y='Ln (MR)', title = '')+
  theme_bw() + theme(panel.border = element_blank())+
  theme(axis.line = element_line(color = 'black'),axis.text.x = element_text(angle = 40, vjust = 0.5, hjust=0.5)) + theme(legend.position = 'none') +
  geom_point(position = position_jitter(seed = 1, width = 0.08), size=1, alpha= 0.25)+ 
  scale_x_discrete(labels=c('IM1' = "GTC/GTC",'NM1' = "GTC/GTC",'PM1' = "GTC/GTC",'UM1' = "GTC/GTC",
                            "IM2" = "GTC/GTT","NM2" = "GTC/GTT","PM2" = "GTC/GTT","UM2" = "GTC/GTT",
                            'IM3'= 'GTC/GCC','NM3'= 'GTC/GCC','PM3'= 'GTC/GCC','UM3'= 'GTC/GCC',
                            'IM4'= 'GTC/TTC', 'NM4'= 'GTC/TTC', 'NM5'= 'GTT/GTT', 'NM6'= 'GTT/GCC',
                            'UM6'= 'GTT/GCC', 'NM7'= 'GTT/TTC')) +
  scale_fill_manual(values=c("white", "indianred2", "white", 'royalblue2','white', 'indianred2', 
                                    'white', 'royalblue2','white','white',
                                    'white','white','indianred2','white', 'white','indianred2',
                                    'white','white'))+
                                      geom_signif(comparisons = list(c("NM1", "NM2")), annotations="*",textsize = 7,
                                                  map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("NM1", "NM4")), annotations="**",textsize = 7,y= 8.5,
              map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("IM1", "IM2")), annotations="***",textsize = 7,
              map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("UM1", "UM6")), annotations="*",textsize = 7,
              map_signif_level=TRUE)

###############################################
#Linear Mixed Models for the TTC haplotype in NM
LnMR= NM$lnMR2
TAD= NM$TAD_impute
ID= NM$GeneticID
TTC= NM$TTC_Presence

LMM= lme(LnMR ~ TTC + TAD, random = ~1|ID)
summary(LMM)

ggplot(NM, aes(x= TTC_Presence, 
                 y=lnMR2, fill=TTC_Presence ))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x= 'Presecnce of the TTC Haplotype', y='Ln (MR)', title = '')+
  theme_bw() + theme(panel.border = element_blank(),text = element_text(size=20))+
  theme(axis.line = element_line(color = 'black')) + theme(legend.position = 'none') +
  scale_x_discrete(labels=c( 'NMTRUE'= 'TRUE','NMFALSE'= 'FALSE')) +
  scale_fill_manual(values=c("white", 'royalblue2'))+
  geom_point(position = position_jitter(seed = 1, width = 0.2), size=1, alpha= 0.2)+
  geom_signif(comparisons = list(c("NMTRUE", "NMFALSE")), annotations="*",textsize = 9,
              y= 7.7, map_signif_level=TRUE)


###############################################
#Linear Mixed Models for the TTC haplotype in NM
LnMR= NM$lnMR2
TAD= NM$TAD_impute
ID= NM$GeneticID
GTT= NM$GTT_Presence

LMM= lme(LnMR ~ TTC + TAD, random = ~1|ID)
summary(LMM)

ggplot(NM, aes(x= GTT_Presence, 
               y=lnMR2, fill=GTT_Presence ))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x= 'Presecnce of the GTT Haplotype', y='Ln (MR)', title = '')+
  theme_bw() + theme(panel.border = element_blank(),text = element_text(size=20))+
  theme(axis.line = element_line(color = 'black')) + theme(legend.position = 'none') +
  scale_x_discrete(labels=c( 'NMTRUE'= 'TRUE','NMFALSE'= 'FALSE')) +
  scale_fill_manual(values=c("white", 'indianred2'))+
  geom_point(position = position_jitter(seed = 1, width = 0.2), size=1, alpha= 0.2)+
  geom_signif(comparisons = list(c("NMTRUE", "NMFALSE")), annotations="*",textsize = 9,
              y= 7.7, map_signif_level=TRUE)



