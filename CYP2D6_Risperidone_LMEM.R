#Script for Making the CYP2D6 variants dataset and constructing a Linear Mixed Effects Model
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
Age= NM$Age
NFIB= as.factor(NM$NFIB)
CYP3A4= as.factor(NM$CYP3A4)
ID= NM$GeneticID

#Linear Mixed Model Analysis for the first 3 SNPs
library(nlme)
library(lme4)
SNP3[SNP3=='TT']= 'TC'
SNP4[SNP4=='GA']= 'GG'
SNP5[SNP5=='CC']= 'CT'
SNP7[SNP7=='CC']= 'CG'

LMM= lme(LnMR ~ SNP1 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP2 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP3 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

##################
#Plotting the results
LnMR_SNP1_T= NM$LnMR_Age_TAD_Corrected[NM$rs28842514== 'TG']
LnMR_SNP1_G= NM$LnMR_Age_TAD_Corrected[NM$rs28842514== 'GG']

LnMR_SNP2_T= NM$LnMR_Age_TAD_Corrected[NM$rs116907259== 'TT']
LnMR_SNP2_C= NM$LnMR_Age_TAD_Corrected[NM$rs116907259== 'CT']

LnMR_SNP3_T= NM$LnMR_Age_TAD_Corrected[NM$rs9623522== 'TC'| NM$rs9623522== 'TT']
LnMR_SNP3_C= NM$LnMR_Age_TAD_Corrected[NM$rs9623522== 'CC']

boxplot(LnMR_SNP1_G,LnMR_SNP1_T,LnMR_SNP2_T,LnMR_SNP2_C,LnMR_SNP3_C,LnMR_SNP3_T, xlab= 'SNP Genotype',
        names= c('SNP1G/G','SNP1T/G', 'SNP2T/T', 'SNP2C/T', 'SNP3C/C','SNP3T/C'),
        main='', 
        ylab= 'Ln(MR)', col= c('green3', 'green', 'hotpink', 'pink', 'deepskyblue', 'deepskyblue3'))
abline(0,0)

#Linear Mixed Model Analysis for the other 4 SNPs
LMM= lme(LnMR ~ SNP4 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP5 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP6 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

LMM= lme(LnMR ~ SNP7 + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)


#Plotting the results
LnMR_SNP4_A= NM$LnMR_Age_TAD_Corrected[NM$rs2839709== 'AA']
LnMR_SNP4_G= NM$LnMR_Age_TAD_Corrected[NM$rs2839709== 'GG'| NM$rs2839709== 'GA']

LnMR_SNP5_T= NM$LnMR_Age_TAD_Corrected[NM$rs71329131== 'TT']
LnMR_SNP5_C= NM$LnMR_Age_TAD_Corrected[NM$rs71329131== 'CC'| NM$rs71329131== 'CT']

LnMR_SNP6_C= NM$LnMR_Age_TAD_Corrected[NM$rs118017929== 'CC']
LnMR_SNP6_G= NM$LnMR_Age_TAD_Corrected[NM$rs118017929== 'GC']

LnMR_SNP7_G= NM$LnMR_Age_TAD_Corrected[NM$rs111988926== 'GG']
LnMR_SNP7_C= NM$LnMR_Age_TAD_Corrected[NM$rs111988926== 'CG'| NM$rs111988926== 'CC']

boxplot(LnMR_SNP4_A,LnMR_SNP4_G,LnMR_SNP5_T,LnMR_SNP5_C,LnMR_SNP6_C,LnMR_SNP6_G,LnMR_SNP7_G
        ,LnMR_SNP7_C, xlab= 'Allele',names= c('SNP4A','SNP4G', 'SNP5T',
                                                    'SNP5C', 'SNP6C','SNP6G',
                                                    'SNP7G', 'SNP7C'),
        main='', 
        ylab= 'Ln(MR)', col= c('royalblue1', 'royalblue3','salmon2', 'salmon3', 'springgreen',
                                       'springgreen3', 'purple1', 'purple3' ))
abline(1,0)

################################################################
#Linear Mixed Effect Models for Haplotypes of the first 3 SNPs in NM only

Haplotype= as.factor(NM$Diplotype_groups)
CYP2D6_1_2= as.factor(NM$CYP2D6_1_2)
CYP2D6_1_3= as.factor(NM$CYP2D6_1_3)
CYP2D6_2_3= as.factor(NM$CYP2D6_2_3)

#All 3 SNPs
LMM= lme(LnMR ~ Haplotype + Age+ TAD, random = ~1|ID)
summary(LMM)

plot(NM$Diplotype_groups, NM$LnMR_TAD_Age_Corrected, xlab= 'Diplotype', 
     ylab= 'Ln(MR)', main= '',
     col= c('yellow2','hotpink', 'green1', 'purple3', 'deepskyblue1','red3','chocolate1'))
abline(0.0188,0)

vioplot(NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==1], 
        NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==2],
        NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==3],
        NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==4],
        NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==5],
        NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==6],
        NM$LnMR_TAD_Age_Corrected[NM$Diplotype_groups==7],xlab= 'Diplotype', 
        ylab= 'Ln(MR)', main= '',
        col= c('yellow2','hotpink', 'green1', 'purple3', 'deepskyblue1','red3','chocolate1'))
abline(0.0188,0)

#SNP 1 & 2
LMM= lme(LnMR ~ CYP2D6_1_2 + Age+ TAD, random = ~1|ID)
summary(LMM)

plot(NM$CYP2D6_1_2, NM$LnMR_TAD_Age_Corrected, xlab= 'Diplotype', 
     ylab= 'Ln(MR)', main= '',
     col= c('yellow2','hotpink', 'green1'))
abline(0.1156,0)

vioplot(NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_2==1], 
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_2==2],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_2==3],xlab= 'Diplotype', ylab= 'Ln(MR)', 
        main= '',col= c('yellow2','hotpink', 'green1'))
abline(0.1156,0)

#SNP 1 & 3
LMM= lme(LnMR ~ CYP2D6_1_3 + Age+ TAD, random = ~1|ID)
summary(LMM)

plot(NM$CYP2D6_1_3, NM$LnMR_TAD_Age_Corrected, xlab= 'Diplotype', 
     ylab= 'Ln(MR)', main= '',
     col= c('yellow2','hotpink', 'green1', 'royalblue1', 'purple3'))
abline(-0.004224,0)

vioplot(NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_3==1], 
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_3==2],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_3==3],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_3==4],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_1_3==5],
        xlab= 'Diplotype', ylab= 'Ln(MR)', 
        main= '',col= c('yellow2','hotpink', 'green1', 'royalblue1', 'purple3'))
abline(-0.004224,0)

#SNP 2 & 3
LMM= lme(LnMR ~ CYP2D6_2_3 + Age+ TAD, random = ~1|ID)
summary(LMM)

plot(NM$CYP2D6_2_3, NM$LnMR_TAD_Age_Corrected, xlab= 'Diplotype', 
     ylab= 'Ln(MR)', main= '',
     col= c('deepskyblue1','red3','yellow2','hotpink', 'chocolate1'))
abline(-0.0216,0)

vioplot(NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_2_3==1], 
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_2_3==2],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_2_3==3],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_2_3==4],
        NM$LnMR_TAD_Age_Corrected[NM$CYP2D6_2_3==5],
        xlab= 'Diplotype', ylab= 'Ln(MR)', 
        main= '',col= c('deepskyblue1','red3','yellow2','hotpink', 'chocolate1'))
abline(-0.0216,0)

#########################################################################
# Linear Mixed Effects Models for the first 3 SNps in All Phenotypes

Haplotype= Risp$Haplotype_All_Phenotypes_1_2_3
LnMR= Risp$lnMR2
age= Risp$Age
TAD= Risp$TAD_impute
ID= Risp$GeneticID
CYP= Risp$CYP2D6_short2
TTC= Risp$TTC_Presence
LMM= lme(LnMR ~ Haplotype + age+ TAD + CYP, random = ~1|ID)
summary(LMM)

boxplot(Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==1], 
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==2],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==3],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==4],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==5],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==6],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==7],xlab= 'Diplotype', 
        ylab= 'Ln(MR)', main= '',
        col= c('purple', 'green1', 'deepskyblue1','red','yellow1', 'salmon2', 'violet'))
abline(-0.1205,0)

vioplot(Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==1], 
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==2],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==3],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==4],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==5],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==6],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==7],xlab= 'Diplotype', 
        ylab= 'Ln(MR)', main= '',
        col= c('purple', 'green1', 'deepskyblue1','red','yellow1', 'salmon2', 'violet'))
abline(-0.1205,0)
##########################

#linear Mixed Models for haplotypes in UM only 
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

Haplotype1= UM$LnMR_Age_TAD_Corrected[UM$Haplotype_All_Phenotypes_1_2_3==1]
Haplotype2= UM$LnMR_Age_TAD_Corrected[UM$Haplotype_All_Phenotypes_1_2_3==2]

Haplotype3= UM$LnMR_Age_TAD_Corrected[UM$Haplotype_All_Phenotypes_1_2_3==3]
Haplotype6= UM$LnMR_Age_TAD_Corrected[UM$Haplotype_All_Phenotypes_1_2_3==6]

boxplot(Haplotype1,Haplotype2,Haplotype3,Haplotype6,xlab= 'Allele',
        main='', 
        ylab= 'Ln(MR)', col= c('red', 'aquamarine', 'purple', 'deepskyblue3', 'violet'
                                           , 'green3','yellow', 'yellow3', 'aquamarine3'))
abline(2.0166,0)

vioplot(Haplotype1,Haplotype2,Haplotype3,Haplotype6,xlab= 'Allele',
        main='', 
        ylab= 'Ln(MR)', col= c('red', 'aquamarine', 'purple', 'deepskyblue3', 'violet'
                                    , 'green3','yellow', 'yellow3', 'aquamarine3'))
abline(2.0166,0)

##################################################
#linear Mixed Models for haplotypes in IM only 


IM= Risp[Risp$CYP2D6_ph1=='IM',]

LnMR= IM$lnMR2
TAD= IM$TAD_impute
Age= IM$Age
ID= IM$GeneticID
Haplotype= IM$Haplotype_All_Phenotypes_1_2_3

library(nlme)
library(lme4)

LMM= lme(LnMR ~ Haplotype + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

Haplotype1= IM$LnMR_Age_TAD_Corrected[IM$Haplotype_All_Phenotypes_1_2_3==1]
Haplotype2= IM$LnMR_Age_TAD_Corrected[IM$Haplotype_All_Phenotypes_1_2_3==2]

Haplotype3= IM$LnMR_Age_TAD_Corrected[IM$Haplotype_All_Phenotypes_1_2_3==3]
Haplotype4= IM$LnMR_Age_TAD_Corrected[IM$Haplotype_All_Phenotypes_1_2_3==4]

boxplot(Haplotype1,Haplotype2,Haplotype3,Haplotype4,xlab= 'Diplotype', main='', 
        ylab= 'Ln(MR)', col= c('deepskyblue', 'violet', 'green3', 'red', 'deepskyblue3','purple', 'violet', 'green3',
                                            'green2', 'indianred1','yellow3', 'yellow',
                                            'aquamarine3', 'orange'))
abline(-0.7784,0)

vioplot(Haplotype1,Haplotype2,Haplotype3,Haplotype4,xlab= 'Diplotype', main='', 
        ylab= 'Ln(MR)', col= c('deepskyblue', 'violet', 'green3', 'red', 'deepskyblue3','purple', 'violet', 'green3',
                                            'green2', 'indianred1','yellow3', 'yellow',
                                            'aquamarine3', 'orange'))
abline(-0.7784,0)
############################################
#linear Mixed Models for haplotypes in PM only 

PM= Risp[Risp$CYP2D6_ph1=='PM',]

LnMR= PM$lnMR2
TAD= PM$TAD_impute
Age= PM$Age
ID= PM$GeneticID
Haplotype= PM$Haplotype_All_Phenotypes_1_2_3

LMM= lme(LnMR ~ Haplotype + TAD + Age, random = ~1|ID, method = 'ML')
summary(LMM)

Haplotype1= PM$LnMR_Age_TAD_Corrected[PM$Haplotype_All_Phenotypes_1_2_3==1]
Haplotype2= PM$LnMR_Age_TAD_Corrected[PM$Haplotype_All_Phenotypes_1_2_3==2]

Haplotype3= PM$LnMR_Age_TAD_Corrected[PM$Haplotype_All_Phenotypes_1_2_3==3]
Haplotype4= PM$LnMR_Age_TAD_Corrected[PM$Haplotype_All_Phenotypes_1_2_3==4]

boxplot(Haplotype1,Haplotype2,Haplotype3,xlab= 'Diplotype', main='', 
        ylab= 'Ln(MR)', col= c('deepskyblue', 'violet', 'green3', 'red', 'deepskyblue3','purple', 'violet', 'green3',
                                            'green2', 'indianred1','yellow3', 'yellow',
                                            'aquamarine3', 'orange'))
abline(-0.7784,0)
                                            
vioplot(Haplotype1,Haplotype2,Haplotype3,xlab= 'Diplotype', main='', 
      ylab= 'Ln(MR)', col= c('deepskyblue', 'violet', 'green3', 'red', 'deepskyblue3'
                                          ,'purple', 'violet', 'green3','green2', 
                                          'indianred1','yellow3', 'yellow',
                                          'aquamarine3', 'orange'))
abline(-0.7784,0)
####################################################
#Linear Mixed Models for all phenotypes

Risp$DummyVariable_CYP2D6_Short2= rep(NA)
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*1/*1']= 1
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*1/*1x2']= 2
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*1/*41']= 3
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*1/*9-10']= 4
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*1/Nonf']= 5
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*41/*41']= 6
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*9-10/*41']= 7
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='*9-10/*9-10']= 8
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='Nonf/*41']= 9
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='Nonf/*9-10']= 10
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='Nonf/Nonf']= 11
Risp$DummyVariable_CYP2D6_Short2[Risp$CYP2D6_short2=='Unknown dup']= 12
cov1= Risp[,c(5,17,33)]

write.table(cov1, file = "File_with_age_TAD_and_CYP2D6_genotypes_as_covariates.txt"
            ,sep='\t', col.names = F, row.names = F)

LnMR_correct_TAD_Age_CYP2D6Short2=read.delim("File_with_corrected_LnMR.txt"
                                             ,sep='\t',header = F)

LnMR= Risp$lnMR2
age= Risp$Age
TAD= Risp$TAD_impute
ID= Risp$GeneticID
CYP= Risp$CYP2D6_short2
Haplotype= Risp$Haplotype_All_Phenotypes_1_2_3

LMM= lme(LnMR ~ Haplotype + age+ TAD + CYP, random = ~1|ID)
summary(LMM)

boxplot(Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==1], 
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==2],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==3],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==4],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==5],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==6],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==7],xlab= 'Diplotype', 
        ylab= 'Ln(MR)', main= '',
        col= c('purple', 'green1', 'deepskyblue1','red','yellow1', 'salmon2', 'violet'))
abline(-0.01907,0)

vioplot(Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==1], 
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==2],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==3],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==4],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==5],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==6],
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$Haplotype_All_Phenotypes_1_2_3==7],xlab= 'Diplotype', 
        ylab= 'Ln(MR)', main= '',
        col= c('purple', 'green1', 'deepskyblue1','red','yellow1', 'salmon2', 'violet'))
abline(-0.01907,0)

###############################################
#Linear Mixed Models for the TTC haplotype in all phenotypes

LnMR= Risp$lnMR2
age= Risp$Age
TAD= Risp$TAD_impute
ID= Risp$GeneticID
CYP= Risp$CYP2D6_short2
TTC= Risp$TTC_Presence
LMM= lme(LnMR ~ TTC + age+ TAD + CYP, random = ~1|ID)
summary(LMM)

vioplot(Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$TTC_Presence==F], 
        Risp$LnMR_Age_TAD_CYP2D6Short2_Corrected[Risp$TTC_Presence==T],
        xlab= 'TTC Presence', 
        ylab= 'Ln(MR)', main= '',
        col= c('purple', 'aquamarine'))
abline(-0.02685,0)





