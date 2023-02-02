PHASE= data.frame(matrix(NA, ncol= 3, nrow= nrow(data_final)*3))
colnames(PHASE)<- c('rs1135612', 'rs2286822', 'rs2302432')
PHASE[c(TRUE, FALSE, FALSE),1]=paste0('#',c(1:309))
PHASE[c(F, T, T),1]=0
PHASE[c(F, T, T),2]=0
PHASE[c(F, T, T),3]=0
PHASE[is.na(PHASE)]=''

AGSNP1= grep('A/G', data_final$POR_rs1135612)
AGSNP1= AGSNP1*3
PHASE[AGSNP1,1]=1

GGSNP1= grep('G/G',data_final$POR_rs1135612)
GGSNP1=GGSNP1*3
PHASE[GGSNP1,1]=1
GGSNP1=GGSNP1-1
PHASE[GGSNP1,1]=1

CTSNP2= grep('C/T', data_final$POR_rs2286822)
CTSNP2= CTSNP2*3
PHASE[CTSNP2,2]=1

TTSNP2= grep('T/T',data_final$POR_rs2286822)
TTSNP2=TTSNP2*3
PHASE[TTSNP2,2]=1
TTSNP2=TTSNP2-1
PHASE[TTSNP2,2]=1

TGSNP3= grep('T/G', data_final$POR_rs2302432)
TGSNP3= TGSNP3*3-1
PHASE[TGSNP3,3]=1

TTSNP3= grep('T/T', data_final$POR_rs2302432)
TTSNP3=TTSNP3*3
PHASE[TTSNP3,3]=1
TTSNP3=TTSNP3-1
PHASE[TTSNP3,3]=1
write.table(PHASE, file = 'C:/Users/mahmo/Downloads/MMIT/RPII/PHASE/POR_SNPs.txt', col.names = F
            ,row.names=F, sep= '\t')

