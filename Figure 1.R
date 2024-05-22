#Storing periodwise data
gene1_periodwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 1_Periodwise.csv")
gene2_periodwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 2_Periodwise.csv")
gene3_periodwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 3_Periodwise.csv")

#Storing treatmentwise data
gene1_treatmentwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 1_Treatmentwise.csv")
gene2_treatmentwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 2_Treatmentwise.csv")
gene3_treatmentwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 3_Treatmentwise.csv")

#Storing subjectwise data
gene1_subjectwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 1_Subjectwise.csv")
gene2_subjectwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 2_Subjectwise.csv")
gene3_subjectwise=read.csv("D:/Paper 2 Gene Leaker1/Gene 3_Subjectwise.csv")


library(dplyr)
library(tidyr)
library(ggplot2)

#For Period effect
Gene=as.factor(c("Gene 1", "Gene 2", "Gene 3"))

P1=c(sum(gene1_periodwise$Period.1)/16, sum(gene2_periodwise$Period.1)/16, sum(gene3_periodwise$Period.1)/16)
P2=c(sum(gene1_periodwise$Period.2)/16, sum(gene2_periodwise$Period.2)/16, sum(gene3_periodwise$Period.2)/16)
P3=c(sum(gene1_periodwise$Period.3)/16, sum(gene2_periodwise$Period.3)/16, sum(gene3_periodwise$Period.3)/16)

df=data.frame(Gene,P1,P2,P3)
colnames(df)=c("Gene", "1", "2", "3")

df1=df %>% gather(var, val , -Gene)
ggplot(df1)+geom_line(aes(var, val, group=Gene, linetype = Gene),lwd=1)+xlab("Period") + ylab("Response")


#For Subject effect
Gene=as.factor(c("Gene 1", "Gene 2", "Gene 3"))

S1=c(sum(gene1_subjectwise$Subject.1)/3,sum(gene2_subjectwise$Subject.1)/3, sum(gene3_subjectwise$Subject.1)/3)
S2=c(sum(gene1_subjectwise$Subject.2)/3,sum(gene2_subjectwise$Subject.2)/3, sum(gene3_subjectwise$Subject.2)/3)
S3=c(sum(gene1_subjectwise$Subject.3)/3,sum(gene2_subjectwise$Subject.3)/3, sum(gene3_subjectwise$Subject.3)/3)
S4=c(sum(gene1_subjectwise$Subject.4)/3,sum(gene2_subjectwise$Subject.4)/3, sum(gene3_subjectwise$Subject.4)/3)
S5=c(sum(gene1_subjectwise$Subject.5)/3,sum(gene2_subjectwise$Subject.5)/3, sum(gene3_subjectwise$Subject.5)/3)
S6=c(sum(gene1_subjectwise$Subject.6)/3,sum(gene2_subjectwise$Subject.6)/3, sum(gene3_subjectwise$Subject.6)/3)
S7=c(sum(gene1_subjectwise$Subject.7)/3,sum(gene2_subjectwise$Subject.7)/3, sum(gene3_subjectwise$Subject.7)/3)
S8=c(sum(gene1_subjectwise$Subject.8)/3,sum(gene2_subjectwise$Subject.8)/3, sum(gene3_subjectwise$Subject.8)/3)
S9=c(sum(gene1_subjectwise$Subject.9)/3,sum(gene2_subjectwise$Subject.9)/3, sum(gene3_subjectwise$Subject.9)/3)
S10=c(sum(gene1_subjectwise$Subject.10)/3,sum(gene2_subjectwise$Subject.10)/3, sum(gene3_subjectwise$Subject.10)/3)
S11=c(sum(gene1_subjectwise$Subject.11)/3,sum(gene2_subjectwise$Subject.11)/3, sum(gene3_subjectwise$Subject.11)/3)
S12=c(sum(gene1_subjectwise$Subject.12)/3,sum(gene2_subjectwise$Subject.12)/3, sum(gene3_subjectwise$Subject.12)/3)
S13=c(sum(gene1_subjectwise$Subject.13)/3,sum(gene2_subjectwise$Subject.13)/3, sum(gene3_subjectwise$Subject.13)/3)
S14=c(sum(gene1_subjectwise$Subject.14)/3,sum(gene2_subjectwise$Subject.14)/3, sum(gene3_subjectwise$Subject.14)/3)
S15=c(sum(gene1_subjectwise$Subject.15)/3,sum(gene2_subjectwise$Subject.15)/3, sum(gene3_subjectwise$Subject.15)/3)
S16=c(sum(gene1_subjectwise$Subject.16)/3,sum(gene2_subjectwise$Subject.16)/3, sum(gene3_subjectwise$Subject.16)/3)

df=data.frame(Gene,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16)
colnames(df) = c("Gene", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")

df1=df %>% gather(var, val , -Gene)
ggplot(df1)+geom_line(aes(var, val, group=Gene, linetype = Gene),lwd=1)+xlab("Subject") + ylab("Response")+theme(axis.text.x = element_text(size = 6, angle = 0, hjust = .5, vjust = .5, face = "plain"))


#For Treatment effect
Gene=as.factor(c("Gene 1", "Gene 2", "Gene 3"))

A=c(sum(gene1_treatmentwise$X10mg[1:12])/12, sum(gene2_treatmentwise$X10mg[1:12])/12, sum(gene3_treatmentwise$X10mg[1:12])/12)
B=c(sum(gene1_treatmentwise$placebo)/18, sum(gene2_treatmentwise$placebo)/18, sum(gene3_treatmentwise$placebo)/18)
C=c(sum(gene1_treatmentwise$X25mg)/18, sum(gene2_treatmentwise$X25mg)/18, sum(gene3_treatmentwise$X25mg)/18)

df=data.frame(Gene,A,B,C)
colnames(df) = c("Gene", "A", "B", "C")
df1=df %>% gather(var, val , -Gene)
ggplot(df1)+geom_line(aes(var, val, group=Gene, linetype = Gene),lwd=1)+xlab("Treatment") + ylab("Response")


#For Overall effect
Gene=as.factor(c("Gene 1", "Gene 2", "Gene 3"))
gene_only=read.csv("D:/Paper 2 Gene Leaker1/Gene only.csv")

gene_overall=c(sum(gene_only$Gene.1)/48,sum(gene_only$Gene.2)/48,sum(gene_only$Gene.3)/48)

df = data.frame(Gene, gene_overall)

ggplot(df, aes(Gene, gene_overall, ymax = gene_overall, ymin = 0)) +
  geom_pointrange(color = "black")+xlab("Gene") + ylab("Response")

