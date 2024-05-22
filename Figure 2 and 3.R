library(readxl)

G1=read_xlsx("D:/Paper 2 Gene Leaker1/Gene 1.xlsx")
G2=read_xlsx("D:/Paper 2 Gene Leaker1/Gene 2.xlsx")
G3=read_xlsx("D:/Paper 2 Gene Leaker1/Gene 3.xlsx")

#Checking Randomness
library(snpar)

#Gene 1
runs.test(array(unlist(G1[1,])), exact=TRUE)
runs.test(array(unlist(G1[2,])), exact=TRUE)
runs.test(array(unlist(G1[3,])), exact=TRUE)

#Gene 2
runs.test(array(unlist(G2[1,])), exact=TRUE)
runs.test(array(unlist(G2[2,])), exact=TRUE)
runs.test(array(unlist(G2[3,])), exact=TRUE)

#Gene 3
runs.test(array(unlist(G3[1,])), exact=TRUE)
runs.test(array(unlist(G3[2,])), exact=TRUE)
runs.test(array(unlist(G3[3,])), exact=TRUE)

#The data on genes is randomly distributed at 5% level of significance

#Checking normality for each period

#Gene 1
shapiro.test(array(unlist(G1[1,])))
shapiro.test(array(unlist(G1[2,])))
shapiro.test(array(unlist(G1[3,])))

#Gene 2
shapiro.test(array(unlist(G2[1,])))
shapiro.test(array(unlist(G2[2,])))
shapiro.test(array(unlist(G2[3,])))

#Gene 3
shapiro.test(array(unlist(G3[1,])))
shapiro.test(array(unlist(G3[2,])))
shapiro.test(array(unlist(G3[3,])))

#The data on genes is normally distributed at 5% level of significance

g1_p1=array(unlist(G1[1,]))
g1_p2=array(unlist(G1[2,]))
g1_p3=array(unlist(G1[3,]))


g2_p1=array(unlist(G2[1,]))
g2_p2=array(unlist(G2[2,]))          
g2_p3=array(unlist(G2[3,]))


g3_p1=array(unlist(G3[1,]))
g3_p2=array(unlist(G3[2,]))          
g3_p3=array(unlist(G3[3,]))   

#Preparing dataframe for checking bivariate normality....Same gene different periods

G1_P1P2=data.frame(g1_p1, g1_p2)
G1_P1P3=data.frame(g1_p1, g1_p3)
G1_P2P3=data.frame(g1_p2, g1_p3)

G2_P1P2=data.frame(g2_p1, g2_p2)
G2_P1P3=data.frame(g2_p1, g2_p3)
G2_P2P3=data.frame(g2_p2, g2_p3)

G3_P1P2=data.frame(g3_p1, g3_p2)
G3_P1P3=data.frame(g3_p1, g3_p3)
G3_P2P3=data.frame(g3_p2, g3_p3)

#Checking Bivariate Normality....Same gene different periods

library(QuantPsyc)
mult.norm(G1_P1P2)$mult.test
mult.norm(G2_P1P2)$mult.test
mult.norm(G3_P1P2)$mult.test

mult.norm(G1_P1P3)$mult.test
mult.norm(G2_P1P3)$mult.test
mult.norm(G3_P1P3)$mult.test

mult.norm(G1_P2P3)$mult.test
mult.norm(G2_P2P3)$mult.test
mult.norm(G3_P2P3)$mult.test

#The data corresponding to same gene but different periods is bivariate normally distributed at 5% level of significance

#Preparing dataframe for checking bivariate normality....Different gene same and different periods

#G1-G2

G1G2_P1=data.frame(g1_p1,g2_p1)
G1G2_P2=data.frame(g1_p2,g2_p2)
G1G2_P3=data.frame(g1_p3,g2_p3)
G1G2_P1P2=data.frame(g1_p1,g2_p2)
G1G2_P1P3=data.frame(g1_p1,g2_p3)
G1G2_P2P1=data.frame(g1_p2,g2_p1)
G1G2_P2P3=data.frame(g1_p2,g2_p3)
G1G2_P3P1=data.frame(g1_p3,g2_p1)
G1G2_P3P2=data.frame(g1_p3,g2_p2)

#G1-G3

G1G3_P1=data.frame(g1_p1,g3_p1)
G1G3_P2=data.frame(g1_p2,g3_p2)
G1G3_P3=data.frame(g1_p3,g3_p3)
G1G3_P1P2=data.frame(g1_p1,g3_p2)
G1G3_P1P3=data.frame(g1_p1,g3_p3)
G1G3_P2P1=data.frame(g1_p2,g3_p1)
G1G3_P2P3=data.frame(g1_p2,g3_p3)
G1G3_P3P1=data.frame(g1_p3,g3_p1)
G1G3_P3P2=data.frame(g1_p3,g3_p2)

#G2-G3

G2G3_P1=data.frame(g2_p1,g3_p1)
G2G3_P2=data.frame(g2_p2,g3_p2)
G2G3_P3=data.frame(g2_p3,g3_p3)
G2G3_P1P2=data.frame(g2_p1,g3_p2)
G2G3_P1P3=data.frame(g2_p1,g3_p3)
G2G3_P2P1=data.frame(g2_p2,g3_p1)
G2G3_P2P3=data.frame(g2_p2,g3_p3)
G2G3_P3P1=data.frame(g2_p3,g3_p1)
G2G3_P3P2=data.frame(g2_p3,g3_p2)

#Checking bivariate normality....Different gene same and different periods

#G1-G2

mult.norm(G1G2_P1)$mult.test
mult.norm(G1G2_P2)$mult.test
mult.norm(G1G2_P3)$mult.test
mult.norm(G1G2_P1P2)$mult.test
mult.norm(G1G2_P1P3)$mult.test
mult.norm(G1G2_P2P1)$mult.test
mult.norm(G1G2_P2P3)$mult.test
mult.norm(G1G2_P3P1)$mult.test
mult.norm(G1G2_P3P2)$mult.test

#G1-G3

mult.norm(G1G3_P1)$mult.test
mult.norm(G1G3_P2)$mult.test
mult.norm(G1G3_P3)$mult.test
mult.norm(G1G3_P1P2)$mult.test
mult.norm(G1G3_P1P3)$mult.test
mult.norm(G1G3_P2P1)$mult.test
mult.norm(G1G3_P2P3)$mult.test
mult.norm(G1G3_P3P1)$mult.test
mult.norm(G1G3_P3P2)$mult.test

#G2-G3

mult.norm(G2G3_P1)$mult.test
mult.norm(G2G3_P2)$mult.test
mult.norm(G2G3_P3)$mult.test
mult.norm(G2G3_P1P2)$mult.test
mult.norm(G2G3_P1P3)$mult.test
mult.norm(G2G3_P2P1)$mult.test
mult.norm(G2G3_P2P3)$mult.test
mult.norm(G2G3_P3P1)$mult.test
mult.norm(G2G3_P3P2)$mult.test

#The data corresponding to different gene but same and different periods is bivariate normally distributed at 5% level of significance

#Correlation test and plot....Same gene different period

library(stats)
library(corrplot)

#G1

g1=data.frame(cbind(g1_p1,g1_p2,g1_p3))
colnames(g1)=c("Period 1","Period 2","Period 3")
cor_mat_g1 = cor(g1)

testRes1=cor.mtest(g1)$p
corrplot(cor_mat_g1, p.mat = testRes1, insig = 'p-value', sig.level = -1)


#G2

g2=data.frame(cbind(g2_p1,g2_p2,g2_p3))
colnames(g2)=c("Period 1","Period 2","Period 3")
cor_mat_g2 = cor(g2)

testRes2=cor.mtest(g2)$p
corrplot(cor_mat_g2, p.mat = testRes2, insig = 'p-value', sig.level = -1)


#G3

g3=data.frame(cbind(g3_p1,g3_p2,g3_p3))
colnames(g3)=c("Period 1","Period 2","Period 3")
cor_mat_g3 = cor(g3)

testRes3=cor.mtest(g3)$p
corrplot(cor_mat_g3, p.mat = testRes3, insig = 'p-value', sig.level = -1)

#Correlation test and plot....Different gene same and different period

#G1-G2

cor.test(g1_p1, g2_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g1_p1, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g2_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g2_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor_mat_g1g2 = cor(g1,g2)
testRes12=matrix(c(0.01711,0.2017,0.03479,0.9029,0.07621,0.8801,0.7112,0.007014,0.5464),nrow=3,ncol=3,byrow=TRUE)
rownames(testRes12)=c("Period 1","Period 2","Period 3")
colnames(testRes12)=c("Period 1","Period 2","Period 3")

corrplot(cor_mat_g1g2, p.mat = testRes12, insig = 'p-value', sig.level = -1)

#G1-G3

cor.test(g1_p1, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g1_p1, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor_mat_g1g3 = cor(g1,g3)
testRes13=matrix(c(0.2191,0.9664,0.3615,0.4846,0.4818,0.2814,0.2613,0.2842,0.2344),nrow=3,ncol=3,byrow=TRUE)
rownames(testRes13)=c("Period 1","Period 2","Period 3")
colnames(testRes13)=c("Period 1","Period 2","Period 3")

corrplot(cor_mat_g1g3, p.mat = testRes13, insig = 'p-value', sig.level = -1)

#G2-G3

cor.test(g2_p1, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g2_p1, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p1, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)


cor_mat_g2g3 = cor(g2,g3)
testRes23=matrix(c(0.3027,0.9257,0.454,0.2415,0.3283,0.6428,0.006652,0.4312,0.2856),nrow=3,ncol=3,byrow=TRUE)
rownames(testRes23)=c("Period 1","Period 2","Period 3")
colnames(testRes23)=c("Period 1","Period 2","Period 3")

corrplot(cor_mat_g2g3, p.mat = testRes23, insig = 'p-value', sig.level = -1)

