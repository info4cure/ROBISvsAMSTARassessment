######################################################################################################################################
# --- 
# title: "Most systematic reviews of high methodological quality on psoriasis interventions are classified as high risk of bias using ROBIS tool" 
# author: "Juan Ruano" 
# date: "16 May 2017" 
# institutions: Department of Dermatology, IMIBIC/Reina Sofia University Hospital/University of Cordoba, Cordoba, Spain
# analysis: Statistical analyses: PCA, Likert scales, correlation matrix and network, and Radviz plots
# --- 
# 
# R version 3.3.1 (2016-06-21)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.9.5 (Mavericks)
# 
######################################################################################################################################

########  read .csv files ----------------

file2 <- read.csv2("ROBIS_vs_AMSTAR_pca.csv", 
                 sep = ";", 
                 dec= ".", 
                 stringsAsFactors = TRUE,
                 header = TRUE)
DB2<-as.data.frame(file2)
names(DB2)
head(DB2)

df_amstar <- DB2[, 5:15]
res.pca <- PCA(df, graph = FALSE, scale=T)
round(res.pca$ind$contrib, 4)

DB2$AMSTAR_levels_2 <- factor(DB2$AMSTAR_levels_2, 
                              levels = c("high_quality", "moderate_quality","low_quality" ))


############################################################################
######## 1: PCA ----------------
######## packages ----------------

library("ggfortify")
library("cluster")
library("reshape2")
library("FactoMineR")
library("factoextra")

######## PCA amstar ----------------

DB2_PCA_amstar <- na.omit(DB2[c(30:40,41,42, 3)])
names(DB2_PCA_amstar)
amstar.pca <- prcomp(DB2_PCA_amstar[c(1:11)])
fviz_screeplot(amstar.pca, ncp=10)

plot_pca_amstar <- autoplot(prcomp(DB2_PCA_amstar[c(1:11)]), 
                            data = DB2_PCA_amstar, 
                            colour = 'AMSTAR_levels_2', 
                            shape="RoB_ROBIS", 
                            loadings = TRUE, 
                            loadings.label = TRUE, 
                            frame = TRUE)
plot_pca_amstar +
  geom_point(aes(shape = factor(RoB_ROBIS), 
                 colour = AMSTAR_levels_2)) +
  scale_shape_manual(values=c(24,25))+
  geom_jitter(aes(shape = factor(RoB_ROBIS), 
                  colour = AMSTAR_levels_2), 
                  position=position_jitter(height=0.009,width = 0.009))

#### plot AMSTAR variable contribution axes

fviz_pca_var(amstar.pca, col.var="contrib") +
  scale_color_gradient2(low="white", mid="blue",high="red", midpoint=10) + 
  theme_minimal()

#### plot top10 AMSTAR variables contribution to PC1/PC2

fviz_pca_contrib(amstar.pca, 
                 choice = "var", 
                 axes = 1:2,  
                 fill="coral2", 
                 color = "grey") +
  theme(axis.text.x = element_text(size=10, angle=90)) +
  scale_y_continuous(limits = c(0.0, 25.0))

#### plot top10 reviews contribution to AMSTAR PC1/PC2

fviz_pca_contrib(amstar.pca, 
                 choice = "ind", 
                 axes = 1:2,  
                 top=50,
                 fill="seagreen4", 
                 color = "grey") +
  theme(axis.text.x = element_text( size=10, angle=90)) +
  scale_y_continuous(limits = c(0.0, 3.0))

#### plot review contributions to AMSTAR PC1/PC2

fviz_pca_ind(amstar.pca, 
             col.ind="contrib", 
             jitter = list(what = "label", width = NULL, height = NULL)) +
  scale_color_gradient2(low="white", mid="blue", high="red", midpoint=1) + 
  theme_minimal() +
  geom_point(colour="white") 


######## PCA robis ----------------

DB2_PCA_robis <- na.omit(DB2[c(4:24,41,42, 3)])

robis.pca <- prcomp(DB2_PCA_robis[c(1:21)])
fviz_screeplot(robis.pca, ncp=10)

plot_pca_robis <- autoplot(prcomp(DB2_PCA_robis[c(1:21)]), 
                           data = DB2_PCA_robis, 
                           colour = 'RoB_ROBIS', 
                           shape="AMSTAR_levels_2", 
                           loadings = TRUE, 
                           loadings.label = TRUE, 
                           frame = TRUE)
plot_pca_robis +
  geom_point(aes(shape = factor(AMSTAR_levels_2), 
                 colour = RoB_ROBIS)) +
  scale_shape_manual(values=c(24,4,25))+
  geom_jitter(aes(shape = factor(AMSTAR_levels_2), 
                  colour = RoB_ROBIS), 
                  position=position_jitter(height=0.009,width = 0.009))

#### plot ROBIS variable contribution axes

fviz_pca_var(robis.pca, col.var="contrib") +
  scale_color_gradient2(low="white", mid="blue", high="red", midpoint=10) + 
  theme_minimal()

#### plot top10 ROBIS variables contribution to PC1/PC2

fviz_pca_contrib(robis.pca, choice = "var", axes = 1:2,  fill="coral2", color = "grey") +
  theme(axis.text.x = element_text( size=10, angle=90)) +
  scale_y_continuous(limits = c(0.0, 25.0))

#### plot top10 individuals contribution to ROBIS PC1/PC2

fviz_pca_contrib(robis.pca, 
                 choice = "ind", 
                 axes = 1:2,  
                 top=50,fill="seagreen4", 
                 color = "grey") +
  theme(axis.text.x = element_text( size=10, angle=90)) +
  scale_y_continuous(limits = c(0.0, 3.0))

#### plot review contributions to ROBIS PC1/PC2
fviz_pca_ind(robis.pca, 
             col.ind="contrib", 
             jitter = list(what = "label", width = NULL, height = NULL)) +
scale_color_gradient2(low="white", mid="blue", high="red", midpoint=1) + 
             theme_minimal() +
             geom_point(colour="white") 


################################
#### Plot features that contribute to the classification only ROBIS
robis.pca

df_out_r <- as.data.frame(robis.pca$rotation)
df_out_r$feature <- row.names(df_out_r)

df_out_r

p<-ggplot(df_out_r,aes(x=PC1,y=PC2,label=feature))
p<-p+geom_point(colour="white")+ geom_text(size=3, check_overlap = TRUE)
p

##### mosaic plot ROBIS vs AMSTAR
mosaicplot(DB2$AMSTAR_levels_2~DB2$RoB_ROBIS,  data=DB2, color = c(23,18), las = 1, xlab="Methodological quality (AMSTAR)", ylab="Risk of Bias (ROBIS)", main="")



############################################################################
########## 2: Likert scales ALL ----------------

######## packages --------------------

require(likert)
require(grid)
require(lattice)
require(latticeExtra)
require(HH)
library(plyr)


######## tyding data ----------------

sgbar.likert <- DB2[,c(30:40, 42,3)]
sgbar.likert$Q1_AMSTAR <- as.factor(sgbar.likert$Q1_AMSTAR)
sgbar.likert$Q1_AMSTAR <- revalue(sgbar.likert$Q1_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q2_AMSTAR <- as.factor(sgbar.likert$Q2_AMSTAR)
sgbar.likert$Q2_AMSTAR <- revalue(sgbar.likert$Q2_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q3_AMSTAR <- as.factor(sgbar.likert$Q3_AMSTAR)
sgbar.likert$Q3_AMSTAR <- revalue(sgbar.likert$Q3_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q4_AMSTAR <- as.factor(sgbar.likert$Q4_AMSTAR)
sgbar.likert$Q4_AMSTAR <- revalue(sgbar.likert$Q4_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q5_AMSTAR <- as.factor(sgbar.likert$Q5_AMSTAR)
sgbar.likert$Q5_AMSTAR <- revalue(sgbar.likert$Q5_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q6_AMSTAR <- as.factor(sgbar.likert$Q6_AMSTAR)
sgbar.likert$Q6_AMSTAR <- revalue(sgbar.likert$Q6_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q7_AMSTAR <- as.factor(sgbar.likert$Q7_AMSTAR)
sgbar.likert$Q7_AMSTAR <- revalue(sgbar.likert$Q7_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q8_AMSTAR <- as.factor(sgbar.likert$Q8_AMSTAR)
sgbar.likert$Q8_AMSTAR <- revalue(sgbar.likert$Q8_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q9_AMSTAR <- as.factor(sgbar.likert$Q9_AMSTAR)
sgbar.likert$Q9_AMSTAR <- revalue(sgbar.likert$Q9_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q10_AMSTAR <- as.factor(sgbar.likert$Q10_AMSTAR)
sgbar.likert$Q10_AMSTAR <- revalue(sgbar.likert$Q10_AMSTAR, c("0"="no", "1"="yes"))

sgbar.likert$Q11_AMSTAR <- as.factor(sgbar.likert$Q11_AMSTAR)
sgbar.likert$Q11_AMSTAR <- revalue(sgbar.likert$Q11_AMSTAR, c("0"="no", "1"="yes"))




sgbar.likert_ROBIS <- DB2[,c(4:24, 42,3)]
sgbar.likert_ROBIS$Q11 <- as.factor(sgbar.likert_ROBIS$Q11)
sgbar.likert_ROBIS$Q11 <- revalue(sgbar.likert_ROBIS$Q11, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q12 <- as.factor(sgbar.likert_ROBIS$Q12)
sgbar.likert_ROBIS$Q12 <- revalue(sgbar.likert_ROBIS$Q12, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q13 <- as.factor(sgbar.likert_ROBIS$Q13)
sgbar.likert_ROBIS$Q13 <- revalue(sgbar.likert_ROBIS$Q13, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q14 <- as.factor(sgbar.likert_ROBIS$Q14)
sgbar.likert_ROBIS$Q14 <- revalue(sgbar.likert_ROBIS$Q14, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q15 <- as.factor(sgbar.likert_ROBIS$Q15)
sgbar.likert_ROBIS$Q15 <- revalue(sgbar.likert_ROBIS$Q15, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q21 <- as.factor(sgbar.likert_ROBIS$Q21)
sgbar.likert_ROBIS$Q21 <- revalue(sgbar.likert_ROBIS$Q21, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q22 <- as.factor(sgbar.likert_ROBIS$Q22)
sgbar.likert_ROBIS$Q22 <- revalue(sgbar.likert_ROBIS$Q22, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q23 <- as.factor(sgbar.likert_ROBIS$Q23)
sgbar.likert_ROBIS$Q23 <- revalue(sgbar.likert_ROBIS$Q23, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q24 <- as.factor(sgbar.likert_ROBIS$Q24)
sgbar.likert_ROBIS$Q24 <- revalue(sgbar.likert_ROBIS$Q24, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q25 <- as.factor(sgbar.likert_ROBIS$Q25)
sgbar.likert_ROBIS$Q25 <- revalue(sgbar.likert_ROBIS$Q25, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q31 <- as.factor(sgbar.likert_ROBIS$Q31)
sgbar.likert_ROBIS$Q31 <- revalue(sgbar.likert_ROBIS$Q31, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q32 <- as.factor(sgbar.likert_ROBIS$Q32)
sgbar.likert_ROBIS$Q32 <- revalue(sgbar.likert_ROBIS$Q32, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q33 <- as.factor(sgbar.likert_ROBIS$Q33)
sgbar.likert_ROBIS$Q33 <- revalue(sgbar.likert_ROBIS$Q33, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q34 <- as.factor(sgbar.likert_ROBIS$Q34)
sgbar.likert_ROBIS$Q34 <- revalue(sgbar.likert_ROBIS$Q34, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q35 <- as.factor(sgbar.likert_ROBIS$Q35)
sgbar.likert_ROBIS$Q35 <- revalue(sgbar.likert_ROBIS$Q35, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q41 <- as.factor(sgbar.likert_ROBIS$Q41)
sgbar.likert_ROBIS$Q41 <- revalue(sgbar.likert_ROBIS$Q41, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q42 <- as.factor(sgbar.likert_ROBIS$Q42)
sgbar.likert_ROBIS$Q42 <- revalue(sgbar.likert_ROBIS$Q42, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q43 <- as.factor(sgbar.likert_ROBIS$Q43)
sgbar.likert_ROBIS$Q43 <- revalue(sgbar.likert_ROBIS$Q43, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q44 <- as.factor(sgbar.likert_ROBIS$Q44)
sgbar.likert_ROBIS$Q44 <- revalue(sgbar.likert_ROBIS$Q44, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q45 <- as.factor(sgbar.likert_ROBIS$Q45)
sgbar.likert_ROBIS$Q45 <- revalue(sgbar.likert_ROBIS$Q45, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))

sgbar.likert_ROBIS$Q46 <- as.factor(sgbar.likert_ROBIS$Q46)
sgbar.likert_ROBIS$Q46 <- revalue(sgbar.likert_ROBIS$Q46, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))


######## plotting Likerts scales ----------------
##### likert scale AMSTAR

plot(likert(sgbar.likert[1:11], 
            grouping=sgbar.likert$RoB_ROBIS), 
            main="", 
            ylab="AMSTAR items")

##### likert scale ROBIS

plot(likert(sgbar.likert_ROBIS[1:21], 
            grouping=sgbar.likert_ROBIS$AMSTAR_levels_2), 
            main="", 
            ylab="ROBIS items")

##### ONLY high quality subgroup Likert scales

DB2_high<-subset(DB2, AMSTAR_levels_2=="high_quality")
DB2_high_robis_H<-subset(DB2_high, RoB_ROBIS=="HIGH")
DB2_high_robis_L<-subset(DB2_high, RoB_ROBIS=="LOW")

### ROBIS HIGH

sgbar.likert_ROBIS_high_H <- DB2_high_robis_H[,c(4:24)]
desired.order <- c("Low", "Probably low","Probably high", "High")
sgbar.likert_ROBIS_high_H$Q11 <- as.factor(sgbar.likert_ROBIS_high_H$Q11)
sgbar.likert_ROBIS_high_H$Q11 <- revalue(sgbar.likert_ROBIS_high_H$Q11, c("0"="Low", "1"="Probably low","2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q11 <- factor(sgbar.likert_ROBIS_high_H$Q11, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q12 <- as.factor(sgbar.likert_ROBIS_high_H$Q12)
sgbar.likert_ROBIS_high_H$Q12 <- revalue(sgbar.likert_ROBIS_high_H$Q12, c("0"="Low","1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q12 <- factor(sgbar.likert_ROBIS_high_H$Q12, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q13 <- as.factor(sgbar.likert_ROBIS_high_H$Q13)
sgbar.likert_ROBIS_high_H$Q13 <- revalue(sgbar.likert_ROBIS_high_H$Q13, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q13 <- factor(sgbar.likert_ROBIS_high_H$Q13, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q14 <- as.factor(sgbar.likert_ROBIS_high_H$Q14)
sgbar.likert_ROBIS_high_H$Q14 <- revalue(sgbar.likert_ROBIS_high_H$Q14, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q14 <- factor(sgbar.likert_ROBIS_high_H$Q14, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q15 <- as.factor(sgbar.likert_ROBIS_high_H$Q15)
sgbar.likert_ROBIS_high_H$Q15 <- revalue(sgbar.likert_ROBIS_high_H$Q15, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q15 <- factor(sgbar.likert_ROBIS_high_H$Q15, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q21 <- as.factor(sgbar.likert_ROBIS_high_H$Q21)
sgbar.likert_ROBIS_high_H$Q21 <- revalue(sgbar.likert_ROBIS_high_H$Q21, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q21 <- factor(sgbar.likert_ROBIS_high_H$Q21, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q22 <- as.factor(sgbar.likert_ROBIS_high_H$Q22)
sgbar.likert_ROBIS_high_H$Q22 <- revalue(sgbar.likert_ROBIS_high_H$Q22, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q22 <- factor(sgbar.likert_ROBIS_high_H$Q22, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q23 <- as.factor(sgbar.likert_ROBIS_high_H$Q23)
sgbar.likert_ROBIS_high_H$Q23 <- revalue(sgbar.likert_ROBIS_high_H$Q23, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q23 <- factor(sgbar.likert_ROBIS_high_H$Q23, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q24 <- as.factor(sgbar.likert_ROBIS_high_H$Q24)
sgbar.likert_ROBIS_high_H$Q24 <- revalue(sgbar.likert_ROBIS_high_H$Q24, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q24 <- factor(sgbar.likert_ROBIS_high_H$Q24, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q25 <- as.factor(sgbar.likert_ROBIS_high_H$Q25)
sgbar.likert_ROBIS_high_H$Q25 <- revalue(sgbar.likert_ROBIS_high_H$Q25, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q25 <- factor(sgbar.likert_ROBIS_high_H$Q25, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q31 <- as.factor(sgbar.likert_ROBIS_high_H$Q31)
sgbar.likert_ROBIS_high_H$Q31 <- revalue(sgbar.likert_ROBIS_high_H$Q31, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q31 <- factor(sgbar.likert_ROBIS_high_H$Q31, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q32 <- as.factor(sgbar.likert_ROBIS_high_H$Q32)
sgbar.likert_ROBIS_high_H$Q32 <- revalue(sgbar.likert_ROBIS_high_H$Q32, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q32 <- factor(sgbar.likert_ROBIS_high_H$Q32, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q33 <- as.factor(sgbar.likert_ROBIS_high_H$Q33)
sgbar.likert_ROBIS_high_H$Q33 <- revalue(sgbar.likert_ROBIS_high_H$Q33, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q33 <- factor(sgbar.likert_ROBIS_high_H$Q33, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q34 <- as.factor(sgbar.likert_ROBIS_high_H$Q34)
sgbar.likert_ROBIS_high_H$Q34 <- revalue(sgbar.likert_ROBIS_high_H$Q34, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q34 <- factor(sgbar.likert_ROBIS_high_H$Q34, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q35 <- as.factor(sgbar.likert_ROBIS_high_H$Q35)
sgbar.likert_ROBIS_high_H$Q35 <- revalue(sgbar.likert_ROBIS_high_H$Q35, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q35 <- factor(sgbar.likert_ROBIS_high_H$Q35, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q41 <- as.factor(sgbar.likert_ROBIS_high_H$Q41)
sgbar.likert_ROBIS_high_H$Q41 <- revalue(sgbar.likert_ROBIS_high_H$Q41, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q41 <- factor(sgbar.likert_ROBIS_high_H$Q41, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q42 <- as.factor(sgbar.likert_ROBIS_high_H$Q42)
sgbar.likert_ROBIS_high_H$Q42 <- revalue(sgbar.likert_ROBIS_high_H$Q42, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q42 <- factor(sgbar.likert_ROBIS_high_H$Q42, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q43 <- as.factor(sgbar.likert_ROBIS_high_H$Q43)
sgbar.likert_ROBIS_high_H$Q43 <- revalue(sgbar.likert_ROBIS_high_H$Q43, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q43 <- factor(sgbar.likert_ROBIS_high_H$Q43, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q44 <- as.factor(sgbar.likert_ROBIS_high_H$Q44)
sgbar.likert_ROBIS_high_H$Q44 <- revalue(sgbar.likert_ROBIS_high_H$Q44, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q44 <- factor(sgbar.likert_ROBIS_high_H$Q44, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q45 <- as.factor(sgbar.likert_ROBIS_high_H$Q45)
sgbar.likert_ROBIS_high_H$Q45 <- revalue(sgbar.likert_ROBIS_high_H$Q45, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q45 <- factor(sgbar.likert_ROBIS_high_H$Q45, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_H$Q46 <- as.factor(sgbar.likert_ROBIS_high_H$Q46)
sgbar.likert_ROBIS_high_H$Q46 <- revalue(sgbar.likert_ROBIS_high_H$Q46, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_H$Q46 <- factor(sgbar.likert_ROBIS_high_H$Q46, levels=desired.order, ordered=TRUE)

attach(sgbar.likert_ROBIS_high_H)

######## likert scale ROBIS subgrupo high ----------------

plot(likert(sgbar.likert_ROBIS_high_H[1:21]),  
     main="High RoB (n=19)", 
     ylab="ROBIS items")


######## ROBIS LOW ----------------

sgbar.likert_ROBIS_high_L <- DB2_high_robis_L[,c(4:24)]
desired.order <- c("Low", "Probably low","Probably high", "High")
sgbar.likert_ROBIS_high_L$Q11 <- as.factor(sgbar.likert_ROBIS_high_L$Q11)
sgbar.likert_ROBIS_high_L$Q11 <- revalue(sgbar.likert_ROBIS_high_L$Q11, c("0"="Low", "1"="Probably low","2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q11 <- factor(sgbar.likert_ROBIS_high_L$Q11, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q12 <- as.factor(sgbar.likert_ROBIS_high_L$Q12)
sgbar.likert_ROBIS_high_L$Q12 <- revalue(sgbar.likert_ROBIS_high_L$Q12, c("0"="Low","1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q12 <- factor(sgbar.likert_ROBIS_high_L$Q12, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q13 <- as.factor(sgbar.likert_ROBIS_high_L$Q13)
sgbar.likert_ROBIS_high_L$Q13 <- revalue(sgbar.likert_ROBIS_high_L$Q13, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q13 <- factor(sgbar.likert_ROBIS_high_L$Q13, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q14 <- as.factor(sgbar.likert_ROBIS_high_L$Q14)
sgbar.likert_ROBIS_high_L$Q14 <- revalue(sgbar.likert_ROBIS_high_L$Q14, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q14 <- factor(sgbar.likert_ROBIS_high_L$Q14, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q15 <- as.factor(sgbar.likert_ROBIS_high_L$Q15)
sgbar.likert_ROBIS_high_L$Q15 <- revalue(sgbar.likert_ROBIS_high_L$Q15, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q15 <- factor(sgbar.likert_ROBIS_high_L$Q15, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q21 <- as.factor(sgbar.likert_ROBIS_high_L$Q21)
sgbar.likert_ROBIS_high_L$Q21 <- revalue(sgbar.likert_ROBIS_high_L$Q21, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q21 <- factor(sgbar.likert_ROBIS_high_L$Q21, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q22 <- as.factor(sgbar.likert_ROBIS_high_L$Q22)
sgbar.likert_ROBIS_high_L$Q22 <- revalue(sgbar.likert_ROBIS_high_L$Q22, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q22 <- factor(sgbar.likert_ROBIS_high_L$Q22, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q23 <- as.factor(sgbar.likert_ROBIS_high_L$Q23)
sgbar.likert_ROBIS_high_L$Q23 <- revalue(sgbar.likert_ROBIS_high_L$Q23, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q23 <- factor(sgbar.likert_ROBIS_high_L$Q23, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q24 <- as.factor(sgbar.likert_ROBIS_high_L$Q24)
sgbar.likert_ROBIS_high_L$Q24 <- revalue(sgbar.likert_ROBIS_high_L$Q24, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q24 <- factor(sgbar.likert_ROBIS_high_L$Q24, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q25 <- as.factor(sgbar.likert_ROBIS_high_L$Q25)
sgbar.likert_ROBIS_high_L$Q25 <- revalue(sgbar.likert_ROBIS_high_L$Q25, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q25 <- factor(sgbar.likert_ROBIS_high_L$Q25, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q31 <- as.factor(sgbar.likert_ROBIS_high_L$Q31)
sgbar.likert_ROBIS_high_L$Q31 <- revalue(sgbar.likert_ROBIS_high_L$Q31, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q31 <- factor(sgbar.likert_ROBIS_high_L$Q31, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q32 <- as.factor(sgbar.likert_ROBIS_high_L$Q32)
sgbar.likert_ROBIS_high_L$Q32 <- revalue(sgbar.likert_ROBIS_high_L$Q32, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q32 <- factor(sgbar.likert_ROBIS_high_L$Q32, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q33 <- as.factor(sgbar.likert_ROBIS_high_L$Q33)
sgbar.likert_ROBIS_high_L$Q33 <- revalue(sgbar.likert_ROBIS_high_L$Q33, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q33 <- factor(sgbar.likert_ROBIS_high_L$Q33, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q34 <- as.factor(sgbar.likert_ROBIS_high_L$Q34)
sgbar.likert_ROBIS_high_L$Q34 <- revalue(sgbar.likert_ROBIS_high_L$Q34, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q34 <- factor(sgbar.likert_ROBIS_high_L$Q34, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q35 <- as.factor(sgbar.likert_ROBIS_high_L$Q35)
sgbar.likert_ROBIS_high_L$Q35 <- revalue(sgbar.likert_ROBIS_high_L$Q35, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q35 <- factor(sgbar.likert_ROBIS_high_L$Q35, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q41 <- as.factor(sgbar.likert_ROBIS_high_L$Q41)
sgbar.likert_ROBIS_high_L$Q41 <- revalue(sgbar.likert_ROBIS_high_L$Q41, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q41 <- factor(sgbar.likert_ROBIS_high_L$Q41, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q42 <- as.factor(sgbar.likert_ROBIS_high_L$Q42)
sgbar.likert_ROBIS_high_L$Q42 <- revalue(sgbar.likert_ROBIS_high_L$Q42, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q42 <- factor(sgbar.likert_ROBIS_high_L$Q42, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q43 <- as.factor(sgbar.likert_ROBIS_high_L$Q43)
sgbar.likert_ROBIS_high_L$Q43 <- revalue(sgbar.likert_ROBIS_high_L$Q43, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q43 <- factor(sgbar.likert_ROBIS_high_L$Q43, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q44 <- as.factor(sgbar.likert_ROBIS_high_L$Q44)
sgbar.likert_ROBIS_high_L$Q44 <- revalue(sgbar.likert_ROBIS_high_L$Q44, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q44 <- factor(sgbar.likert_ROBIS_high_L$Q44, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q45 <- as.factor(sgbar.likert_ROBIS_high_L$Q45)
sgbar.likert_ROBIS_high_L$Q45 <- revalue(sgbar.likert_ROBIS_high_L$Q45, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q45 <- factor(sgbar.likert_ROBIS_high_L$Q45, levels=desired.order, ordered=TRUE)

sgbar.likert_ROBIS_high_L$Q46 <- as.factor(sgbar.likert_ROBIS_high_L$Q46)
sgbar.likert_ROBIS_high_L$Q46 <- revalue(sgbar.likert_ROBIS_high_L$Q46, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
sgbar.likert_ROBIS_high_L$Q46 <- factor(sgbar.likert_ROBIS_high_L$Q46, levels=desired.order, ordered=TRUE)

attach(sgbar.likert_ROBIS_high_L)

########  Likert scale ROBIS high soubgroup --------------------

plot(likert(sgbar.likert_ROBIS_high_L[1:21]), main="Low RoB (n=12)", ylab="ROBIS items")
b<-likert(sgbar.likert_ROBIS_high[1:21])


########  ROBIS ALL --------------------

DB4_all<- DB4[,c(12:32)]
desired.order <- c("Low", "Probably low","Probably high", "High")
DB4_all$Q11 <- as.factor(DB4_all$Q11)
DB4_all$Q11 <- revalue(DB4_all$Q11, c("0"="Low", "1"="Probably low","2"="Probably high", "3"="High"))
DB4_all$Q11 <- factor(DB4_all$Q11, levels=desired.order, ordered=TRUE)

DB4_all$Q12 <- as.factor(DB4_all$Q12)
DB4_all$Q12 <- revalue(DB4_all$Q12, c("0"="Low","1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q12 <- factor(DB4_all$Q12, levels=desired.order, ordered=TRUE)

DB4_all$Q13 <- as.factor(DB4_all$Q13)
DB4_all$Q13 <- revalue(DB4_all$Q13, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q13 <- factor(DB4_all$Q13, levels=desired.order, ordered=TRUE)

DB4_all$Q14 <- as.factor(DB4_all$Q14)
DB4_all$Q14 <- revalue(DB4_all$Q14, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q14 <- factor(DB4_all$Q14, levels=desired.order, ordered=TRUE)

DB4_all$Q15 <- as.factor(DB4_all$Q15)
DB4_all$Q15 <- revalue(DB4_all$Q15, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q15 <- factor(DB4_all$Q15, levels=desired.order, ordered=TRUE)

DB4_all$Q21 <- as.factor(DB4_all$Q21)
DB4_all$Q21 <- revalue(DB4_all$Q21, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q21 <- factor(DB4_all$Q21, levels=desired.order, ordered=TRUE)

DB4_all$Q22 <- as.factor(DB4_all$Q22)
DB4_all$Q22 <- revalue(DB4_all$Q22, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q22 <- factor(DB4_all$Q22, levels=desired.order, ordered=TRUE)

DB4_all$Q23 <- as.factor(DB4_all$Q23)
DB4_all$Q23 <- revalue(DB4_all$Q23, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q23 <- factor(DB4_all$Q23, levels=desired.order, ordered=TRUE)

DB4_all$Q24 <- as.factor(DB4_all$Q24)
DB4_all$Q24 <- revalue(DB4_all$Q24, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q24 <- factor(DB4_all$Q24, levels=desired.order, ordered=TRUE)

DB4_all$Q25 <- as.factor(DB4_all$Q25)
DB4_all$Q25 <- revalue(DB4_all$Q25, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q25 <- factor(DB4_all$Q25, levels=desired.order, ordered=TRUE)

DB4_all$Q31 <- as.factor(DB4_all$Q31)
DB4_all$Q31 <- revalue(DB4_all$Q31, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q31 <- factor(DB4_all$Q31, levels=desired.order, ordered=TRUE)

DB4_all$Q32 <- as.factor(DB4_all$Q32)
DB4_all$Q32 <- revalue(DB4_all$Q32, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q32 <- factor(DB4_all$Q32, levels=desired.order, ordered=TRUE)

DB4_all$Q33 <- as.factor(DB4_all$Q33)
DB4_all$Q33 <- revalue(DB4_all$Q33, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q33 <- factor(DB4_all$Q33, levels=desired.order, ordered=TRUE)

DB4_all$Q34 <- as.factor(DB4_all$Q34)
DB4_all$Q34 <- revalue(DB4_all$Q34, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q34 <- factor(DB4_all$Q34, levels=desired.order, ordered=TRUE)

DB4_all$Q35 <- as.factor(DB4_all$Q35)
DB4_all$Q35 <- revalue(DB4_all$Q35, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q35 <- factor(DB4_all$Q35, levels=desired.order, ordered=TRUE)

DB4_all$Q41 <- as.factor(DB4_all$Q41)
DB4_all$Q41 <- revalue(DB4_all$Q41, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q41 <- factor(DB4_all$Q41, levels=desired.order, ordered=TRUE)

DB4_all$Q42 <- as.factor(DB4_all$Q42)
DB4_all$Q42 <- revalue(DB4_all$Q42, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q42 <- factor(DB4_all$Q42, levels=desired.order, ordered=TRUE)

DB4_all$Q43 <- as.factor(DB4_all$Q43)
DB4_all$Q43 <- revalue(DB4_all$Q43, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q43 <- factor(DB4_all$Q43, levels=desired.order, ordered=TRUE)

DB4_all$Q44 <- as.factor(DB4_all$Q44)
DB4_all$Q44 <- revalue(DB4_all$Q44, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q44 <- factor(DB4_all$Q44, levels=desired.order, ordered=TRUE)

DB4_all$Q45 <- as.factor(DB4_all$Q45)
DB4_all$Q45 <- revalue(DB4_all$Q45, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q45 <- factor(DB4_all$Q45, levels=desired.order, ordered=TRUE)

DB4_all$Q46 <- as.factor(DB4_all$Q46)
DB4_all$Q46 <- revalue(DB4_all$Q46, c("0"="Low", "1"="Probably low", "2"="Probably high", "3"="High"))
DB4_all$Q46 <- factor(DB4_all$Q46, levels=desired.order, ordered=TRUE)

attach(DB4_all)

##### likert scale ROBIS subgrupo high

plot(likert(DB4_all[1:21]), 
     main="All SRs (n=139)", 
     ylab="ROBIS items")
likert.bar.plot(likert(DB4_all[1:21]), 
                centered = FALSE, 
                main="All SRs (n=139)", 
                ylab="ROBIS items", 
                legend="", 
                legend.position = "right", 
                ordered=FALSE, 
                low.color="lightsalmon1", 
                high.color="skyblue3", 
                neutral.color="seagreen3")

########### with DIMENSIONS

file4<-read.csv2("plot_paper_ROBIS.csv", 
                 sep = ";", 
                 dec= ".", 
                 stringsAsFactors = TRUE,
                 header = TRUE)
DB4<-as.data.frame(file4)
DB4$D1_sum <- (DB4$Q11+DB4$Q12+DB4$Q13+DB4$Q14+DB4$Q15)
DB4$D2_sum <- (DB4$Q21+DB4$Q22+DB4$Q23+DB4$Q24+DB4$Q25)
DB4$D3_sum <- (DB4$Q31+DB4$Q32+DB4$Q33+DB4$Q34+DB4$Q35)
DB4$D4_sum <- (DB4$Q41+DB4$Q42+DB4$Q43+DB4$Q44+DB4$Q45+DB4$Q46)
DB4$radius<-(D1_sum+D2_sum+D3_sum+D4_sum)
desired.order <- c("High","Low", "Unclear")

DB4$D1 <- as.factor(DB4$D1)
DB4$D1 <- revalue(DB4$D1, c("LOW"="Low", "HIGH"="High", "UNCLEAR"="Unclear"))
DB4$D1 <- factor(DB4$D1, levels=desired.order, ordered=TRUE)

DB4$D2 <- as.factor(DB4$D2)
DB4$D2 <- revalue(DB4$D2, c("LOW"="Low", "HIGH"="High", "UNCLEAR"="Unclear"))
DB4$D2 <- factor(DB4$D2, levels=desired.order, ordered=TRUE)

DB4$D3 <- as.factor(DB4$D3)
DB4$D3 <- revalue(DB4$D3, c("LOW"="Low", "HIGH"="High", "UNCLEAR"="Unclear"))
DB4$D3 <- factor(DB4$D3, levels=desired.order, ordered=TRUE)

DB4$D4 <- as.factor(DB4$D4)
DB4$D4 <- revalue(DB4$D4, c("LOW"="Low", "HIGH"="High", "UNCLEAR"="Unclear"))
DB4$D4 <- factor(DB4$D4, levels=desired.order, ordered=TRUE)

DB4$RoB <- as.factor(DB4$RoB)
DB4$RoB <- revalue(DB4$RoB, c("LOW"="Low", "HIGH"="High", "UNCLEAR"="Unclear"))
DB4$RoB <- factor(DB4$RoB, levels=desired.order, ordered=TRUE)

names(DB4)[3] <- "1. Study elegibility criteria"
names(DB4)[4] <- "2. Identification and selection of studies"
names(DB4)[5] <- "3. Data collection and study appraisal"
names(DB4)[6] <- "4. Synthesis and findings"
names(DB4)[7] <- "Risk of Bias in the review"

likert.bar.plot(likert(DB4[3:7]), 
                centered = FALSE, 
                legend="", 
                legend.position = "right", 
                ordered=FALSE, 
                low.color="lightsalmon1", 
                high.color="skyblue3", 
                neutral.color="seagreen3")


############################################################################
######## 3: correlation matrix of ROBIS items ----------------

corr_robis_melt <- melt(cor(DB2[,4:24]))
corr_robis_melt$value <-round(corr_robis_melt$value, 2) 
corr_robis_melt<- subset(corr_robis_melt, value < 1) 
ggplot(corr_robis_melt, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label=value), size = 3, fontface = "bold") + 
  #scale_fill_gradient(low = "lightgreen", high = "red") +
  scale_fill_gradientn(colours=c("darkgreen","lightgreen","orange", "red"), na.value = "grey98", limits = c(-0.16, 1), breaks=seq(min(corr_robis_melt$value), max(corr_robis_melt$value), by=10)) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 12, face = "bold"))


############################################################################
######## 4: correlation network of ROBIS items ----------------

#### packages

library(ggraph)
library(igraph)

#### subsetting coef spearman >0.5

corr_robis_melt_0.6 <- subset(corr_robis_melt, value>0.5)

####  network plotting

graph <- graph_from_data_frame(corr_robis_melt_0.6, directed = FALSE)
p <- ggraph(graph, layout = 'kk') + 
    geom_edge_link() + 
    geom_node_point(show.legend = TRUE) + 
    ggtitle('') +
    geom_node_text(aes(label = name), repel = TRUE) +
    geom_edge_link(aes(edge_alpha = abs(value), edge_width = abs(value), color = value)) +
    guides(edge_alpha = "none", edge_width = "none") +
    scale_edge_colour_gradientn(limits = c(0.55, 0.9), colors = c("orange","firebrick2")) +
    geom_node_point(color = "darkblue", size = 3) 
p + theme_graph()
dodgerblue2


############################################################################
######### 5: Visualizing Multivariate Data with Radviz ---------------------


library ("Radviz")

#########  ROBIS signaling questions --------------------

######## normalizing the data
norm <- apply(DB2[,4:24],2,do.L,fun=function(x) quantile(x,c(0.025,0.975)))

######## defining the anchors
ct.S <- make.S(dimnames(DB2[,4:24])[[2]])

######## computing the similarity matrix (distances between columns and rows)
ct.sim <- cosine(norm)

######## getting the current radviz-independent measure of projection efficiency
in.da(ct.S,ct.sim)

######## getting the current radviz-dependent measure of projection efficiency
rv.da(ct.S,ct.sim)

######## The radviz-independent score should be maximal when the dimensional anchor positions are optimal. 
######## The radviz-dependent score should be minimal when the dimensional anchor positions are optimal.
######## Optimization procedure:

optim.ct <- do.optim(ct.S,ct.sim,iter=100,n=1000)
ct.S <- make.S(tail(optim.ct$best,1)[[1]])

######## Getting final projections ----------------
######## The do.radviz function will then use the normalized values and the Springs 
######## to project each thing in a 2D space:
ct.rv <- do.radviz(norm,ct.S)

######## Visualizing the results  ----------------
########  There is a S3 plot function defined for radviz; using the default will give the following result

#### AMSTAR
sub.rv_low_quality <- subset(ct.rv,DB2$AMSTAR_levels_2=="low_quality")
sub.rv_moderate_quality <- subset(ct.rv,DB2$AMSTAR_levels_2=="moderate_quality")
sub.rv_high_quality <- subset(ct.rv,DB2$AMSTAR_levels_2=="high_quality")
plot(ct.rv,point.shape=19, point.color=c("coral", "seagreen3", "dodgerblue2")[DB2$AMSTAR_levels_2])
contour(sub.rv_low_quality, add=T, contour.color = "dodgerblue2")
contour(sub.rv_moderate_quality, add=T, contour.color = "seagreen3")
contour(sub.rv_high_quality, add=T, contour.color = "coral")

#### ROBIS
sub.rv_low <- subset(ct.rv,DB2$RoB_ROBIS=="LOW")
sub.rv_high <- subset(ct.rv,DB2$RoB_ROBIS=="HIGH")
plot(ct.rv,point.shape=19, point.color=c("coral", "turquoise3")[DB2$RoB_ROBIS])
contour(sub.rv_low, add=T, contour.color = "turquoise3")
contour(sub.rv_high, add=T, contour.color = "coral")

####  ROBIS domains
####  normalizing the data
DB2$D1<-DB2[,4]+DB2[,5]+DB2[,6]+DB2[,7]+DB2[,8]
DB2$D2<-DB2[,9]+DB2[,10]+DB2[,11]+DB2[,12]+DB2[,13]
DB2$D3<-DB2[,14]+DB2[,15]+DB2[,16]+DB2[,17]+DB2[,18]
DB2$D4<-DB2[,19]+DB2[,20]+DB2[,21]+DB2[,22]+DB2[,23]+DB2[,24]
attach(DB2)
norm <- apply(DB2[,43:46],2,do.L,fun=function(x) quantile(x,c(0.025,0.975)))

ct.S <- make.S(dimnames(DB2[,43:46])[[2]])
ct.sim <- cosine(norm)
in.da(ct.S,ct.sim)
rv.da(ct.S,ct.sim)
optim.ct <- do.optim(ct.S,ct.sim,iter=100,n=1000)
ct.S <- make.S(tail(optim.ct$best,1)[[1]])
ct.rv <- do.radviz(norm,ct.S)

# Visualizing the results:
## ROBIS
sub.rv_low <- subset(ct.rv,DB2$RoB_ROBIS=="LOW")
sub.rv_high <- subset(ct.rv,DB2$RoB_ROBIS=="HIGH")
plot(ct.rv,point.shape=19, point.color=c("coral", "turquoise3")[DB2$RoB_ROBIS])
contour(sub.rv_low, add=T, contour.color = "turquoise3")
contour(sub.rv_high, add=T, contour.color = "coral")

## AMSTAR
sub.rv_low_quality <- subset(ct.rv,DB2$AMSTAR_levels_2=="low_quality")
sub.rv_moderate_quality <- subset(ct.rv,DB2$AMSTAR_levels_2=="moderate_quality")
sub.rv_high_quality <- subset(ct.rv,DB2$AMSTAR_levels_2=="high_quality")
plot(ct.rv,point.shape=19, point.color=c("coral", "seagreen3", "dodgerblue2")[DB2$AMSTAR_levels_2])
contour(sub.rv_low_quality, add=T, contour.color = "dodgerblue2")
contour(sub.rv_moderate_quality, add=T, contour.color = "seagreen3")
contour(sub.rv_high_quality, add=T, contour.color = "coral")



