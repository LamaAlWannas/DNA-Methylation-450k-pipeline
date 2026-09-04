## ----setup, include=FALSE--------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)


## ----config----------------------------------------------------------------------------------
rm(list=ls())
suppressMessages(library(minfi))
address       <- 18743318      
detP_threshold <- 0.01         



## ----load_data-------------------------------------------------------------------------------
# Read the sample sheet
SampleSheet <- read.csv("Input_Data/SampleSheet_Report_II.csv", header=T, stringsAsFactors=T)
SampleSheet

# minfi reads the sample sheet and builds the path to the IDAT files
baseDir <- "Input_Data"
targets <- read.metharray.sheet(baseDir)

# Our IDAT files are named SampleID_Sentrix_ID_Sentrix_Position
# (e.g. GSM5319592_200121140049_R01C02_Grn.idat), so we build the Basename
# (the path without the _Grn.idat / _Red.idat ending) accordingly:
targets$Basename <- file.path(baseDir, paste(targets$SampleID, targets$Slide, targets$Array, sep="_")) 

# Create the RGChannelSet object
RGset <- read.metharray.exp(targets = targets)
save(RGset, file="RGset.RData")
RGset


## ----red_green-------------------------------------------------------------------------------
Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)

Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)


## ----address---------------------------------------------------------------------------------
# Red and Green fluorescence for the address assigned to the group
Red[rownames(Red)==address, ]
Green[rownames(Green)==address, ]


## ----address_optional------------------------------------------------------------------------
# Optional: check in the manifest whether the address belongs to a Type I or
# Type II probe and, for Type I probes, report the colour.
load("Illumina450Manifest_clean.RData")
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID==address |
                          Illumina450Manifest_clean$AddressB_ID==address,
                          c("IlmnID","AddressA_ID","AddressB_ID","Infinium_Design_Type","Color_Channel")]


## ----mset------------------------------------------------------------------------------------
MSet.raw <- preprocessRaw(RGset)
save(MSet.raw, file="MSet_raw.RData")
MSet.raw


## ----qcplot----------------------------------------------------------------------------------
qc <- getQC(MSet.raw)
plotQC(qc)


## ----controls--------------------------------------------------------------------------------
controlStripPlot(RGset, controls="NEGATIVE")


## ----detP------------------------------------------------------------------------------------
detP <- detectionP(RGset)
save(detP, file="detP.RData")

# How many probes fail the detection p-value threshold in each sample?
failed <- detP > detP_threshold
table(failed)

# number of failed positions per sample
n_failed_per_sample <- colSums(failed)
n_failed_per_sample


## ----beta_M_raw------------------------------------------------------------------------------
beta <- getBeta(MSet.raw)
M    <- getM(MSet.raw)

# subset the matrices to retain CTRL or DIS subjects
beta_CTRL <- beta[, SampleSheet$Group=="CTRL"]
beta_DIS  <- beta[, SampleSheet$Group=="DIS"]
M_CTRL    <- M[, SampleSheet$Group=="CTRL"]
M_DIS     <- M[, SampleSheet$Group=="DIS"]

# mean methylation value of each probe, in the two groups
mean_beta_CTRL <- apply(beta_CTRL, 1, mean, na.rm=T)
mean_beta_DIS  <- apply(beta_DIS,  1, mean, na.rm=T)
mean_M_CTRL    <- apply(M_CTRL,    1, mean, na.rm=T)
mean_M_DIS     <- apply(M_DIS,     1, mean, na.rm=T)

d_mean_beta_CTRL <- density(mean_beta_CTRL, na.rm=T)
d_mean_beta_DIS  <- density(mean_beta_DIS,  na.rm=T)
d_mean_M_CTRL    <- density(mean_M_CTRL,    na.rm=T)
d_mean_M_DIS     <- density(mean_M_DIS,     na.rm=T)

par(mfrow=c(1,2))
plot(d_mean_beta_CTRL, main="Density of mean Beta", col="orange")
lines(d_mean_beta_DIS, col="purple")
legend("top", legend=c("CTRL","DIS"), col=c("orange","purple"), lty=1)

plot(d_mean_M_CTRL, main="Density of mean M", col="orange")
lines(d_mean_M_DIS, col="purple")
legend("top", legend=c("CTRL","DIS"), col=c("orange","purple"), lty=1)


## ----chemistry-------------------------------------------------------------------------------
# Split the probes according to their chemistry (Type I / Type II)
dfI  <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI  <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)


## ----normalize-------------------------------------------------------------------------------
set.seed(123)   # SWAN has a random component: fix the seed for reproducibility
RGset_normalized <- preprocessSWAN(RGset)         
beta_normalized <- getBeta(RGset_normalized)
save(beta_normalized, file="beta_normalized.RData")


## ----six_panel, fig.width=15, fig.height=8---------------------------------------------------
# RAW beta, split by chemistry
beta_I  <- beta[rownames(beta) %in% dfI$IlmnID, ]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]
mean_beta_I  <- apply(beta_I, 1, mean, na.rm=T)
mean_beta_II <- apply(beta_II,1, mean, na.rm=T)
sd_beta_I    <- apply(beta_I, 1, sd,   na.rm=T)
sd_beta_II   <- apply(beta_II,1, sd,   na.rm=T)
d_mean_beta_I  <- density(mean_beta_I,  na.rm=T)
d_mean_beta_II <- density(mean_beta_II, na.rm=T)
d_sd_beta_I    <- density(sd_beta_I,    na.rm=T)
d_sd_beta_II   <- density(sd_beta_II,   na.rm=T)

# NORMALIZED beta, split by chemistry
beta_norm_I  <- beta_normalized[rownames(beta_normalized) %in% dfI$IlmnID, ]
beta_norm_II <- beta_normalized[rownames(beta_normalized) %in% dfII$IlmnID, ]
mean_beta_norm_I  <- apply(beta_norm_I, 1, mean, na.rm=T)
mean_beta_norm_II <- apply(beta_norm_II,1, mean, na.rm=T)
sd_beta_norm_I    <- apply(beta_norm_I, 1, sd,   na.rm=T)
sd_beta_norm_II   <- apply(beta_norm_II,1, sd,   na.rm=T)
d_mean_beta_norm_I  <- density(mean_beta_norm_I,  na.rm=T)
d_mean_beta_norm_II <- density(mean_beta_norm_II, na.rm=T)
d_sd_beta_norm_I    <- density(sd_beta_norm_I,    na.rm=T)
d_sd_beta_norm_II   <- density(sd_beta_norm_II,   na.rm=T)

par(mfrow=c(2,3))
plot(d_mean_beta_I, col="blue", main="raw beta")
lines(d_mean_beta_II, col="red")
plot(d_sd_beta_I, col="blue", main="raw sd")
lines(d_sd_beta_II, col="red")
boxplot(beta, main="raw beta")

plot(d_mean_beta_norm_I, col="blue", main="normalized beta")
lines(d_mean_beta_norm_II, col="red")
plot(d_sd_beta_norm_I, col="blue", main="normalized sd")
lines(d_sd_beta_norm_II, col="red")
boxplot(beta_normalized, main="normalized beta")


## ----boxplot_by_group, fig.width=12, fig.height=5--------------------------------------------
# Optional: colour the boxplots according to the group (CTRL / DIS)
palette(c("orange","purple"))   # CTRL, DIS  (levels are alphabetical: CTRL, DIS)
par(mfrow=c(1,2))
boxplot(beta,            col=SampleSheet$Group, main="raw beta - by group",        ylim=c(0,1))
boxplot(beta_normalized, col=SampleSheet$Group, main="normalized beta - by group", ylim=c(0,1))


## ----pca-------------------------------------------------------------------------------------
pca_results <- prcomp(t(beta_normalized), scale.=T)
print(summary(pca_results))


## ----pca_plots, fig.width=15, fig.height=5---------------------------------------------------
par(mfrow=c(1,3))

# by Group
palette(c("orange","purple"))
plot(pca_results$x[,1], pca_results$x[,2], cex=2, pch=2, col=SampleSheet$Group,
     xlab="PC1", ylab="PC2", main="PCA - Group")
text(pca_results$x[,1], pca_results$x[,2], labels=rownames(pca_results$x), cex=0.5, pos=1)
legend("bottomright", legend=levels(SampleSheet$Group), col=1:nlevels(SampleSheet$Group), pch=2)

# by Sex
palette(c("pink","blue"))
plot(pca_results$x[,1], pca_results$x[,2], cex=2, pch=2, col=SampleSheet$Sex,
     xlab="PC1", ylab="PC2", main="PCA - Sex")
text(pca_results$x[,1], pca_results$x[,2], labels=rownames(pca_results$x), cex=0.5, pos=1)
legend("bottomright", legend=levels(SampleSheet$Sex), col=1:nlevels(SampleSheet$Sex), pch=2)

# by batch (Sentrix_ID)
SampleSheet$Sentrix_ID <- factor(SampleSheet$Sentrix_ID)
palette(rainbow(nlevels(SampleSheet$Sentrix_ID)))
plot(pca_results$x[,1], pca_results$x[,2], cex=2, pch=2, col=SampleSheet$Sentrix_ID,
     xlab="PC1", ylab="PC2", main="PCA - Batch (Sentrix_ID)")
text(pca_results$x[,1], pca_results$x[,2], labels=rownames(pca_results$x), cex=0.5, pos=1)
legend("bottomright", legend=levels(SampleSheet$Sentrix_ID), col=1:nlevels(SampleSheet$Sentrix_ID), pch=2)


## ----dm_test---------------------------------------------------------------------------------

dmp <- dmpFinder(beta_normalized, pheno=SampleSheet$Group, type="categorical")
dmp <- dmp[rownames(beta_normalized), ]        
pValues <- dmp$pval
names(pValues) <- rownames(beta_normalized)

# table of probes with their nominal p-value
final <- data.frame(beta_normalized, pValue = pValues)
final <- final[order(final$pValue), ]
head(final)


## ----correction------------------------------------------------------------------------------
corrected_pValues_BH   <- p.adjust(final$pValue, "BH")
corrected_pValues_Bonf <- p.adjust(final$pValue, "bonferroni")
final <- data.frame(final, corrected_pValues_BH, corrected_pValues_Bonf)

# how many differentially methylated probes (threshold 0.05)?
nominal    <- nrow(final[final$pValue <= 0.05, ])
bonferroni <- nrow(final[final$corrected_pValues_Bonf <= 0.05, ])
BH         <- nrow(final[final$corrected_pValues_BH <= 0.05, ])

data.frame(nominal, Bonferroni=bonferroni, BH)


## ----volcano---------------------------------------------------------------------------------
# delta = mean(DIS) - mean(CTRL), aligned to the rows of beta_normalized
mean_beta_norm_CTRL <- apply(beta_normalized[, SampleSheet$Group=="CTRL"], 1, mean, na.rm=T)
mean_beta_norm_DIS  <- apply(beta_normalized[, SampleSheet$Group=="DIS"],  1, mean, na.rm=T)
delta <- mean_beta_norm_DIS - mean_beta_norm_CTRL

toVolcPlot <- data.frame(delta, -log10(pValues))

plot(toVolcPlot[,1], toVolcPlot[,2], pch=16, cex=0.5,
     xlab="delta (DIS - CTRL)", ylab="-log10(p-value)", main="Volcano plot")
abline(h=-log10(0.01), col="red")
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)), ]
points(toHighlight[,1], toHighlight[,2], pch=16, cex=0.7, col="red")


## ----manhattan, fig.width=12, fig.height=6---------------------------------------------------
library(qqman)

# annotate the results with the manifest (CHR and MAPINFO)
final_annot <- data.frame(IlmnID=rownames(final), final)
final_annot <- merge(final_annot, Illumina450Manifest_clean, by="IlmnID")

input_Manhattan <- final_annot[, c("IlmnID","CHR","MAPINFO","pValue")]

# put the chromosomes in the right order and convert to numeric (as in class)
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14",
               "15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR, levels=order_chr)
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)

manhattan(input_Manhattan, snp="IlmnID", chr="CHR", bp="MAPINFO", p="pValue",
          col=rainbow(24))


## ----heatmap, fig.width=9, fig.height=9------------------------------------------------------
library(gplots)

# the beta values are the first columns of "final" (one per sample)
n_samples <- ncol(beta_normalized)
input_heatmap <- as.matrix(final[1:100, 1:n_samples])

# colour bar for the samples, according to the group (CTRL / DIS)
colorbar <- ifelse(SampleSheet$Group=="CTRL", "green", "orange")

# green-black-red palette, as in the lesson's final heatmap
col2 <- colorRampPalette(c("green", "black", "red"))(100)

heatmap.2(input_heatmap, col=col2, Rowv=T, Colv=T,
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram="both", key=T, ColSideColors=colorbar,
          density.info="none", trace="none", scale="none", symm=F,
          main="Top 100 differentially methylated probes - Average linkage")

legend("topright", legend=c("CTRL","DIS"),
       fill=c("green","orange"), border="black", bty="n", cex=0.8)

