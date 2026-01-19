library(readxl)
library(writexl)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(xlsx)


####HEATMAP####
#### filtering & formatting ####
setwd("C:/files/CHH_pan20_orthogroup/")
data <- read_excel("pan20_orthogroup_CHH.xlsx")
clean_data <- data[,c(-2,-3)] 
filtered_data <- clean_data %>%
  filter(rowSums(is.na(.)) <= 2)
print(filtered_data)
write_xlsx(filtered_data, "pan20_merged_plot_2NA_CHH.xlsx")

data <- read_excel("pan20_merged_plot_2NA_CHH.xlsx")
data[,-1] <- sapply(data[,-1], as.numeric)
data_matrix <- as.matrix(data[,-1])
head(data_matrix)
dim(data_matrix)


#### creating complex heatmap ####
colors <- viridis(100, option="D", direction=-1)
colors <- cividis(100, direction=-1)
pdf("pan20_region_cluster_CHH_noraster_civ_5_5.pdf", width = 5, height = 5)
heatmap_DMO_pan20_all_CHH <- Heatmap(data_matrix, name = "CHH Methylation Level %",
                                     row_title = "Orthogroup", column_title = "Accession",
                                     col = colors, na_col = "grey",
                                     cluster_rows = TRUE, cluster_columns = TRUE,
                                     show_row_names = FALSE, show_column_names = TRUE,
                                     use_raster = FALSE
)
plot(heatmap_DMO_pan20_all_CHH)
dev.off()

browseURL("C:/files/CHH_pan20_orthogroup/pan20_region_cluster_CHH_noraster_civ_5_5.pdf")

####MANHATTAN PLOT : compare CHG meth. across Asian and non-Asian groups####

library(readxl)
library(writexl)

setwd("C:/files/CHH_pan20_orthogroup/")
data <- read_xlsx("pan20_merged_plot_2NA_CHH.xlsx")
asian_group <- c("FT11", "Du_Li_Huang", "HOR_7552", "OUN333", "Chi_Ba_Damai", "Akashinriki", "HOR_21599")
non_asian_group <- c("Barke", "Golden_Promise", "Hockett", "HOR_3081", "HOR_3365", "HOR_8148", "HOR_9043", "HOR_10350", "HOR_13821", "HOR_13942", "Igri", "Morex", "RGT_Planet")
head(data)

calculate_stats <- function(row) {
  asian_values <- as.numeric(row[asian_group])
  non_asian_values <- as.numeric(row[non_asian_group])
  t_test <- t.test(asian_values, non_asian_values)
  p_value <- t_test$p.value
  fold_change <- mean(asian_values, na.rm = TRUE) / mean(non_asian_values, na.rm = TRUE)
  return(c(p_value, fold_change))
}

results <- t(data.frame(apply(data[,-1], 1, calculate_stats)))       #remove col 1(HOG_ID) and calculate t & p
head(results)
colnames(results) <- c("p_value", "fold_change")

final_results <- cbind(data[, 1], results)  #give HoG_ID to each t & p
colnames(final_results)[1] <- "HOG"
head(final_results)

write_xlsx(final_results, "pan20_orthogroup_CHH_2NA_ttest_AsianVsNonAsian.xlsx")

print(final_results)


##combine Morex gene ID and coordinates (TSS)##
library(readxl)
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_AsianVsNonAsian.xlsx")
file2 <- read_excel("pan20_single_copy_core.xlsx")
file2_sub <- file2[, c("HOG", "Morex")]
merged <- merge(file1, file2_sub, by = "HOG", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_AsianVsNonAsian.xlsx")

##Trimming out transcript .1 , .2 and keeping gene name only##
merged[,4] <- sub("\\.\\d+$", "", merged[,4])       #inclusive removal of . and after . from 4th Column
head(merged)
write_xlsx(merged,"pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_AsianVsNonAsian.xlsx" )

##three orthologs on contigs were removed, transcripts '.1' and '.2' were removed###
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_AsianVsNonAsian.xlsx")
file2 <- read.csv("Morex_genes_tss.txt", header = F, sep = " ")
file2 <- as.data.frame(file2)
head(file1)
head(file2)
colnames(file2) = c("CHR", "Morex", "TSS")
tail(file2)
merged <- merge(file1, file2, by = "Morex", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_MorexTSS_AsianVsNonAsian.xlsx")
##three orthologs on contigs were removed for plotting##


###plot p value
library(qqman)
library(ggplot2)
df <- read_excel("pan20_orthogroup_CHH_2NA_ttest_MorexTSS_AsianVsNonAsian.xlsx")

#Keep only standard barley chromosomes (chr1H-chr7H)
tail(df)
df <- df[grepl("^chr[1-7]H$", df$CHR), ]
tail(df)
df$CHR <- as.numeric(gsub("chr|H", "", df$CHR))   #subsitute chr_H with blank so CHR is all number
tail(df)
any(is.na(df$CHR))


manhattan_df <- data.frame(
  SNP = df$HOG,
  CHR = df$CHR,
  BP  = df$TSS,
  P   = df$p_value
)



tail(manhattan_df)
pdf("CHH_AsianVsNonAsian_DMO.pdf", width = 8, height = 4)
manhattan(manhattan_df,
          main = "Differentially CHH-methylated orthologs Asian vs Non-Asian",
          ylim = c(0, max(-log10(manhattan_df$P), na.rm=TRUE) + 1),
          genomewideline = FALSE,
          suggestiveline = FALSE,    
          cex = 1, cex.axis = 1.2,
          col = c("#1e32b4", "grey"))
abline(h = -log10(0.05/22783), col = "black", lty = 2)          

dev.off()


####MANHATTAN PLOT: CHH Cluster 1 vs Cluster 2####

data <- read_xlsx("pan20_merged_plot_2NA_CHH.xlsx")
cluster1 <- c("Barke", "Golden_Promise", "Hockett", "HOR_3081", "HOR_3365", "HOR_8148", "HOR_9043", "HOR_10350", "HOR_13821", "HOR_13942", "HOR_21599","FT11", "Du_Li_Huang", "HOR_7552", "OUN333", "Igri", "Morex", "RGT_Planet")
cluster2 <- c("Chi_Ba_Damai", "Akashinriki")
head(data)

calculate_stats <- function(row) {
  cluster1 <- as.numeric(row[cluster1])
  cluster2 <- as.numeric(row[cluster2])
  t_test <- t.test(cluster1, cluster2)
  p_value <- t_test$p.value
  fold_change <- mean(cluster1, na.rm = TRUE) / mean(cluster2, na.rm = TRUE)
  return(c(p_value, fold_change))
}

results <- t(data.frame(apply(data[,-1], 1, calculate_stats)))       #remove col 1(HOG_ID) and calculate t & p
head(results)
colnames(results) <- c("p_value", "fold_change")

final_results <- cbind(data[, 1], results)  #give HoG_ID to each t & p
colnames(final_results)[1] <- "HOG"
head(final_results)

write_xlsx(final_results, "pan20_orthogroup_CHH_2NA_ttest_cluster1vs2.xlsx")

print(final_results)


##combine Morex gene ID and coordinates (TSS)##
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_cluster1vs2.xlsx")
file2 <- read_excel("pan20_single_copy_core.xlsx")
file2_sub <- file2[, c("HOG", "Morex")]
merged <- merge(file1, file2_sub, by = "HOG", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_cluster1vs2.xlsx")

##Trimming out transcript .1 , .2 and keeping gene name only##
merged[,4] <- sub("\\.\\d+$", "", merged[,4])       #inclusive removal of . and after . from 4th Column
head(merged)
write_xlsx(merged,"pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_cluster1vs2.xlsx" )

##three orthologs on contigs were removed, transcripts '.1' and '.2' were removed###
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_cluster1vs2.xlsx")
file2 <- read.csv("Morex_genes_tss.txt", header = F, sep = " ")
file2 <- as.data.frame(file2)
head(file1)
head(file2)
colnames(file2) = c("CHR", "Morex", "TSS")
tail(file2)
merged <- merge(file1, file2, by = "Morex", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_MorexTSS_cluster1vs2.xlsx")
##three orthologs on contigs were removed for plotting##


###plot p value
library(qqman)
library(ggplot2)
df <- read_excel("pan20_orthogroup_CHH_2NA_ttest_MorexTSS_cluster1vs2.xlsx")

#Keep only standard barley chromosomes (chr1H-chr7H)
tail(df)
df <- df[grepl("^chr[1-7]H$", df$CHR), ]
tail(df)
df$CHR <- as.numeric(gsub("chr|H", "", df$CHR))   #subsitute chr_H with blank so CHR is all number
tail(df)
any(is.na(df$CHR))


manhattan_df <- data.frame(
  SNP = df$HOG,
  CHR = df$CHR,
  BP  = df$TSS,
  P   = df$p_value
)


tail(manhattan_df)
pdf("CHH_Cluster1vs2_DMO.pdf", width = 8, height = 4)
manhattan(manhattan_df,
          main = "Differentially CHH methylated orthologs between cluster 1 & 2 groups",
          ylim = c(0, max(-log10(manhattan_df$P), na.rm=TRUE) + 1),
          genomewideline = FALSE,
          suggestiveline = FALSE,    
          cex = 1, cex.axis = 1.2,
          col = c("#1e32b4", "grey"))
abline(h = -log10(0.05/22783), col = "black", lty = 2)          

dev.off()


####MANHATTAN PLOT : compare CHG-methylation between 2-rowed and 6-rowed groups####

library(readxl)
library(writexl)

setwd("C:/files/CHH_pan20_orthogroup/")
data <- read_xlsx("pan20_merged_plot_2NA_CHH.xlsx")
two_row_group <- c("Barke", "Igri", "RGT_Planet", "FT11", "HOR_21599", "Hockett", "Golden_Promise", "HOR_13821", "HOR_8148")
six_row_group <- c("Morex", "Akashinriki", "HOR_3081", "OUN333", "Du_Li_Huang", "HOR_3365", "Chi_Ba_Damai", "HOR_10350", "HOR_13942" ,"HOR_7552", "HOR_9043")
head(data)

calculate_stats <- function(row) {
  two_row_values <- as.numeric(row[two_row_group])
  six_row_values <- as.numeric(row[six_row_group])
  t_test <- t.test(two_row_values, six_row_values)
  p_value <- t_test$p.value
  fold_change <- mean(two_row_values, na.rm = TRUE) / mean(six_row_values, na.rm = TRUE)
  return(c(p_value, fold_change))
}

results <- t(data.frame(apply(data[,-1], 1, calculate_stats)))       #remove col 1(HOG_ID) and calculate t & p
head(results)
colnames(results) <- c("p_value", "fold_change")

final_results <- cbind(data[, 1], results)  #give HoG_ID to each t & p
colnames(final_results)[1] <- "HOG"
head(final_results)

write_xlsx(final_results, "pan20_orthogroup_CHH_2NA_ttest_row_type.xlsx")

print(final_results)


##combine Morex gene ID and coordinates (TSS)##
library(readxl)
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_row_type.xlsx")
file2 <- read_excel("pan20_single_copy_core.xlsx")
file2_sub <- file2[, c("HOG", "Morex")]
merged <- merge(file1, file2_sub, by = "HOG", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_row_type.xlsx")

##Trimming out transcript .1 , .2 and keeping gene name only##
merged[,4] <- sub("\\.\\d+$", "", merged[,4])       #inclusive removal of . and after . from 4th Column
head(merged)
write_xlsx(merged,"pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_row_type.xlsx" )

##three orthologs on contigs were removed, transcripts '.1' and '.2' were removed###
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_row_type.xlsx")
file2 <- read.csv("Morex_genes_tss.txt", header = F, sep = " ")
file2 <- as.data.frame(file2)
head(file1)
head(file2)
colnames(file2) = c("CHR", "Morex", "TSS")
tail(file2)
merged <- merge(file1, file2, by = "Morex", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_MorexTSS_row_type.xlsx")
##three orthologs on contigs were removed for plotting##


###plot p value
library(qqman)
library(ggplot2)
df <- read_excel("pan20_orthogroup_CHH_2NA_ttest_MorexTSS_row_type.xlsx")

#Keep only standard barley chromosomes (chr1H-chr7H)
tail(df)
df <- df[grepl("^chr[1-7]H$", df$CHR), ]
tail(df)
df$CHR <- as.numeric(gsub("chr|H", "", df$CHR))   #subsitute chr_H with blank so CHR is all number
tail(df)
any(is.na(df$CHR))


manhattan_df <- data.frame(
  SNP = df$HOG,
  CHR = df$CHR,
  BP  = df$TSS,
  P   = df$p_value
)



tail(manhattan_df)
pdf("CHH_twoVSsix_rowtype_DMO.pdf", width = 8, height = 4)
manhattan(manhattan_df,
          main = "Differentially CHH-methylated orthologs Two vs Six Row Type",
          ylim = c(0, max(-log10(manhattan_df$P), na.rm=TRUE) + 1),
          genomewideline = FALSE,
          suggestiveline = FALSE,    
          cex = 1, cex.axis = 1.2,
          col = c("#1e32b4", "grey"))
abline(h = -log10(0.05/22783), col = "black", lty = 2)          
dev.off()



####MANHATTAN PLOT : compare CHG-methylation between spring and winter groups####

setwd("C:/files/CHH_pan20_orthogroup/")
data <- read_xlsx("pan20_merged_plot_2NA_CHH.xlsx")
spring_group <- c("Morex", "Barke", "RGT_Planet", "Du_Li_Huang", "FT11", "OUN333", "Hockett", "Chi_Ba_Damai", "Golden_Promise", "HOR_10350", "HOR_13821", "HOR_13942" ,"HOR_7552", "HOR_8148", "HOR_9043")
winter_group <- c("Igri", "Akashinriki", "HOR_3081", "HOR_21599", "HOR_3365")
head(data)

calculate_stats <- function(row) {
  spring_values <- as.numeric(row[spring_group])
  winter_values <- as.numeric(row[winter_group])
  t_test <- t.test(spring_values, winter_values)
  p_value <- t_test$p.value
  fold_change <- mean(spring_values, na.rm = TRUE) / mean(winter_values, na.rm = TRUE)
  return(c(p_value, fold_change))
}

results <- t(data.frame(apply(data[,-1], 1, calculate_stats)))       #remove col 1(HOG_ID) and calculate t & p
head(results)
colnames(results) <- c("p_value", "fold_change")

final_results <- cbind(data[, 1], results)  #give HoG_ID to each t & p
colnames(final_results)[1] <- "HOG"
head(final_results)

write_xlsx(final_results, "pan20_orthogroup_CHH_2NA_ttest_springVSwinter.xlsx")

print(final_results)


##combine Morex gene ID and coordinates (TSS)##
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_springVSwinter.xlsx")
file2 <- read_excel("pan20_single_copy_core.xlsx")
file2_sub <- file2[, c("HOG", "Morex")]
merged <- merge(file1, file2_sub, by = "HOG", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_springVSwinter.xlsx")

##Trimming out transcript .1 , .2 and keeping gene name only##
merged[,4] <- sub("\\.\\d+$", "", merged[,4])       #inclusive removal of . and after . from 4th Column
head(merged)
write_xlsx(merged,"pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_springVSwinter.xlsx" )

##three orthologs on contigs were removed, transcripts '.1' and '.2' were removed###
file1 <- read_excel("pan20_orthogroup_CHH_2NA_ttest_GeneID_Morex_trimmed_springVSwinter.xlsx")
file2 <- read.csv("Morex_genes_tss.txt", header = F, sep = " ")
file2 <- as.data.frame(file2)
head(file1)
head(file2)
colnames(file2) = c("CHR", "Morex", "TSS")
tail(file2)
merged <- merge(file1, file2, by = "Morex", all.x = TRUE)
head(merged)
write_xlsx(merged, "pan20_orthogroup_CHH_2NA_ttest_MorexTSS_springVSwinter.xlsx")
##three orthologs on contigs were removed for plotting##


###plot p value
library(qqman)
library(ggplot2)
df <- read_excel("pan20_orthogroup_CHH_2NA_ttest_MorexTSS_springVSwinter.xlsx")

#Keep only standard barley chromosomes (chr1H-chr7H)
tail(df)
df <- df[grepl("^chr[1-7]H$", df$CHR), ]
tail(df)
df$CHR <- as.numeric(gsub("chr|H", "", df$CHR))   #subsitute chr_H with blank so CHR is all number
tail(df)
any(is.na(df$CHR))


manhattan_df <- data.frame(
  SNP = df$HOG,
  CHR = df$CHR,
  BP  = df$TSS,
  P   = df$p_value
)


tail(manhattan_df)
pdf("CHH_springVSwinter_growth_DMO.pdf", width = 8, height = 4)
manhattan(manhattan_df,
          main = "Differentially CHH-methylated orthologs Spring Vs Winter Growth",
          ylim = c(0, max(-log10(manhattan_df$P), na.rm=TRUE) + 1),
          genomewideline = FALSE,
          suggestiveline = FALSE,    
          cex = 1, cex.axis = 1.2,
          col = c("#1e32b4", "grey"))
abline(h = -log10(0.05/22783), col = "black", lty = 2)          
dev.off()



