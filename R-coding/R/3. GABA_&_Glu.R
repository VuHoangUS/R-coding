# Loading the appropriate libraries
library(datasets)
library(ggplot2)
library(multcompView)
library(dplyr)
if (!require(readxl)) install.packages("readxl")
library(readxl)

GABAergic <- read_xlsx("D:/baitapR/means_data.xlsx", sheet = "GABAergic")
Glutamatergic <- read_xlsx("D:/baitapR/means_data.xlsx", sheet = "Glutamatergic")
Nonneuronal <- read_xlsx("D:/baitapR/means_data.xlsx", sheet = "Nonneuronal")

# Get data of 5 research genes
## NNAT
NNAT_GABA <- GABAergic[34884, ]
NNAT_Glu <- Glutamatergic[34884, ]
NNAT_Non <- Nonneuronal[34884, ]

## STMN1
STMN1_GABA <- GABAergic[45460, ]
STMN1_Glu <- Glutamatergic[45460, ]
STMN1_Non <- Nonneuronal[45460, ]

## STMN2
STMN2_GABA <- GABAergic[45463, ]
STMN2_Glu <- Glutamatergic[45463, ]
STMN2_Non <- Nonneuronal[45463, ]

## BASP1
BASP1_GABA <- GABAergic[1934, ]
BASP1_Glu <- Glutamatergic[1934, ]
BASP1_Non <- Nonneuronal[1934, ]

## ITM2B
ITM2B_GABA <- GABAergic[11758, ]
ITM2B_Glu <- Glutamatergic[11758, ]
ITM2B_Non <- Nonneuronal[11758, ]


library(tidyverse)
# Convert from row to column 
NNAT_GABA <- NNAT_GABA %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values") 
NNAT_GABA$feature <- gsub("NNAT", "GABAergic", NNAT_GABA$feature)
NNAT_Glu <- NNAT_Glu %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
NNAT_Glu$feature <- gsub("NNAT", "Glutamatergic", NNAT_Glu$feature)
NNAT_Non <- NNAT_Non %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
NNAT_Non$feature <- gsub("NNAT", "Non-neuronal", NNAT_Non$feature)
# Concatenate data frames
NNAT <- full_join(full_join(NNAT_GABA, NNAT_Glu),NNAT_Non)


# Convert from row to column 
STMN1_GABA <- STMN1_GABA %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
STMN1_GABA$feature <- gsub("STMN1", "GABAergic", STMN1_GABA$feature)
STMN1_Glu <- STMN1_Glu %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
STMN1_Glu$feature <- gsub("STMN1", "Glutamatergic", STMN1_Glu$feature)
STMN1_Non <- STMN1_Non %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
STMN1_Non$feature <- gsub("STMN1", "Non-neuronal", STMN1_Non$feature)
# Concatenate data frames
STMN1 <- full_join(full_join(STMN1_GABA, STMN1_Glu),STMN1_Non)


# Convert from row to column 
STMN2_GABA <- STMN2_GABA %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
STMN2_GABA$feature <- gsub("STMN2", "GABAergic", STMN2_GABA$feature)
STMN2_Glu <- STMN2_Glu %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
STMN2_Glu$feature <- gsub("STMN2", "Glutamatergic", STMN2_Glu$feature)
STMN2_Non <- STMN2_Non %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
STMN2_Non$feature <- gsub("STMN2", "Non-neuronal", STMN2_Non$feature)
# Concatenate data frames
STMN2 <- full_join(full_join(STMN2_GABA, STMN2_Glu),STMN2_Non)


# Convert from row to column 
BASP1_GABA <- BASP1_GABA %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
BASP1_GABA$feature <- gsub("BASP1", "GABAergic", BASP1_GABA$feature)
BASP1_Glu <- BASP1_Glu %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
BASP1_Glu$feature <- gsub("BASP1", "Glutamatergic", BASP1_Glu$feature)
BASP1_Non <- BASP1_Non %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
BASP1_Non$feature <- gsub("BASP1", "Non-neuronal", BASP1_Non$feature)
# Concatenate data frames
BASP1 <- full_join(full_join(BASP1_GABA, BASP1_Glu),BASP1_Non)


# Convert from row to column 
ITM2B_GABA <- ITM2B_GABA %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
ITM2B_GABA$feature <- gsub("ITM2B", "GABAergic", ITM2B_GABA$feature)
ITM2B_Glu <- ITM2B_Glu %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
ITM2B_Glu$feature <- gsub("ITM2B", "Glutamatergic", ITM2B_Glu$feature)
ITM2B_Non <- ITM2B_Non %>% pivot_longer(!feature, names_to = "Cell_Type", values_to = "Values")
ITM2B_Non$feature <- gsub("ITM2B", "Non-neuronal", ITM2B_Non$feature)
# Concatenate data frames
ITM2B <- full_join(full_join(ITM2B_GABA, ITM2B_Glu),ITM2B_Non)



# 1. NNAT
NNAT_GABA_values <- NNAT_GABA[order(NNAT_GABA$Values, decreasing = TRUE), ]
NNAT_Glu_values <- NNAT_Glu[order(NNAT_Glu$Values, decreasing = TRUE), ]
NNAT_Non_values <- NNAT_Non[order(NNAT_Non$Values, decreasing = TRUE), ]
library(dplyr)
NNAT_new <- full_join(full_join(NNAT_GABA_values, NNAT_Glu_values),NNAT_Non_values)

# 2. STMN1
STMN1_GABA_values <- STMN1_GABA[order(STMN1_GABA$Values, decreasing = TRUE), ]
STMN1_Glu_values <- STMN1_Glu[order(STMN1_Glu$Values, decreasing = TRUE), ]
STMN1_Non_values <- STMN1_Non[order(STMN1_Non$Values, decreasing = TRUE), ]
STMN1_new <- full_join(full_join(STMN1_GABA_values, STMN1_Glu_values),STMN1_Non_values)

# 3. STMN2
STMN2_GABA_values <- STMN2_GABA[order(STMN2_GABA$Values, decreasing = TRUE), ]
STMN2_Glu_values <- STMN2_Glu[order(STMN2_Glu$Values, decreasing = TRUE), ]
STMN2_Non_values <- STMN2_Non[order(STMN2_Non$Values, decreasing = TRUE), ]
STMN2_new <- full_join(full_join(STMN2_GABA_values, STMN2_Glu_values),STMN2_Non_values)

# 4. BASP1
BASP1_GABA_values <- BASP1_GABA[order(BASP1_GABA$Values, decreasing = TRUE), ]
BASP1_Glu_values <- BASP1_Glu[order(BASP1_Glu$Values, decreasing = TRUE), ]
BASP1_Non_values <- BASP1_Non[order(BASP1_Non$Values, decreasing = TRUE), ]
BASP1_new <- full_join(full_join(BASP1_GABA_values, BASP1_Glu_values),BASP1_Non_values)

# 5. ITM2B
ITM2B_GABA_values <- ITM2B_GABA[order(ITM2B_GABA$Values, decreasing = TRUE), ]
ITM2B_Glu_values <- ITM2B_Glu[order(ITM2B_Glu$Values, decreasing = TRUE), ]
ITM2B_Non_values <- ITM2B_Non[order(ITM2B_Non$Values, decreasing = TRUE), ]
ITM2B_new <- full_join(full_join(ITM2B_GABA_values, ITM2B_Glu_values),ITM2B_Non_values)



  # Custom graph frame
  old_par <- par(mar = c(10, 6, 4, 0), xpd = NA)
  # Reduce font size by 50%
  par(cex.axis = 0.5)
  # Bar chart ---------------------------------------------------------------
  
# NNAT
  A <- barplot(NNAT_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n",
               names.arg = NNAT_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1)
  # Add graph name and x-axis
  mtext(paste("Neuronatin"), side = 1, 
        line = -30, font = 2, cex = 2)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 1, cex = 1)
  
  
# STMN1
  B <- barplot(STMN1_new$Values, 
              col = "black",
              ylim = c(0, 12),
              ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n",
              names.arg = STMN1_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1)
  # Add graph name and x-axis
  mtext(paste("Stathmin 1"), side = 1, 
        line = -30, font = 2, cex = 2)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 1, cex = 1)
  
  
# STMN2
  C <- barplot(STMN2_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n",
               names.arg = STMN2_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1)
  # Add graph name and x-axis
  mtext(paste("Stathmin like-2"), side = 1, 
        line = -30, font = 2, cex = 2)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 1, cex = 1)
  
  
# BASP1
  D <- barplot(BASP1_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n",
               names.arg = BASP1_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1)
  # Add graph name and x-axis
  mtext(paste("Brain abundant, membrane attached signal protein 1"), side = 1, 
        line = -30, font = 2, cex = 2)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 1, cex = 1)
  
  
# ITM2B
  E <- barplot(ITM2B_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n",
               names.arg = ITM2B_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1)
  # Add graph name and x-axis
  mtext(paste("Integral membrane protein 2B"), side = 1, 
        line = -30, font = 2, cex = 2)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 1, cex = 1)


# Mann-Whitney U-Test:
  # NNAT
  wilcox.test(NNAT_GABA$Values, NNAT_Glu$Values)
  # STMN1
  wilcox.test(STMN1_GABA$Values, STMN1_Glu$Values)
  # STMN2
  wilcox.test(STMN2_GABA$Values, STMN2_Glu$Values)
  # BASP1
  wilcox.test(BASP1_GABA$Values, BASP1_Glu$Values)
  # ITM2B
  wilcox.test(ITM2B_GABA$Values, ITM2B_Glu$Values)
  
  
  
  # Comparative bar chart  ---------------------------------------------------------------
  # NNAT
  old_par <- par(mar = c(6, 6, 4, 2), xpd = NA)
  boxplot(Values~as.factor(feature), col = "white", data = NNAT, 
          ylab = "Average log2 RPKM",
          ylim = c(0, 10), yaxt = "n",
          xlab = "Cell Type Broad Class", cex.axis = 1.25,
          main = "Neuronatin", cex.main = 2, 
          cex.lab = 1.5, col.lab = "black")
  # Adjust angle and text size on y axis
  axis(2, las = 2, cex.axis = 1.25)
  # Add means
  means <-tapply (NNAT$Values, NNAT$feature, mean)
  means
  points (means, col = "red", pch = 18)
  # Sign of statistically significant difference
  text(1, max(NNAT_GABA$Values) + 1, "a", cex=2.5, font = 2)
  text(2, max(NNAT_Glu$Values) + 1, "b", cex=2.5, font = 2)
  
  
  # STMN1
  old_par <- par(mar = c(6, 6, 4, 2), xpd = NA)
  boxplot(Values~as.factor(feature), col = "white", data = STMN1, 
          ylab = "Average log2 RPKM", 
          ylim = c(0, 12), yaxt = "n",
          xlab = "Cell Type Broad Class", cex.axis = 1.25,
          main = "Stathmin 1", cex.main = 2,
          cex.lab = 1.5, col.lab = "black")
  # Adjust angle and text size on y axis
  axis(2, las = 2, cex.axis = 1.25)
  # Add means
  means <-tapply (STMN1$Values, STMN1$feature, mean)
  means
  points (means, col = "red", pch = 18)
  # Sign of statistically significant difference
  text(1, max(STMN1_GABA$Values) + 1, "a", cex=2.5, font = 2)
  text(2, max(STMN1_Glu$Values) + 1, "b", cex=2.5, font = 2)
  
  
  # STMN2
  old_par <- par(mar = c(6, 6, 4, 2), xpd = NA)
  boxplot(Values~as.factor(feature), col = "white", data = STMN2, 
          ylab = "Average log2 RPKM", 
          ylim = c(0, 10), yaxt = "n",
          xlab = "Cell Type Broad Class", cex.axis = 1.25,
          main = "Stathmin like-2", cex.main = 2,
          cex.lab = 1.5, col.lab = "black")
  # Adjust angle and text size on y axis
  axis(2, las = 2, cex.axis = 1.25)
  # Add means
  means <-tapply (STMN2$Values, STMN2$feature, mean)
  means
  points (means, col = "red", pch = 18)
  # Sign of statistically significant difference
  text(1, max(STMN2_GABA$Values) + 1, "a", cex=2.5, font = 2)
  text(2, max(STMN2_Glu$Values) + 1, "b", cex=2.5, font = 2)
  
  
  # BASP1 
  old_par <- par(mar = c(6, 6, 4, 2), xpd = NA)
  boxplot(Values~as.factor(feature), col = "white", data = BASP1, 
          ylab = "Average log2 RPKM", 
          ylim = c(0, 10), yaxt = "n",
          xlab = "Cell Type Broad Class", cex.axis = 1.25,
          main = "Brain abundant, membrane attached signal protein 1", cex.main = 2,
          cex.lab = 1.5, col.lab = "black")
  # Adjust angle and text size on y axis
  axis(2, las = 2, cex.axis = 1.25)
  # Add means
  means <-tapply (BASP1$Values, BASP1$feature, mean)
  means
  points (means, col = "red", pch = 18)      
  # Sign of statistically significant difference
  text(1, max(BASP1_GABA$Values) + 1, "a", cex=2.5, font = 2)
  text(2, max(BASP1_Glu$Values) + 1, "b", cex=2.5, font = 2)

  
  # ITM2B 
  old_par <- par(mar = c(6, 6, 4, 2), xpd = NA)
  boxplot(Values~as.factor(feature), col = "white", data = ITM2B, 
          ylab = "Average log2 RPKM", 
          ylim = c(0, 12), yaxt = "n",
          xlab = "Cell Type Broad Class", cex.axis = 1.25,
          main = "Integral membrane protein 2B", cex.main = 2,
          cex.lab = 1.5, col.lab = "black")
  # Adjust angle and text size on y axis
  axis(2, las = 2, cex.axis = 1.25)
  # Add means
  means <-tapply (ITM2B$Values, ITM2B$feature, mean)
  means
  points (means, col = "red", pch = 18)
  # Sign of statistically significant difference
  text(1, max(ITM2B_GABA$Values) + 1, "a", cex=2.5, font = 2)
  text(2, max(ITM2B_Glu$Values) + 1, "b", cex=2.5, font = 2)
