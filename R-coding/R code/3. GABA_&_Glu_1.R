# Loading the appropriate libraries
library(datasets)
library(ggplot2)
library(multcompView)
library(dplyr)
library(openxlsx)
library(ggpubr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(rcompanion)

# Loading data
GABAergic1 <- read.csv("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/GABA1.csv")
Glutamatergic1 <- read.csv("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Glu1.csv")
Nonneuronal1 <- read.csv("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Nonneuronal1.csv")

# Get data of 5 research genes
## NNAT
NNAT_GABA <- GABAergic1[34884, ]
NNAT_Glu <- Glutamatergic1[34884, ]
NNAT_Non <- Nonneuronal1[34884, ]

## STMN1
STMN1_GABA <- GABAergic1[45460, ]
STMN1_Glu <- Glutamatergic1[45460, ]
STMN1_Non <- Nonneuronal1[45460, ]

## STMN2
STMN2_GABA <- GABAergic1[45463, ]
STMN2_Glu <- Glutamatergic1[45463, ]
STMN2_Non <- Nonneuronal1[45463, ]

## BASP1
BASP1_GABA <- GABAergic1[1934, ]
BASP1_Glu <- Glutamatergic1[1934, ]
BASP1_Non <- Nonneuronal1[1934, ]

## ITM2B
ITM2B_GABA <- GABAergic1[11758, ]
ITM2B_Glu <- Glutamatergic1[11758, ]
ITM2B_Non <- Nonneuronal1[11758, ]


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
  old_par <- par(mar = c(10, 6, 6, 0), xpd = NA)
  # Reduce font size by 50%
  par(cex.axis = 0.5)
  # Bar chart ---------------------------------------------------------------
  
#1. NNAT ----
  A <- barplot(NNAT_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n", cex.lab = 1.5, font.lab = 2, font = 2,
               names.arg = NNAT_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1.5, lwd = 2, font = 2)
  # Add graph name and x-axis
  mtext(paste("NNAT"), side = 1, 
        line = -29, font = 2, cex = 3.5)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 2, cex = 1.5)
  
  
#2. STMN1 ----
  B <- barplot(STMN1_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n", cex.lab = 1.5, font.lab = 2, font = 2,
               names.arg = NNAT_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1.5, lwd = 2, font = 2)
  # Add graph name and x-axis
  mtext(paste("STMN1"), side = 1, 
        line = -29, font = 2, cex = 3.5)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 2, cex = 1.5)
  
  
#3. STMN2 ----
  C <- barplot(STMN2_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n", cex.lab = 1.5, font.lab = 2, font = 2,
               names.arg = NNAT_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1.5, lwd = 2, font = 2)
  # Add graph name and x-axis
  mtext(paste("STMN2"), side = 1, 
        line = -29, font = 2, cex = 3.5)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 2, cex = 1.5)
  
  
#4. BASP1 ----
  D <- barplot(BASP1_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n", cex.lab = 1.5, font.lab = 2, font = 2,
               names.arg = NNAT_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1.5, lwd = 2, font = 2)
  # Add graph name and x-axis
  mtext(paste("BASP1"), side = 1, 
        line = -29, font = 2, cex = 3.5)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 2, cex = 1.5)
  
   
#5. ITM2B ----
  E <- barplot(ITM2B_new$Values, 
               col = "black",
               ylim = c(0, 12),
               ylab = "Trimmed Mean Expression (log2 RPKM)", main ="", yaxt = "n", cex.lab = 1.5, font.lab = 2, font = 2,
               names.arg = NNAT_new$Cell_Type, las=2)
  # Add y-axis
  axis(2, las = 2, cex.axis = 1.5, lwd = 2, font = 2)
  # Add graph name and x-axis
  mtext(paste("ITM2B"), side = 1, 
        line = -29, font = 2, cex = 3.5)
  mtext(paste("Cell Type Label"), side = 1, 
        line = 8, font = 2, cex = 1.5)
  
  
  # Krusral - Wallis -----------------------------------------------------------
#1. NNAT ----
  # Krusral - Wallis test
  res_NNAT <- kruskal.test(Values ~ feature, data = NNAT)
  res_NNAT
  # Effect size
  kruskal_effsize(Values ~ feature, data = NNAT)
  # Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
  NNATdunn <- NNAT %>% 
    dunn_test(Values ~ feature, p.adjust.method = "bonferroni") 
  NNATdunn
  
  NNAT_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                      "p" = c(NNATdunn[1,8], NNATdunn[2,8], NNATdunn[3,8]))
  cldList(p ~ Comparison,NNAT_dunn)
  
  
#2. STMN1 ----
  # Krusral - Wallis test
  res_STMN1 <- kruskal.test(Values ~ feature, data = STMN1)
  res_STMN1
  # Effect size
  kruskal_effsize(Values ~ feature, data = STMN1)
  # Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
  STMN1dunn <- STMN1 %>% 
    dunn_test(Values ~ feature, p.adjust.method = "bonferroni") 
  STMN1dunn
  
  STMN1_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                      "p" = c(STMN1dunn[1,8], STMN1dunn[2,8], STMN1dunn[3,8]))
  cldList(p ~ Comparison,STMN1_dunn)
  
  
#3. STMN2 ----
  # Krusral - Wallis test
  res_STMN2 <- kruskal.test(Values ~ feature, data = STMN2)
  res_STMN2
  # Effect size
  kruskal_effsize(Values ~ feature, data = STMN2)
  # Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
  STMN2dunn <- STMN2 %>% 
    dunn_test(Values ~ feature, p.adjust.method = "bonferroni") 
  STMN2dunn
  
  STMN2_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                      "p" = c(STMN2dunn[1,8], STMN2dunn[2,8], STMN2dunn[3,8]))
  cldList(p ~ Comparison,STMN2_dunn)

  
#4. BASP1 ----
  # Krusral - Wallis test
  res_BASP1 <- kruskal.test(Values ~ feature, data = BASP1)
  res_BASP1
  # Effect size
  kruskal_effsize(Values ~ feature, data = BASP1)
  # Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
  BASP1dunn <- BASP1 %>% 
    dunn_test(Values ~ feature, p.adjust.method = "bonferroni") 
  BASP1dunn
  
  BASP1_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                      "p" = c(BASP1dunn[1,8], BASP1dunn[2,8], BASP1dunn[3,8]))
  cldList(p ~ Comparison,BASP1_dunn)
  
  
  #5. ITM2B ----
  # Krusral - Wallis test
  res_ITM2B <- kruskal.test(Values ~ feature, data = ITM2B)
  res_ITM2B
  # Effect size
  kruskal_effsize(Values ~ feature, data = ITM2B)
  # Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
  ITM2Bdunn <- ITM2B %>% 
    dunn_test(Values ~ feature, p.adjust.method = "bonferroni") 
  ITM2Bdunn
  
  ITM2B_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                      "p" = c(ITM2Bdunn[1,8], ITM2Bdunn[2,8], ITM2Bdunn[3,8]))
  cldList(p ~ Comparison,ITM2B_dunn)

  
# Mann-Whitney U-Test ----------------------------------------------------------
  #1. NNAT
  wilcox.test(NNAT_GABA$Values, NNAT_Glu$Values)
  wilcox.test(NNAT_GABA$Values, NNAT_Non$Values)
  wilcox.test(NNAT_Glu$Values, NNAT_Non$Values)
  #2. STMN1
  wilcox.test(STMN1_GABA$Values, STMN1_Glu$Values)
  wilcox.test(STMN1_GABA$Values, STMN1_Non$Values)
  wilcox.test(STMN1_Glu$Values, STMN1_Non$Values)
  #3. STMN2
  wilcox.test(STMN2_GABA$Values, STMN2_Glu$Values)
  wilcox.test(STMN2_GABA$Values, STMN2_Non$Values)
  wilcox.test(STMN2_Glu$Values, STMN2_Non$Values)
  #4. BASP1
  wilcox.test(BASP1_GABA$Values, BASP1_Glu$Values)
  wilcox.test(BASP1_GABA$Values, BASP1_Non$Values)
  wilcox.test(BASP1_Glu$Values, BASP1_Non$Values)
  #5. ITM2B
  wilcox.test(ITM2B_GABA$Values, ITM2B_Glu$Values)
  wilcox.test(ITM2B_GABA$Values, ITM2B_Non$Values)
  wilcox.test(ITM2B_Glu$Values, ITM2B_Non$Values)
  
  
  # Comparative boxplot  ---------------------------------------------------------------
  #1. NNAT ----
ggplot(NNAT, aes(x = feature, y = Values, shape = feature)) +
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour = NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 3) +
  theme_pubr() +
  labs(
    title = "NNAT",
    x = "Cell Type Broad Class",
    y = "Trimmed Mean Expression (log2 RPKM)",
    subtitle = ifelse(res_NNAT$p.value < 0.0001,
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =",
                            format(round(res_NNAT$p.value, 4), 
                                   scientific = FALSE))) 
    ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 10),
                     breaks = seq(0, 10, 2)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 35,
                              face = "bold",
                              hjust = 0.5),
    plot.subtitle = element_text(size = 18,
                                 face = "bold.italic",
                                 hjust = 1),
    axis.text = element_text(size = 20,
                             face =  "bold"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15),
                                size = 23,
                                face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10),
                                size = 18,
                                face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2)
  ) + 
    annotate("text", x = 1.3, y = 0.75, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 2.3, y = 4.25, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 3.3, y = 0.75, label = "a", color = "#FF0033", fontface = "bold", size = 14)

  
  # 2. STMN1 ----
  ggplot(STMN1, aes(x = feature, y = Values, shape = feature)) +
    scale_shape_manual(values = c(19, 15, 17)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.3,
                 size = 1.2) +
    geom_boxplot(outlier.colour = NA,
                 linewidth = 1.2) +
    geom_jitter(position = position_jitter(0.2),
                size = 3) +
    theme_pubr() +
    labs(
      title = "STMN1",
      x = "Cell Type Broad Class",
      y = "Trimmed Mean Expression (log2 RPKM)",
      subtitle = ifelse(res_STMN1$p.value < 0.0001,
                        "Krusral - Wallis, p < 0.0001",
                        paste("Krusral - Wallis, p =",
                              format(round(res_STMN1$p.value, 4), 
                                     scientific = FALSE))) 
    ) +
    scale_y_continuous(expand = expansion(0),
                       limits = c(0, 12),
                       breaks = seq(0, 12, 2)) +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(size = 35,
                                face = "bold",
                                hjust = 0.5),
      plot.subtitle = element_text(size = 18,
                                   face = "bold.italic",
                                   hjust = 1),
      axis.text = element_text(size = 20,
                               face =  "bold"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title.x = element_text(margin = margin(t = 15),
                                  size = 23,
                                  face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10),
                                  size = 18,
                                  face = "bold"),
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2)
    ) + 
    annotate("text", x = 1.3, y = 9.25, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 2.3, y = 10, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 3.3, y = 1, label = "c", color = "#FF0033", fontface = "bold", size = 14)
  
  
  # 3. STMN2 ----
  ggplot(STMN2, aes(x = feature, y = Values, shape = feature)) +
    scale_shape_manual(values = c(19, 15, 17)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.3,
                 size = 1.2) +
    geom_boxplot(outlier.colour = NA,
                 linewidth = 1.2) +
    geom_jitter(position = position_jitter(0.2),
                size = 3) +
    theme_pubr() +
    labs(
      title = "STMN2",
      x = "Cell Type Broad Class",
      y = "Trimmed Mean Expression (log2 RPKM)",
      subtitle = ifelse(res_STMN2$p.value < 0.0001,
                        "Krusral - Wallis, p < 0.0001",
                        paste("Krusral - Wallis, p =",
                              format(round(res_STMN2$p.value, 4), 
                                     scientific = FALSE))) 
    ) +
    scale_y_continuous(expand = expansion(0),
                       limits = c(0, 10),
                       breaks = seq(0, 10, 2)) +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(size = 35,
                                face = "bold",
                                hjust = 0.5),
      plot.subtitle = element_text(size = 18,
                                   face = "bold.italic",
                                   hjust = 1),
      axis.text = element_text(size = 20,
                               face =  "bold"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title.x = element_text(margin = margin(t = 15),
                                  size = 23,
                                  face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10),
                                  size = 18,
                                  face = "bold"),
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2)
    ) + 
    annotate("text", x = 1.3, y = 8, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 2.3, y = 8.85, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 3.3, y = 0.75, label = "c", color = "#FF0033", fontface = "bold", size = 14)
  
  
  # 4. BASP1 ----
  ggplot(BASP1, aes(x = feature, y = Values, shape = feature)) +
    scale_shape_manual(values = c(19, 15, 17)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.3,
                 size = 1.2) +
    geom_boxplot(outlier.colour = NA,
                 linewidth = 1.2) +
    geom_jitter(position = position_jitter(0.2),
                size = 3) +
    theme_pubr() +
    labs(
      title = "BASP1",
      x = "Cell Type Broad Class",
      y = "Trimmed Mean Expression (log2 RPKM)",
      subtitle = ifelse(res_BASP1$p.value < 0.0001,
                        "Krusral - Wallis, p < 0.0001",
                        paste("Krusral - Wallis, p =",
                              format(round(res_BASP1$p.value, 4), 
                                     scientific = FALSE))) 
    ) +
    scale_y_continuous(expand = expansion(0),
                       limits = c(0, 10),
                       breaks = seq(0, 10, 2)) +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(size = 35,
                                face = "bold",
                                hjust = 0.5),
      plot.subtitle = element_text(size = 18,
                                   face = "bold.italic",
                                   hjust = 1),
      axis.text = element_text(size = 20,
                               face =  "bold"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title.x = element_text(margin = margin(t = 15),
                                  size = 23,
                                  face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10),
                                  size = 18,
                                  face = "bold"),
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2)
    ) + 
    annotate("text", x = 1.3, y = 6.8, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 2.3, y = 8, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 3.3, y = 0.75, label = "c", color = "#FF0033", fontface = "bold", size = 14)

  
  # 5. ITM2B ---- 
  ggplot(ITM2B, aes(x = feature, y = Values, shape = feature)) +
    scale_shape_manual(values = c(19, 15, 17)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.3,
                 size = 1.2) +
    geom_boxplot(outlier.colour = NA,
                 linewidth = 1.2) +
    geom_jitter(position = position_jitter(0.2),
                size = 3) +
    theme_pubr() +
    labs(
      title = "ITM2B",
      x = "Cell Type Broad Class",
      y = "Trimmed Mean Expression (log2 RPKM)",
      subtitle = ifelse(res_ITM2B$p.value < 0.0001,
                        "Krusral - Wallis, p < 0.0001",
                        paste("Krusral - Wallis, p =",
                              format(round(res_ITM2B$p.value, 4), 
                                     scientific = FALSE))) 
    ) +
    scale_y_continuous(expand = expansion(0),
                       limits = c(0, 10),
                       breaks = seq(0, 10, 2)) +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(size = 35,
                                face = "bold",
                                hjust = 0.5),
      plot.subtitle = element_text(size = 18,
                                   face = "bold.italic",
                                   hjust = 1),
      axis.text = element_text(size = 20,
                               face =  "bold"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title.x = element_text(margin = margin(t = 15),
                                  size = 23,
                                  face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10),
                                  size = 18,
                                  face = "bold"),
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2)
    ) + 
    annotate("text", x = 1.3, y = 9.6, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 2.3, y = 8.6, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
    annotate("text", x = 3.3, y = 8.25, label = "b", color = "#FF0033", fontface = "bold", size = 14)
