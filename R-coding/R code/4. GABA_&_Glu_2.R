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
GABAergic2 <- read.xlsx("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/GABA2.csv")
Glutamatergic2 <- read.xlsx("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Glu2.csv")
Nonneuronal2 <- read.xlsx("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Nonneuronal2.csv")

write.csv(GABAergic, "D:/baitapR/GABA2.csv", row.names = FALSE)
write.csv(Glutamatergic, "D:/baitapR/Glu2.csv", row.names = FALSE)
write.csv(Nonneuronal, "D:/baitapR/Nonneuronal2.csv", row.names = FALSE)

GABAergic <- read.csv("D:/baitapR/GABA2.csv")
Glutamatergic <- read.csv("D:/baitapR/Glu2.csv")
Nonneuronal <- read.csv("D:/baitapR/Nonneuronal2.csv")

# Get data of 5 research genes 
## NNAT
NNAT_GABA <- GABAergic[3, ]
NNAT_Glu <- Glutamatergic[3, ]
NNAT_Non <- Nonneuronal[3, ]


## STMN1
STMN1_GABA <- GABAergic[4, ]
STMN1_Glu <- Glutamatergic[4, ]
STMN1_Non <- Nonneuronal[4, ]

## STMN2
STMN2_GABA <- GABAergic[5, ]
STMN2_Glu <- Glutamatergic[5, ]
STMN2_Non <- Nonneuronal[5, ]

## BASP1
BASP1_GABA <- GABAergic[1, ]
BASP1_Glu <- Glutamatergic[1, ]
BASP1_Non <- Nonneuronal[1, ]

## ITM2B
ITM2B_GABA <- GABAergic[2, ]
ITM2B_Glu <- Glutamatergic[2, ]
ITM2B_Non <- Nonneuronal[2, ]



# Convert from row to column 
library(tidyverse)
  # NNAT
NNAT_GABA <- NNAT_GABA %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values") 
NNAT_GABA$symbol <- gsub("Nnat", "GABAergic", NNAT_GABA$symbol)
NNAT_Glu <- NNAT_Glu %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
NNAT_Glu$symbol <- gsub("Nnat", "Glutamatergic", NNAT_Glu$symbol)
NNAT_Non <- NNAT_Non %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
NNAT_Non$symbol <- gsub("Nnat", "Non-neuronal", NNAT_Non$symbol)
# Concatenate data frames
NNAT <- full_join(full_join(NNAT_GABA, NNAT_Glu),NNAT_Non)


  # STMN1
STMN1_GABA <- STMN1_GABA %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values") 
STMN1_GABA$symbol <- gsub("Stmn1", "GABAergic", STMN1_GABA$symbol)
STMN1_Glu <- STMN1_Glu %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
STMN1_Glu$symbol <- gsub("Stmn1", "Glutamatergic", STMN1_Glu$symbol)
STMN1_Non <- STMN1_Non %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
STMN1_Non$symbol <- gsub("Stmn1", "Non-neuronal", STMN1_Non$symbol)
# Concatenate data frames
STMN1 <- full_join(full_join(STMN1_GABA, STMN1_Glu),STMN1_Non)


# STMN2
STMN2_GABA <- STMN2_GABA %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values") 
STMN2_GABA$symbol <- gsub("Stmn2", "GABAergic", STMN2_GABA$symbol)
STMN2_Glu <- STMN2_Glu %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
STMN2_Glu$symbol <- gsub("Stmn2", "Glutamatergic", STMN2_Glu$symbol)
STMN2_Non <- STMN2_Non %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
STMN2_Non$symbol <- gsub("Stmn2", "Non-neuronal", STMN2_Non$symbol)
# Concatenate data frames
STMN2 <- full_join(full_join(STMN2_GABA, STMN2_Glu),STMN2_Non)


# BASP1
BASP1_GABA <- BASP1_GABA %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values") 
BASP1_GABA$symbol <- gsub("Basp1", "GABAergic", BASP1_GABA$symbol)
BASP1_Glu <- BASP1_Glu %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
BASP1_Glu$symbol <- gsub("Basp1", "Glutamatergic", BASP1_Glu$symbol)
BASP1_Non <- BASP1_Non %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
BASP1_Non$symbol <- gsub("Basp1", "Non-neuronal", BASP1_Non$symbol)
# Concatenate data frames
BASP1 <- full_join(full_join(BASP1_GABA, BASP1_Glu),BASP1_Non)


# ITM2B
ITM2B_GABA <- ITM2B_GABA %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values") 
ITM2B_GABA$symbol <- gsub("Itm2b", "GABAergic", ITM2B_GABA$symbol)
ITM2B_Glu <- ITM2B_Glu %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
ITM2B_Glu$symbol <- gsub("Itm2b", "Glutamatergic", ITM2B_Glu$symbol)
ITM2B_Non <- ITM2B_Non %>% pivot_longer(!symbol, names_to = "Cell_Type", values_to = "Values")
ITM2B_Non$symbol <- gsub("Itm2b", "Non-neuronal", ITM2B_Non$symbol)
# Concatenate data frames
ITM2B <- full_join(full_join(ITM2B_GABA, ITM2B_Glu),ITM2B_Non)


  # Krusral - Wallis -----------------------------------------------------------
  # 1. NNAT ----
# Krusral - Wallis test
res_NNAT <- kruskal.test(Values ~ symbol, data = NNAT)
res_NNAT$p.value
# Effect size
kruskal_effsize(Values ~ symbol, data = NNAT)
# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
NNATdunn <- NNAT %>% 
  dunn_test(Values ~ symbol, p.adjust.method = "bonferroni") 
NNATdunn

NNAT_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                      "p" = c(NNATdunn[1,8], NNATdunn[2,8], NNATdunn[3,8]))
cldList(p ~ Comparison,NNAT_dunn)


# 2. STMN1 ----
# Krusral - Wallis test
res_STMN1 <- kruskal.test(Values ~ symbol, data = STMN1)
res_STMN1
# Effect size
kruskal_effsize(Values ~ symbol, data = STMN1)
# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
STMN1dunn <- STMN1 %>% 
  dunn_test(Values ~ symbol, p.adjust.method = "bonferroni") 
STMN1dunn

STMN1_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                    "p" = c(STMN1dunn[1,8], STMN1dunn[2,8], STMN1dunn[3,8]))
cldList(p ~ Comparison,STMN1_dunn)


# 3. STMN2 ----
# Krusral - Wallis test
res_STMN2 <- kruskal.test(Values ~ symbol, data = STMN2)
res_STMN2
# Effect size
kruskal_effsize(Values ~ symbol, data = STMN2)
# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
STMN2dunn <- STMN2 %>% 
  dunn_test(Values ~ symbol, p.adjust.method = "bonferroni") 
STMN2dunn

STMN2_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                    "p" = c(STMN2dunn[1,8], STMN2dunn[2,8], STMN2dunn[3,8]))
cldList(p ~ Comparison,STMN2_dunn)


# 4. BASP1 ----
# Krusral - Wallis test
res_BASP1 <- kruskal.test(Values ~ symbol, data = BASP1)
res_BASP1
# Effect size
kruskal_effsize(Values ~ symbol, data = BASP1)
# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
BASP1dunn <- BASP1 %>% 
  dunn_test(Values ~ symbol, p.adjust.method = "bonferroni") 
BASP1dunn

BASP1_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                    "p" = c(BASP1dunn[1,8], BASP1dunn[2,8], BASP1dunn[3,8]))
cldList(p ~ Comparison,BASP1_dunn)


# 5. ITM2B ----
# Krusral - Wallis test
res_ITM2B <- kruskal.test(Values ~ symbol, data = ITM2B)
res_ITM2B
# Effect size
kruskal_effsize(Values ~ symbol, data = ITM2B)
# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
ITM2Bdunn <- ITM2B %>% 
  dunn_test(Values ~ symbol, p.adjust.method = "bonferroni") 
ITM2Bdunn

ITM2Bdunn$p.adj

ITM2B_dunn <- tibble("Comparison" = c("GABA - Glu", "GABA - Non", "Glu - Non"),
                    "p" = c(ITM2Bdunn[1,8], ITM2Bdunn[2,8], ITM2Bdunn[3,8]))
cldList(p ~ Comparison,ITM2B_dunn)


  # Mann-Whitney U-Test:--------------------------------------------------------
# NNAT
wilcox.test(NNAT_GABA$Values, NNAT_Glu$Values)
wilcox.test(NNAT_GABA$Values, NNAT_Non$Values)
wilcox.test(NNAT_Glu$Values, NNAT_Non$Values)
# STMN1
wilcox.test(STMN1_GABA$Values, STMN1_Glu$Values)
wilcox.test(STMN1_GABA$Values, STMN1_Non$Values)
wilcox.test(STMN1_Glu$Values, STMN1_Non$Values)
# STMN2
wilcox.test(STMN2_GABA$Values, STMN2_Glu$Values)
wilcox.test(STMN2_GABA$Values, STMN2_Non$Values)
wilcox.test(STMN2_Glu$Values, STMN2_Non$Values)
 # BASP1
wilcox.test(BASP1_GABA$Values, BASP1_Glu$Values)
wilcox.test(BASP1_GABA$Values, BASP1_Non$Values)
wilcox.test(BASP1_Glu$Values, BASP1_Non$Values)
# ITM2B
wilcox.test(ITM2B_GABA$Values, ITM2B_Glu$Values)
wilcox.test(ITM2B_GABA$Values, ITM2B_Non$Values)
wilcox.test(ITM2B_Glu$Values, ITM2B_Non$Values)



  # Comparative boxplot  -----------------------------------------------------
#1. NNAT ------------------
ggplot(NNAT, aes(x = symbol, y = Values, shape = symbol)) +
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour=NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 2) +
  theme_pubr() +
  labs(
    title = "NNAT",
    x = "Cell Type Broad Class",
    y = "log2 RPKM",
    subtitles = ifelse(a$p.value < 0.0001, 
                       "Krusral - Wallis, p < 0.0001",
                       paste("Krusral - Wallis, p =",
                             format(round(res_NNAT$p.value, 4),
                                    scientific = FALSE)))
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 5000),
                     breaks = seq(0, 5000, 1000)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title = element_text(size = 35,
                              face = "bold",
                              hjust = 0.5),
    plot.subtitle = element_text(size = 18,
                                 face = "bold.italic",
                                 hjust = 1),
    axis.text = element_text(size = 20,
                             face = "bold"),
    axis.title = element_text(size = 23,
                              face = "bold"),
    axis.text.x = element_text(margin = margin(t = 15)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x=1.3, y=950, label="a", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=2.3, y=950, label="a", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=3.3, y=700, label="a", color="#FF0033", fontface="bold", size=14)



#2. STMN1 ------------------
ggplot(STMN1, aes(x = symbol, y = Values, shape = symbol)) +
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour=NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 2) +
  theme_pubr() +
  labs(
    title = "STMN1",
    x = "Cell Type Broad Class",
    y = "log2 RPKM",
    subtitle = ifelse(b$p.value < 0.0001, 
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =", 
                            format(round(res_STMN1$p.value, 4),
                                   scientific = FALSE)))
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 1500)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title = element_text(size = 35,
                              face = "bold",
                              hjust = 0.5),
    plot.subtitle = element_text(size = 18,
                                 face = "bold.italic",
                                 hjust = 1),
    axis.text = element_text(size = 20,
                             face = "bold"),
    axis.title = element_text(size = 23,
                              face = "bold"),
    axis.text.x = element_text(margin = margin(t = 15)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x=1.3, y=375, label="a", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=2.3, y=525, label="b", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=3.3, y=475, label="a", color="#FF0033", fontface="bold", size=14)



#3. STMN2 ------------------
ggplot(STMN2, aes(x = symbol, y = Values, shape = symbol)) +
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour=NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 2) +
  theme_pubr() +
  labs(
    title = "STMN2",
    x = "Cell Type Broad Class",
    y = "log2 RPKM",
    subtitles = ifelse(c$p.value < 0.0001, 
                       "Krusral - Wallis, p < 0.0001",
                       paste("Krusral - Wallis, p =", 
                             format(round(res_STMN2$p.value, 4),
                                    scientific = FALSE)))
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 3000)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title = element_text(size = 35,
                              face = "bold",
                              hjust = 0.5),
    plot.subtitle = element_text(size = 18,
                                 face = "bold.italic",
                                 hjust = 1),
    axis.text = element_text(size = 20,
                             face = "bold"),
    axis.title = element_text(size = 23,
                              face = "bold"),
    axis.text.x = element_text(margin = margin(t = 15)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x=1.3, y=900, label="a", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=2.3, y=1150, label="b", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=3.3, y=1375, label="a", color="#FF0033", fontface="bold", size=14)



#4. BASP1 ------------------
ggplot(BASP1, aes(x = symbol, y = Values, shape = symbol)) +
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour=NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 2) +
  theme_pubr() +
  labs(
    title = "BASP1",
    x = "Cell Type Broad Class",
    y = "log2 RPKM",
    subtitles = ifelse(d$p.value < 0.0001, 
                       "Krusral - Wallis, p < 0.0001",
                       paste("Krusral - Wallis, p =", 
                             format(round(res_BASP1$p.value, 4),
                                    scientific = FALSE)))
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 3000)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title = element_text(size = 35,
                              face = "bold",
                              hjust = 0.5),
    plot.subtitle = element_text(size = 18,
                                 face = "bold.italic",
                                 hjust = 1),
    axis.text = element_text(size = 20,
                             face = "bold"),
    axis.title = element_text(size = 23,
                              face = "bold"),
    axis.text.x = element_text(margin = margin(t = 15)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x=1.3, y=900, label="a", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=2.3, y=1150, label="b", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=3.3, y=975, label="a", color="#FF0033", fontface="bold", size=14)


#5. ITM2B ------------------
ggplot(ITM2B, aes(x = symbol, y = Values, shape = symbol)) +
  scale_shape_manual(values = c(19, 15, 17)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour=NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 2) +
  theme_pubr() +
  labs(
    title = "ITM2B",
    x = "Cell Type Broad Class",
    y = "log2 RPKM",
    subtitles = ifelse(e$p.value < 0.0001, 
                       "Krusral - Wallis, p < 0.0001",
                       paste("Krusral - Wallis, p =", 
                             format(round(res_ITM2B$p.value, 4),
                                    scientific = FALSE)))
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 5000)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title = element_text(size = 35,
                              face = "bold",
                              hjust = 0.5),
    axis.text = element_text(size = 20,
                             face = "bold"),
    axis.title = element_text(size = 23,
                              face = "bold"),
    plot.subtitle = element_text(size = 18,
                                 face = "bold.italic",
                                 hjust = 1),
    axis.text.x = element_text(margin = margin(t = 15)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(size = 1.1),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x=1.3, y=1250, label="a", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=2.35, y=1000, label="ab", color="#FF0033", fontface="bold", size=14) +
  annotate("text", x=3.3, y=1175, label="b", color="#FF0033", fontface="bold", size=14)


