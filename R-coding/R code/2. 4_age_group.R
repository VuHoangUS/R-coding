# Loading the appropriate libraries
library(data.table)
library(dplyr)
library(ggpubr)
library(rstatix)
library(rcompanion)

# Loading data
row_data <- read.csv("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Rows.csv")
columns_data <- read.csv("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Columns.csv")
expression_data <- read.csv("https://raw.githubusercontent.com/VuHoangUS/R-coding/refs/heads/main/R-coding/Data/Expression.csv", , header = FALSE)

# Remove the first column (header) of expression_data
expression_data <- expression_data[,-1]

# Convert "donor_age" column to list of rows
column_names <- columns_data$donor_age

# Add "donor_age" row to expression_data
expression_data_1 <- rbind(column_names, expression_data)

# Assign the first row as header
setnames(expression_data_1, as.character(expression_data_1[1,]))

# Remove row number 1
expression_data_1 <- expression_data_1[-1,]


# Add the "gene.name" column of row_data to expression_data_1
expression_data_1 <- cbind(row_data$gene.name, expression_data_1)

# Set row numbering starting from 1
rownames(expression_data_1) <- 1:nrow(expression_data_1)

# Get all values of 5 genes
NNAT <- expression_data_1[759, ]
STMN1 <- expression_data_1[1101, ]
STMN2 <- expression_data_1[1102, ]
BASP1 <- expression_data_1[84, ]
ITM2B <- expression_data_1[580, ]

# Create data tables -----------------------------------------------------------
  #1. NNAT ----
values_prenatal_NNAT <- sapply(NNAT[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(NNAT))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_childhood_NNAT <- sapply(NNAT[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(NNAT))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adolescent_NNAT <- sapply(NNAT[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(NNAT))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adulthood_NNAT <- sapply(NNAT[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(NNAT))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))

# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df_NNAT <- data.frame(age = vector_prenatal, value = values_prenatal_NNAT) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df_NNAT <- data.frame(age = vector_childhood, value = values_childhood_NNAT)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df_NNAT <- data.frame(age = vector_adolescent, value = values_adolescent_NNAT)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df_NNAT <- data.frame(age = vector_adulthood, value = values_adulthood_NNAT)

# Concatenate data frames
NNAT <- full_join(prenal_df_NNAT, childhood_df_NNAT)
NNAT <- full_join(NNAT, adolescent_df_NNAT)
NNAT <- full_join(NNAT, adulthood_df_NNAT)


  #2. STMN1 ----
values_prenatal_STMN1 <- sapply(STMN1[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(STMN1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_childhood_STMN1 <- sapply(STMN1[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(STMN1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adolescent_STMN1 <- sapply(STMN1[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(STMN1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adulthood_STMN1 <- sapply(STMN1[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(STMN1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))

# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df_STMN1 <- data.frame(age = vector_prenatal, value = values_prenatal_STMN1) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df_STMN1 <- data.frame(age = vector_childhood, value = values_childhood_STMN1)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df_STMN1 <- data.frame(age = vector_adolescent, value = values_adolescent_STMN1)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df_STMN1 <- data.frame(age = vector_adulthood, value = values_adulthood_STMN1)

# Concatenate data frames
STMN1 <- full_join(prenal_df_STMN1, childhood_df_STMN1)
STMN1 <- full_join(STMN1, adolescent_df_STMN1)
STMN1 <- full_join(STMN1, adulthood_df_STMN1)


  #3. STMN2 ----
values_prenatal_STMN2 <- sapply(STMN2[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(STMN2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_childhood_STMN2 <- sapply(STMN2[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(STMN2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adolescent_STMN2 <- sapply(STMN2[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(STMN2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adulthood_STMN2 <- sapply(STMN2[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(STMN2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))

# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df_STMN2 <- data.frame(age = vector_prenatal, value = values_prenatal_STMN2) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df_STMN2 <- data.frame(age = vector_childhood, value = values_childhood_STMN2)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df_STMN2 <- data.frame(age = vector_adolescent, value = values_adolescent_STMN2)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df_STMN2 <- data.frame(age = vector_adulthood, value = values_adulthood_STMN2)

# Concatenate data frames
STMN2 <- full_join(prenal_df_STMN2, childhood_df_STMN2)
STMN2 <- full_join(STMN2, adolescent_df_STMN2)
STMN2 <- full_join(STMN2, adulthood_df_STMN2)


  #4. BASP1 ----
values_prenatal_BASP1 <- sapply(BASP1[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_childhood_BASP1 <- sapply(BASP1[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adolescent_BASP1 <- sapply(BASP1[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adulthood_BASP1 <- sapply(BASP1[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))

# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df_BASP1 <- data.frame(age = vector_prenatal, value = values_prenatal_BASP1) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df_BASP1 <- data.frame(age = vector_childhood, value = values_childhood_BASP1)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df_BASP1 <- data.frame(age = vector_adolescent, value = values_adolescent_BASP1)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df_BASP1 <- data.frame(age = vector_adulthood, value = values_adulthood_BASP1)

# Concatenate data frames
BASP1 <- full_join(prenal_df_BASP1, childhood_df_BASP1)
BASP1 <- full_join(BASP1, adolescent_df_BASP1)
BASP1 <- full_join(BASP1, adulthood_df_BASP1)


  #5. ITM2B ----
values_prenatal_ITM2B <- sapply(ITM2B[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_childhood_ITM2B <- sapply(ITM2B[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adolescent_ITM2B <- sapply(ITM2B[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_adulthood_ITM2B <- sapply(ITM2B[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))

# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df_ITM2B <- data.frame(age = vector_prenatal, value = values_prenatal_ITM2B) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df_ITM2B <- data.frame(age = vector_childhood, value = values_childhood_ITM2B)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df_ITM2B <- data.frame(age = vector_adolescent, value = values_adolescent_ITM2B)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df_ITM2B <- data.frame(age = vector_adulthood, value = values_adulthood_ITM2B)

# Concatenate data frames
ITM2B <- full_join(prenal_df_ITM2B, childhood_df_ITM2B)
ITM2B <- full_join(ITM2B, adolescent_df_ITM2B)
ITM2B <- full_join(ITM2B, adulthood_df_ITM2B)



## Krusral - Wallis ------------------------------------------------------------
  #1. NNAT ----
# Krusral - Wallis test
res_NNAT <- kruskal.test(value ~ age, data = NNAT)
res_NNAT
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
NNAT$age <- factor(NNAT$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
NNATdunn <- NNAT %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
NNATdunn

NNAT_dunn <- tibble("Comparison" = c("Prenatal - Childhood", 
                                      "Prenatal - Adolescent", 
                                      "Prenatal - Adulthood",
                                      "Childhood - Adolescent",
                                      "Childhood - Adulthood",
                                      "Adolescent - Adulthood"),
                     "p" = c(NNATdunn[1,8],
                             NNATdunn[2,8],
                             NNATdunn[3,8],
                             NNATdunn[4,8],
                             NNATdunn[5,8],
                             NNATdunn[6,8]
                             ))
cldList(p ~ Comparison,NNAT_dunn)



  #2. STMN1 ----
# Krusral - Wallis test
res_STMN1 <- kruskal.test(value ~ age, data = STMN1)
res_STMN1
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
STMN1$age <- factor(STMN1$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
STMN1dunn <- STMN1 %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
STMN1dunn

STMN1_dunn <- tibble("Comparison" = c("Prenatal - Childhood", 
                                     "Prenatal - Adolescent", 
                                     "Prenatal - Adulthood",
                                     "Childhood - Adolescent",
                                     "Childhood - Adulthood",
                                     "Adolescent - Adulthood"),
                    "p" = c(STMN1dunn[1,8],
                            STMN1dunn[2,8],
                            STMN1dunn[3,8],
                            STMN1dunn[4,8],
                            STMN1dunn[5,8],
                            STMN1dunn[6,8]
                    ))
cldList(p ~ Comparison,STMN1_dunn)


  #3. STMN2 ----
# Krusral - Wallis test
res_STMN2 <- kruskal.test(value ~ age, data = STMN2)
res_STMN2
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
STMN2$age <- factor(STMN2$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
STMN2dunn <- STMN2 %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
STMN2dunn

STMN2_dunn <- tibble("Comparison" = c("Prenatal - Childhood", 
                                     "Prenatal - Adolescent", 
                                     "Prenatal - Adulthood",
                                     "Childhood - Adolescent",
                                     "Childhood - Adulthood",
                                     "Adolescent - Adulthood"),
                    "p" = c(STMN2dunn[1,8],
                            STMN2dunn[2,8],
                            STMN2dunn[3,8],
                            STMN2dunn[4,8],
                            STMN2dunn[5,8],
                            STMN2dunn[6,8]
                    ))
cldList(p ~ Comparison,STMN2_dunn)


  #4. BASP1 ----
# Krusral - Wallis test
res_BASP1 <- kruskal.test(value ~ age, data = BASP1)
res_BASP1
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
BASP1$age <- factor(BASP1$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
BASP1dunn <- BASP1 %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
BASP1dunn

BASP1_dunn <- tibble("Comparison" = c("Prenatal - Childhood", 
                                     "Prenatal - Adolescent", 
                                     "Prenatal - Adulthood",
                                     "Childhood - Adolescent",
                                     "Childhood - Adulthood",
                                     "Adolescent - Adulthood"),
                    "p" = c(BASP1dunn[1,8],
                            BASP1dunn[2,8],
                            BASP1dunn[3,8],
                            BASP1dunn[4,8],
                            BASP1dunn[5,8],
                            BASP1dunn[6,8]
                    ))
cldList(p ~ Comparison,BASP1_dunn)


  #5. ITM2B ----
# Krusral - Wallis test
res_ITM2B <- kruskal.test(value ~ age, data = ITM2B)
res_ITM2B
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
ITM2B$age <- factor(ITM2B$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
ITM2Bdunn <- ITM2B %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
ITM2Bdunn

ITM2B_dunn <- tibble("Comparison" = c("Prenatal - Childhood", 
                                     "Prenatal - Adolescent", 
                                     "Prenatal - Adulthood",
                                     "Childhood - Adolescent",
                                     "Childhood - Adulthood",
                                     "Adolescent - Adulthood"),
                    "p" = c(ITM2Bdunn[1,8],
                            ITM2Bdunn[2,8],
                            ITM2Bdunn[3,8],
                            ITM2Bdunn[4,8],
                            ITM2Bdunn[5,8],
                            ITM2Bdunn[6,8]
                    ))
cldList(p ~ Comparison,ITM2B_dunn)


# Comparative boxplot ----------------------------------------------------------
  #1. NNAT ----
ggplot(NNAT, aes(x = age, y = value, shape = age)) +
  scale_shape_manual(values = c(19, 15, 17, 25)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour = NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 3, fill = "black") +
  theme_pubr() +
  labs(
    title = "NNAT",
    x = "Age group",
    y = "log2 RPKM",
    subtitle = ifelse(res_NNAT$p.value < 0.0001,
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =",
                            format(round(res_NNAT$p.value, 4), 
                                   scientific = FALSE))) 
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(2, 14),
                     breaks = seq(2, 14, 2)) +
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
                                size = 23,
                                face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x = 1.35, y = 11.5, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 2.35, y = 8.9, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 3.35, y = 8, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 4.35, y = 8, label = "b", color = "#FF0033", fontface = "bold", size = 14)



  #2. STMN1 ----
ggplot(STMN1, aes(x = age, y = value, shape = age)) +
  scale_shape_manual(values = c(19, 15, 17, 25)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour = NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 3, fill = "black") +
  theme_pubr() +
  labs(
    title = "STMN1",
    x = "Age group",
    y = "log2 RPKM",
    subtitle = ifelse(res_STMN1$p.value < 0.0001,
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =",
                            format(round(res_STMN1$p.value, 4), 
                                   scientific = FALSE))) 
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(2, 10),
                     breaks = seq(2, 10, 2)) +
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
                                size = 23,
                                face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x = 1.35, y = 9.5, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 2.35, y = 8.25, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 3.35, y = 7.9, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 4.35, y = 7.7, label = "b", color = "#FF0033", fontface = "bold", size = 14)



  #3. STMN2 ----
ggplot(STMN2, aes(x = age, y = value, shape = age)) +
  scale_shape_manual(values = c(19, 15, 17, 25)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour = NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 3, fill = "black") +
  theme_pubr() +
  labs(
    title = "STMN2",
    x = "Age group",
    y = "log2 RPKM",
    subtitle = ifelse(res_STMN2$p.value < 0.0001,
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =",
                            format(round(res_STMN2$p.value, 4), 
                                   scientific = FALSE))) 
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(2, 12),
                     breaks = seq(2, 12, 2)) +
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
                                size = 23,
                                face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x = 1.35, y = 11.3, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 2.35, y = 10, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 3.35, y = 9.5, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 4.35, y = 9.25, label = "b", color = "#FF0033", fontface = "bold", size = 14)


    #4. BASP1 ----
ggplot(BASP1, aes(x = age, y = value, shape = age)) +
  scale_shape_manual(values = c(19, 15, 17, 25)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour = NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 3, fill = "black") +
  theme_pubr() +
  labs(
    title = "BASP1",
    x = "Age group",
    y = "log2 RPKM",
    subtitle = ifelse(res_BASP1$p.value < 0.0001,
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =",
                            format(round(res_BASP1$p.value, 4), 
                                   scientific = FALSE))) 
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(2, 10),
                     breaks = seq(2, 10, 2)) +
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
                                size = 23,
                                face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x = 1.35, y = 9.35, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 2.35, y = 8.7, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 3.35, y = 8.35, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 4.35, y = 8.35, label = "b", color = "#FF0033", fontface = "bold", size = 14)


  #5. ITM2B ----
ggplot(ITM2B, aes(x = age, y = value, shape = age)) +
  scale_shape_manual(values = c(19, 15, 17, 25)) +
  stat_boxplot(geom = "errorbar",
               width = 0.3,
               size = 1.2) +
  geom_boxplot(outlier.colour = NA,
               linewidth = 1.2) +
  geom_jitter(position = position_jitter(0.2),
              size = 3, fill = "black") +
  theme_pubr() +
  labs(
    title = "ITM2B",
    x = "Age group",
    y = "log2 RPKM",
    subtitle = ifelse(res_ITM2B$p.value < 0.0001,
                      "Krusral - Wallis, p < 0.0001",
                      paste("Krusral - Wallis, p =",
                            format(round(res_ITM2B$p.value, 4), 
                                   scientific = FALSE))) 
  ) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(4, 10),
                     breaks = seq(4, 10, 2)) +
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
                                size = 23,
                                face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2)
  ) + 
  annotate("text", x = 1.35, y = 6.9, label = "a", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 2.35, y = 8.5, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 3.35, y = 8.35, label = "b", color = "#FF0033", fontface = "bold", size = 14) +
  annotate("text", x = 4.35, y = 8.35, label = "b", color = "#FF0033", fontface = "bold", size = 14)
