setwd("D:/baitapR")
row_data <- read.csv("rows.csv", header = TRUE, stringsAsFactors = FALSE)
columns_data <- read.csv("columns.csv", header = TRUE, stringsAsFactors = FALSE)
expression_data <- read.csv("expression.csv", header = FALSE, stringsAsFactors = FALSE)

# Remove the first column (header) of expression_data
expression_data <- expression_data[,-1]

# Convert "donor_age" column to list of rows
column_names <- columns_data$donor_age

# Add "donor_age" row to expression_data
expression_data_1 <- rbind(column_names, expression_data)

# Assign the first row as header
library(data.table)
setnames(expression_data_1, as.character(expression_data_1[1,]))

# Remove row number 1
expression_data_1 <- expression_data_1[-1,]


# Add the "gene.name" column of row_data to expression_data_1
expression_data_1 <- cbind(row_data$gene.name, expression_data_1)

# Set row numbering starting from 1
rownames(expression_data_1) <- 1:nrow(expression_data_1)

# Get all values of genes "neuronatin", "stathmin 1", "stathmin-like 2" and "brain abundant, membrane attached signal protein 1"  
neuronatin <- expression_data_1[759, ]
stathmin_1 <- expression_data_1[1101, ]
stathmin.like_2 <- expression_data_1[1102, ]
BASP1 <- expression_data_1[84, ]
ITM2B <- expression_data_1[580, ]


### Neuronatin
# Mean of prenatal age
values_prenatal_neuronatin <- sapply(neuronatin[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(neuronatin))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_prenatal_neuronatin <- mean(values_prenatal_neuronatin, na.rm = TRUE)
mean_prenatal_neuronatin <- round(mean_prenatal_neuronatin, 4)
mean_prenatal_neuronatin
# Standard deviation of prenatal age
sd_prenatal_neuronatin <- sd(values_prenatal_neuronatin, na.rm = TRUE)

# Mean of childhood age
values_childhood_neuronatin <- sapply(neuronatin[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(neuronatin))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_childhood_neuronatin <- mean(values_childhood_neuronatin, na.rm = TRUE)
mean_childhood_neuronatin <- round(mean_childhood_neuronatin, 4)
mean_childhood_neuronatin
# Standard deviation of childhood age
sd_childhood_neuronatin <- sd(values_childhood_neuronatin, na.rm = TRUE)

# Mean of adolescent age
values_adolescent_neuronatin <- sapply(neuronatin[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(neuronatin))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adolescent_neuronatin <- mean(values_adolescent_neuronatin, na.rm = TRUE)
mean_adolescent_neuronatin <- round(mean_adolescent_neuronatin, 4)
mean_adolescent_neuronatin
# Standard deviation of adolescent age
sd_adolescent_neuronatin <- sd(values_adolescent_neuronatin, na.rm = TRUE)

# Mean of adulthood age
values_adulthood_neuronatin <- sapply(neuronatin[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(neuronatin))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adulthood_neuronatin <- mean(values_adulthood_neuronatin, na.rm = TRUE)
mean_adulthood_neuronatin <- round(mean_adulthood_neuronatin, 4)
mean_adulthood_neuronatin
# Standard deviation of adulthood age
sd_adulthood_neuronatin <- sd(values_adulthood_neuronatin, na.rm = TRUE)


# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df <- data.frame(age = vector_prenatal, value = values_prenatal_neuronatin) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df <- data.frame(age = vector_childhood, value = values_childhood_neuronatin)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df <- data.frame(age = vector_adolescent, value = values_adolescent_neuronatin)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df <- data.frame(age = vector_adulthood, value = values_adulthood_neuronatin)

# Concatenate data frames
library(dplyr)
merge_df <- full_join(prenal_df, childhood_df)
merge_df <- full_join(merge_df, adolescent_df)
merge_df <- full_join(merge_df, adulthood_df)

## Krusral - Wallis
library(ggpubr)
library(rstatix)
# Krusral - Wallis test
res.kruskal <- kruskal.test(value ~ age, data = merge_df)
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
merge_df$age <- factor(merge_df$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
pwc <- merge_df %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
pwc

# Bar chart ---------------------------------------------------------------
pwc <- pwc %>% add_xy_position(x = "age")
ggboxplot(merge_df, x = "age", y = "value") +
  scale_y_continuous(limits = c(2, 14),                   # Limit y-axis
                     breaks = seq(2, 14, 2)) +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +               # Sign of statistically significant difference
  labs(
    title = "Neuronatin",
    x = "Age group", y = "log2 RPKM") +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 25, 
                              hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(t = 15)),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20,
                              color = "black",
                              face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10))
  )



#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

### Stathmin 1
# Mean of prenatal age 
values_prenatal_stathmin_1 <- sapply(stathmin_1[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(stathmin_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_prenatal_stathmin_1 <- mean(values_prenatal_stathmin_1, na.rm = TRUE)
mean_prenatal_stathmin_1 <- round(mean_prenatal_stathmin_1, 4)
mean_prenatal_stathmin_1

sd_prenatal_stathmin_1 <- sd(values_prenatal_stathmin_1, na.rm = TRUE)

# Mean of childhood age
values_childhood_stathmin_1 <- sapply(stathmin_1[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(stathmin_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_childhood_stathmin_1 <- mean(values_childhood_stathmin_1, na.rm = TRUE)
mean_childhood_stathmin_1 <- round(mean_childhood_stathmin_1, 4)
mean_childhood_stathmin_1

sd_childhood_stathmin_1 <- sd(values_childhood_stathmin_1, na.rm = TRUE)

# Mean of adolescent age
values_adolescent_stathmin_1 <- sapply(stathmin_1[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(stathmin_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adolescent_stathmin_1 <- mean(values_adolescent_stathmin_1, na.rm = TRUE)
mean_adolescent_stathmin_1 <- round(mean_adolescent_stathmin_1, 4)
mean_adolescent_stathmin_1

sd_adolescent_stathmin_1 <- sd(values_adolescent_stathmin_1, na.rm = TRUE)

# Mean of adulthood age
values_adulthood_stathmin_1 <- sapply(stathmin_1[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(stathmin_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adulthood_stathmin_1 <- mean(values_adulthood_stathmin_1, na.rm = TRUE)
mean_adulthood_stathmin_1 <- round(mean_adulthood_stathmin_1, 4)
mean_adulthood_stathmin_1

sd_adulthood_stathmin_1 <- sd(values_adulthood_stathmin_1, na.rm = TRUE)


values_stathmin_1 <- c(mean_prenatal_stathmin_1, mean_childhood_stathmin_1, mean_adolescent_stathmin_1, mean_adulthood_stathmin_1)
sd_stathmin_1 <- c(sd_prenatal_stathmin_1, sd_childhood_stathmin_1, sd_adolescent_stathmin_1, sd_adulthood_stathmin_1)



value <- c(values_prenatal_stathmin_1, values_childhood_stathmin_1, values_adolescent_stathmin_1, values_adulthood_stathmin_1)
  
# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df <- data.frame(age = vector_prenatal, value = values_prenatal_stathmin_1) 
  
# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df <- data.frame(age = vector_childhood, value = values_childhood_stathmin_1)
  
# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df <- data.frame(age = vector_adolescent, value = values_adolescent_stathmin_1)
  
# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df <- data.frame(age = vector_adulthood, value = values_adulthood_stathmin_1)
  
# Concatenate data frames
library(dplyr)
merge_df <- full_join(prenal_df, childhood_df)
merge_df <- full_join(merge_df, adolescent_df)
merge_df <- full_join(merge_df, adulthood_df)

### Krusral - Wallis
library(ggpubr)
library(rstatix)
# Krusral - Wallis test
res.kruskal <- kruskal.test(value ~ age, data = merge_df)
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
merge_df$age <- factor(merge_df$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
pwc <- merge_df %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
pwc

# Bar chart ---------------------------------------------------------------
pwc <- pwc %>% add_xy_position(x = "age")
ggboxplot(merge_df, x = "age", y = "value") +
  scale_y_continuous(limits = c(3, 11),                    
                     breaks = seq(3, 11, 2)) +
  stat_pvalue_manual(pwc, hide.ns = TRUE) + 
  labs(
    title = "Stathmin 1",
    x = "Age group", y = "log2 RPKM") +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 25, 
                              hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(t = 15)),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20,
                              color = "black",
                              face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10))
  )


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

### Stathmin like_2
# Mean of prenatal age
values_prenatal_stathmin.like_2 <- sapply(stathmin.like_2[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(stathmin.like_2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_prenatal_stathmin.like_2 <- mean(values_prenatal_stathmin.like_2, na.rm = TRUE)
mean_prenatal_stathmin.like_2 <- round(mean_prenatal_stathmin.like_2, 4)
mean_prenatal_stathmin.like_2
# Standard deviation of prenatal age
sd_prenatal_stathmin.like_2 <- sd(values_prenatal_stathmin.like_2, na.rm = TRUE)

# Mean of childhood age 
values_childhood_stathmin.like_2 <- sapply(stathmin.like_2[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(stathmin.like_2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_childhood_stathmin.like_2 <- mean(values_childhood_stathmin.like_2, na.rm = TRUE)
mean_childhood_stathmin.like_2 <- round(mean_childhood_stathmin.like_2, 4)
mean_childhood_stathmin.like_2
# Standard deviation of childhood age
sd_childhood_stathmin.like_2 <- sd(values_childhood_stathmin.like_2, na.rm = TRUE)

# Mean of adolescent age
values_adolescent_stathmin.like_2 <- sapply(stathmin.like_2[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(stathmin.like_2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adolescent_stathmin.like_2 <- mean(values_adolescent_stathmin.like_2, na.rm = TRUE)
mean_adolescent_stathmin.like_2 <- round(mean_adolescent_stathmin.like_2, 4)
mean_adolescent_stathmin.like_2
# Standard deviation of adolescent age

sd_adolescent_stathmin.like_2 <- sd(values_adolescent_stathmin.like_2, na.rm = TRUE)

# Mean of adulthood age
values_adulthood_stathmin.like_2 <- sapply(stathmin.like_2[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(stathmin.like_2))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adulthood_stathmin.like_2 <- mean(values_adulthood_stathmin.like_2, na.rm = TRUE)
mean_adulthood_stathmin.like_2 <- round(mean_adulthood_stathmin.like_2, 4)
mean_adulthood_stathmin.like_2
# Standard deviation of adulthood age
sd_adulthood_stathmin.like_2 <- sd(values_adulthood_stathmin.like_2, na.rm = TRUE)


# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df <- data.frame(age = vector_prenatal, value = values_prenatal_stathmin.like_2) 
  
# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df <- data.frame(age = vector_childhood, value = values_childhood_stathmin.like_2)
  
# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df <- data.frame(age = vector_adolescent, value = values_adolescent_stathmin.like_2)
  
# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df <- data.frame(age = vector_adulthood, value = values_adulthood_stathmin.like_2)
  
# Concatenate data frames
library(dplyr)
merge_df <- full_join(prenal_df, childhood_df)
merge_df <- full_join(merge_df, adolescent_df)
merge_df <- full_join(merge_df, adulthood_df)


### Krusral - Wallis
library(ggpubr)
library(rstatix)
# Krusral - Wallis test
res.kruskal <- kruskal.test(value ~ age, data = merge_df)
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
merge_df$age <- factor(merge_df$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
pwc <- merge_df %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
pwc

# Bar chart ---------------------------------------------------------------
pwc <- pwc %>% add_xy_position(x = "age")
ggboxplot(merge_df, x = "age", y = "value") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) + 
  labs(
    title = "Stathmin like-2",
    x = "Age group", y = "log2 RPKM") +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 25, 
                              hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(t = 15)),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20,
                              color = "black",
                              face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10))
  )

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

### Brain abundant, membrane attached signal protein 1
# Mean of prenatal age
values_prenatal_BASP1 <- sapply(BASP1[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_prenatal_BASP1 <- mean(values_prenatal_BASP1, na.rm = TRUE)
mean_prenatal_BASP1 <- round(mean_prenatal_BASP1, 4)
mean_prenatal_BASP1
# Standard deviation of prenatal age
sd_prenatal_BASP1 <- sd(values_prenatal_BASP1, na.rm = TRUE)


# Mean of childhood age
values_childhood_BASP1 <- sapply(BASP1[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_childhood_BASP1 <- mean(values_childhood_BASP1, na.rm = TRUE)
mean_childhood_BASP1 <- round(mean_childhood_BASP1, 4)
mean_childhood_BASP1
# Standard deviation of childhood age
sd_childhood_BASP1 <- sd(values_childhood_BASP1, na.rm = TRUE)


# Mean of adolescent age
values_adolescent_BASP1 <- sapply(BASP1[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adolescent_BASP1 <- mean(values_adolescent_BASP1, na.rm = TRUE)
mean_adolescent_BASP1 <- round(mean_adolescent_BASP1, 4)
mean_adolescent_BASP1
# Standard deviation of adolescent age
sd_adolescent_BASP1 <- sd(values_adolescent_BASP1, na.rm = TRUE)


# Mean of adulthood age
values_adulthood_BASP1 <- sapply(BASP1[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(BASP1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adulthood_BASP1 <- mean(values_adulthood_BASP1, na.rm = TRUE)
mean_adulthood_BASP1 <- round(mean_adulthood_BASP1, 4)
mean_adulthood_BASP1
# Standard deviation of adulthood age
sd_adulthood_BASP1 <- sd(values_adulthood_BASP1, na.rm = TRUE)


# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df <- data.frame(age = vector_prenatal, value = values_prenatal_BASP1) 
  
# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df <- data.frame(age = vector_childhood, value = values_childhood_BASP1)
  
# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df <- data.frame(age = vector_adolescent, value = values_adolescent_BASP1)
  
# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df <- data.frame(age = vector_adulthood, value = values_adulthood_BASP1)
  
# Concatenate data frames
library(dplyr)
merge_df <- full_join(prenal_df, childhood_df)
merge_df <- full_join(merge_df, adolescent_df)
merge_df <- full_join(merge_df, adulthood_df)

### Krusral - Wallis
library(ggpubr)
library(rstatix)
# Krusral - Wallis test
res.kruskal <- kruskal.test(value ~ age, data = merge_df)
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
merge_df$age <- factor(merge_df$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
pwc <- merge_df %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
pwc

# Bar chart ---------------------------------------------------------------
pwc <- pwc %>% add_xy_position(x = "age")
ggboxplot(merge_df, x = "age", y = "value") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) + 
  labs(
    title = "Brain abundant, membrane attached signal protein 1",
    x = "Age group", y = "log2 RPKM") +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 25, 
                              hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(t = 15)),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20,
                              color = "black",
                              face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10))
    )


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

### Integral membrane protein 2B 
# Mean of prenatal age 
values_prenatal_ITM2B <- sapply(ITM2B[, grepl("^(8.pcw|9.pcw|12.pcw|13.pcw|16.pcw|17.pcw|19.pcw|21.pcw|24.pcw|25.pcw|26.pcw|35.pcw|37.pcw)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_prenatal_ITM2B <- mean(values_prenatal_ITM2B, na.rm = TRUE)
mean_prenatal_ITM2B <- round(mean_prenatal_ITM2B, 4)
mean_prenatal_ITM2B
# Standard deviation of prenatal age
sd_prenatal_ITM2B <- sd(values_prenatal_ITM2B, na.rm = TRUE)

# Mean of childhood age
values_childhood_ITM2B <- sapply(ITM2B[, grepl("^(4.mos|10.mos|1.yrs|2.yrs|3.yrs)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_childhood_ITM2B <- mean(values_childhood_ITM2B, na.rm = TRUE)
mean_childhood_ITM2B <- round(mean_childhood_ITM2B, 4)
mean_childhood_ITM2B
# Standard deviation of childhood age
sd_childhood_ITM2B <- sd(values_childhood_ITM2B, na.rm = TRUE)

# Mean of adolescent age
values_adolescent_ITM2B <- sapply(ITM2B[, grepl("^(11.yrs|13.yrs|15.yrs|18.yrs|19.yrs)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adolescent_ITM2B <- mean(values_adolescent_ITM2B, na.rm = TRUE)
mean_adolescent_ITM2B <- round(mean_adolescent_ITM2B, 4)
mean_adolescent_ITM2B
# Standard deviation of adolescent age
sd_adolescent_ITM2B <- sd(values_adolescent_ITM2B, na.rm = TRUE)

# Mean of adulthood age
values_adulthood_ITM2B <- sapply(ITM2B[, grepl("^(21.yrs|23.yrs|30.yrs|36.yrs|37.yrs|40.yrs)", names(ITM2B))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
mean_adulthood_ITM2B <- mean(values_adulthood_ITM2B, na.rm = TRUE)
mean_adulthood_ITM2B <- round(mean_adulthood_ITM2B, 4)
mean_adulthood_ITM2B
# Standard deviation of adulthood age
sd_adulthood_ITM2B <- sd(values_adulthood_ITM2B, na.rm = TRUE)


# Create df prenatal
vector_prenatal <- rep("Prenatal", 237)
prenal_df <- data.frame(age = vector_prenatal, value = values_prenatal_ITM2B) 

# Create df childhood
vector_childhood <- rep("Childhood", 96)
childhood_df <- data.frame(age = vector_childhood, value = values_childhood_ITM2B)

# Create df adolescent
vector_adolescent <- rep("Adolescent", 64)
adolescent_df <- data.frame(age = vector_adolescent, value = values_adolescent_ITM2B)

# Create df adulthood
vector_adulthood <- rep("Adulthood", 93)
adulthood_df <- data.frame(age = vector_adulthood, value = values_adulthood_ITM2B)

# Concatenate data frames
library(dplyr)
merge_df <- full_join(prenal_df, childhood_df)
merge_df <- full_join(merge_df, adolescent_df)
merge_df <- full_join(merge_df, adulthood_df)

### Krusral - Wallis
library(ggpubr)
library(rstatix)
# Krusral - Wallis test
res.kruskal <- kruskal.test(value ~ age, data = merge_df)
# Effect size
kruskal_effsize(value ~ age, data = merge_df)

# Edit the order of the x-axis in the graph
merge_df$age <- factor(merge_df$age, levels = c("Prenatal", "Childhood", "Adolescent", "Adulthood"))

# Multiple pairwise-comparisons: Pairwise comparisons using Dunn’s test:
pwc <- merge_df %>% 
  dunn_test(value ~ age, p.adjust.method = "bonferroni") 
pwc

# Bar chart ---------------------------------------------------------------
pwc <- pwc %>% add_xy_position(x = "age")
ggboxplot(merge_df, x = "age", y = "value") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    title = "Integral membrane protein 2B",
    x = "Age group", y = "log2 RPKM") +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 25, 
                              hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(t = 15)),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20,
                              color = "black",
                              face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10))
  )
