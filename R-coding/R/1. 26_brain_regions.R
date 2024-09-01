setwd("D:/baitapR")
row_data <- read.csv("rows.csv", header = TRUE, stringsAsFactors = FALSE)
columns_data <- read.csv("columns.csv", header = TRUE, stringsAsFactors = FALSE)
expression_data <- read.csv("expression.csv", header = FALSE, stringsAsFactors = FALSE)

# Remove the first column (header) of the expression_data
expression_data <- expression_data[,-1]

unique_values <- unique(columns_data$structure_abbreviation)
print(unique_values)

unique_values_1 <- unique(columns_data$structure_name)
print(unique_values_1)

# Convert the "structure abbreviation" column into a list of rows
column_names <- columns_data$structure_abbreviation

# Add "structure_abbreviation" row to expression_data
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

# 1. Primary motor cortex
# Get all values from columns starting with "M1C"
values_M1C <- sapply(expression_data_1[, grepl("^M1C", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_M1C <- rowMeans(values_M1C, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_M1C <- head(sort(means_per_row_M1C, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_M1C <- head(order(-means_per_row_M1C), 10)

print(top_10_means_value_M1C)
print(top_10_means_rows_M1C)


# Find the gene name of 10 max M1C and save it to top_10_gen_name_M1C
top_10_gen_name_M1C <- character(10)

for (i in 1:10) {
  gen_max_M1C <- expression_data_1[top_10_means_rows_M1C[i], 1]
  top_10_gen_name_M1C[i] <- gen_max_M1C
  cat("Tên gen hàng", top_10_means_rows_M1C[i], "là", gen_max_M1C, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_TCx
sd_top_10_means_rows_M1C <- apply(values_M1C[top_10_means_rows_M1C, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_M1C)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_M1C, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Primary motor cortex", xaxt = "n",
               names.arg = top_10_gen_name_M1C, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_M1C - sd_top_10_means_rows_M1C, x1 = mid, y1 = top_10_means_value_M1C + sd_top_10_means_rows_M1C, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_M1C) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_M1C,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 2. Dorsolateral prefrontal cortex
# Get all values from columns starting with "DFC"
values_DFC <- sapply(expression_data_1[, grepl("^DFC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))

# Calculate the mean of each row
means_per_row_DFC <- rowMeans(values_DFC, na.rm = TRUE)

# Take the 10 means with the largest value
top_10_means_value_DFC <- head(sort(means_per_row_DFC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_DFC <- head(order(-means_per_row_DFC), 10)

print(top_10_means_value_DFC)
print(top_10_means_rows_DFC)


# Find the gene name of 10 max DFC and save it to top_10_gen_name_DFC
top_10_gen_name_DFC <- character(10)

for (i in 1:10) {
  gen_max_DFC <- expression_data_1[top_10_means_rows_DFC[i], 1]
  top_10_gen_name_DFC[i] <- gen_max_DFC
  cat("Tên gen hàng", top_10_means_rows_DFC[i], "là", gen_max_DFC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_DFC
sd_top_10_means_rows_DFC <- apply(values_DFC[top_10_means_rows_DFC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_DFC)



# Bar chart ---------------------------------------------------------------

gray.scale <- seq(0, 1, length = 10)

colors <- gray(gray.scale)

old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_DFC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Dorsolateral prefrontal cortex", xaxt = "n",
               names.arg = top_10_gen_name_DFC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_DFC - sd_top_10_means_rows_DFC, x1 = mid, y1 = top_10_means_value_DFC + sd_top_10_means_rows_DFC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_DFC) * 1.2) -0.45,
     y = par("usr")[3],
     labels = top_10_gen_name_DFC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

# 3. Orbital frontal cortex
# Get all values from columns starting with "OFC"
values_OFC <- sapply(expression_data_1[, grepl("^OFC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_OFC <- rowMeans(values_OFC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_OFC <- head(sort(means_per_row_OFC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_OFC <- head(order(-means_per_row_OFC), 10)

print(top_10_means_value_OFC)
print(top_10_means_rows_OFC)


# Find the gene name of 10 max OFC and save it to top_10_gen_name_OFC
top_10_gen_name_OFC <- character(10)

for (i in 1:10) {
  gen_max_OFC <- expression_data_1[top_10_means_rows_OFC[i], 1]
  top_10_gen_name_OFC[i] <- gen_max_OFC
  cat("Tên gen hàng", top_10_means_rows_OFC[i], "là", gen_max_OFC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_OFC
sd_top_10_means_rows_OFC <- apply(values_OFC[top_10_means_rows_OFC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_OFC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_OFC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Orbital frontal cortex", xaxt = "n",
               names.arg = top_10_gen_name_OFC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_OFC - sd_top_10_means_rows_OFC, x1 = mid, y1 = top_10_means_value_OFC + sd_top_10_means_rows_OFC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_OFC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_OFC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 4. Anterior cingulate cortex
# Get all values from columns starting with "MFC"
values_MFC <- sapply(expression_data_1[, grepl("^MFC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_MFC <- rowMeans(values_MFC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_MFC <- head(sort(means_per_row_MFC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_MFC <- head(order(-means_per_row_MFC), 10)

print(top_10_means_value_MFC)
print(top_10_means_rows_MFC)


# Find the gene name of 10 max MFC and save it to top_10_gen_name_MFC
top_10_gen_name_MFC <- character(10)

for (i in 1:10) {
  gen_max_MFC <- expression_data_1[top_10_means_rows_MFC[i], 1]
  top_10_gen_name_MFC[i] <- gen_max_MFC
  cat("Tên gen hàng", top_10_means_rows_MFC[i], "là", gen_max_MFC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_MFC
sd_top_10_means_rows_MFC <- apply(values_MFC[top_10_means_rows_MFC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_MFC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_MFC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Anterior cingulate cortex", xaxt = "n",
               names.arg = top_10_gen_name_MFC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_MFC - sd_top_10_means_rows_MFC, x1 = mid, y1 = top_10_means_value_MFC + sd_top_10_means_rows_MFC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_MFC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_MFC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 5. Ventrolateral prefrontal cortex
# Get all values from columns starting with "VFC"
values_VFC <- sapply(expression_data_1[, grepl("^VFC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_VFC <- rowMeans(values_VFC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_VFC <- head(sort(means_per_row_VFC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_VFC <- head(order(-means_per_row_VFC), 10)

print(top_10_means_value_VFC)
print(top_10_means_rows_VFC)


# Find the gene name of 10 max VFC and save it to top_10_gen_name_VFC
top_10_gen_name_VFC <- character(10)

for (i in 1:10) {
  gen_max_VFC <- expression_data_1[top_10_means_rows_VFC[i], 1]
  top_10_gen_name_VFC[i] <- gen_max_VFC
  cat("Tên gen hàng", top_10_means_rows_VFC[i], "là", gen_max_VFC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_VFC
sd_top_10_means_rows_VFC <- apply(values_VFC[top_10_means_rows_VFC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_VFC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_VFC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Ventrolateral prefrontal cortex", xaxt = "n",
               names.arg = top_10_gen_name_VFC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_VFC - sd_top_10_means_rows_VFC, x1 = mid, y1 = top_10_means_value_VFC + sd_top_10_means_rows_VFC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_VFC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_VFC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 6. Primary somatosensory cortex
# Get all values from columns starting with "S1C"
values_S1C <- sapply(expression_data_1[, grepl("^S1C", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_S1C <- rowMeans(values_S1C, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_S1C <- head(sort(means_per_row_S1C, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_S1C <- head(order(-means_per_row_S1C), 10)

print(top_10_means_value_S1C)
print(top_10_means_rows_S1C)


# Find the gene name of 10 max S1C and save it to top_10_gen_name_S1C
top_10_gen_name_S1C <- character(10)

for (i in 1:10) {
  gen_max_S1C <- expression_data_1[top_10_means_rows_S1C[i], 1]
  top_10_gen_name_S1C[i] <- gen_max_S1C
  cat("Tên gen hàng", top_10_means_rows_S1C[i], "là", gen_max_S1C, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_S1C
sd_top_10_means_rows_S1C <- apply(values_S1C[top_10_means_rows_S1C, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_S1C)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_S1C, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Primary somatosensory cortex", xaxt = "n",
               names.arg = top_10_gen_name_S1C, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_S1C - sd_top_10_means_rows_S1C, x1 = mid, y1 = top_10_means_value_S1C + sd_top_10_means_rows_S1C, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_S1C) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_S1C,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 7. Parietal neocortex
# Get all values from columns starting with "PCx"
values_PCx <- sapply(expression_data_1[, grepl("^PCx", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_PCx <- rowMeans(values_PCx, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_PCx <- head(sort(means_per_row_PCx, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_PCx <- head(order(-means_per_row_PCx), 10)

print(top_10_means_value_PCx)
print(top_10_means_rows_PCx)


# Find the gene name of 10 max PCx and save it to top_10_gen_name_PCx
top_10_gen_name_PCx <- character(10)

for (i in 1:10) {
  gen_max_PCx <- expression_data_1[top_10_means_rows_PCx[i], 1]
  top_10_gen_name_PCx[i] <- gen_max_PCx
  cat("Tên gen hàng", top_10_means_rows_PCx[i], "là", gen_max_PCx, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_PCx
sd_top_10_means_rows_PCx <- apply(values_PCx[top_10_means_rows_PCx, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_PCx)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_PCx, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Parietal neocortex", xaxt = "n",
               names.arg = top_10_gen_name_PCx, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_PCx - sd_top_10_means_rows_PCx, x1 = mid, y1 = top_10_means_value_PCx + sd_top_10_means_rows_PCx, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_PCx) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_PCx,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 8. Posteroventral parietal cortex
# Get all values from columns starting with "IPC"
values_IPC <- sapply(expression_data_1[, grepl("^IPC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_IPC <- rowMeans(values_IPC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_IPC <- head(sort(means_per_row_IPC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_IPC <- head(order(-means_per_row_IPC), 10)

print(top_10_means_value_IPC)
print(top_10_means_rows_IPC)


# Find the gene name of 10 max IPC and save it to top_10_gen_name_IPC
top_10_gen_name_IPC <- character(10)

for (i in 1:10) {
  gen_max_IPC <- expression_data_1[top_10_means_rows_IPC[i], 1]
  top_10_gen_name_IPC[i] <- gen_max_IPC
  cat("Tên gen hàng", top_10_means_rows_IPC[i], "là", gen_max_IPC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_IPC
sd_top_10_means_rows_IPC <- apply(values_IPC[top_10_means_rows_IPC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_IPC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_IPC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Posteroventral parietal cortex", xaxt = "n",
               names.arg = top_10_gen_name_IPC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_IPC - sd_top_10_means_rows_IPC, x1 = mid, y1 = top_10_means_value_IPC + sd_top_10_means_rows_IPC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_IPC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_IPC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 9. Inferolateral temporal cortex
# Get all values from columns starting with "ITC"
values_ITC <- sapply(expression_data_1[, grepl("^ITC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_ITC <- rowMeans(values_ITC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_ITC <- head(sort(means_per_row_ITC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_ITC <- head(order(-means_per_row_ITC), 10)

print(top_10_means_value_ITC)
print(top_10_means_rows_ITC)


# Find the gene name of 10 max ITC and save it to top_10_gen_name_ITC
top_10_gen_name_ITC <- character(10)

for (i in 1:10) {
  gen_max_ITC <- expression_data_1[top_10_means_rows_ITC[i], 1]
  top_10_gen_name_ITC[i] <- gen_max_ITC
  cat("Tên gen hàng", top_10_means_rows_ITC[i], "là", gen_max_ITC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_ITC
sd_top_10_means_rows_ITC <- apply(values_ITC[top_10_means_rows_ITC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_ITC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_ITC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Inferolateral temporal cortex", xaxt = "n",
               names.arg = top_10_gen_name_ITC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_ITC - sd_top_10_means_rows_ITC, x1 = mid, y1 = top_10_means_value_ITC + sd_top_10_means_rows_ITC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_ITC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_ITC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 10. Posterior superior temporal cortex
# Get all values from columns starting with "STC"
values_STC <- sapply(expression_data_1[, grepl("^STC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_STC <- rowMeans(values_STC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_STC <- head(sort(means_per_row_STC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_STC <- head(order(-means_per_row_STC), 10)

print(top_10_means_value_STC)
print(top_10_means_rows_STC)


# Find the gene name of 10 max STC and save it to top_10_gen_name_STC
top_10_gen_name_STC <- character(10)

for (i in 1:10) {
  gen_max_STC <- expression_data_1[top_10_means_rows_STC[i], 1]
  top_10_gen_name_STC[i] <- gen_max_STC
  cat("Tên gen hàng", top_10_means_rows_STC[i], "là", gen_max_STC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_STC
sd_top_10_means_rows_STC <- apply(values_STC[top_10_means_rows_STC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_STC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_STC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Posterior superior temporal cortex", xaxt = "n",
               names.arg = top_10_gen_name_STC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_STC - sd_top_10_means_rows_STC, x1 = mid, y1 = top_10_means_value_STC + sd_top_10_means_rows_STC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_STC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_STC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 11. Occipital neocortex
# Get all values from columns starting with "Ocx"
values_Ocx <- sapply(expression_data_1[, grepl("^Ocx", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_Ocx <- rowMeans(values_Ocx, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_Ocx <- head(sort(means_per_row_Ocx, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_Ocx <- head(order(-means_per_row_Ocx), 10)

print(top_10_means_value_Ocx)
print(top_10_means_rows_Ocx)


# Find the gene name of 10 max Ocx and save it to top_10_gen_name_Ocx
top_10_gen_name_Ocx <- character(10)

for (i in 1:10) {
  gen_max_Ocx <- expression_data_1[top_10_means_rows_Ocx[i], 1]
  top_10_gen_name_Ocx[i] <- gen_max_Ocx
  cat("Tên gen hàng", top_10_means_rows_Ocx[i], "là", gen_max_Ocx, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_Ocx
sd_top_10_means_rows_Ocx <- apply(values_Ocx[top_10_means_rows_Ocx, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_Ocx)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_Ocx, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Occipital neocortex", xaxt = "n",
               names.arg = top_10_gen_name_Ocx, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_Ocx - sd_top_10_means_rows_Ocx, x1 = mid, y1 = top_10_means_value_Ocx + sd_top_10_means_rows_Ocx, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_Ocx) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_Ocx,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 12. Temporal neocortex
# Get all values from columns starting with "TCx"
values_TCx <- sapply(expression_data_1[, grepl("^TCx", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))
values_TCx <- unname(values_TCx)
values_TCx <- data.frame(Value = values_TCx)

# Create a vector to store the 10 largest values
top_10_max_values <- numeric(10)

# Create a vector to store the row corresponding to the 10 largest values
top_10_max_rows <- character(10)

# Loop through each row of the data
for (i in 1:nrow(expression_data_1)) {
  # Get the largest value in the current row
  max_value_TCx <- max(values_TCx[i, ], na.rm = TRUE)
  
  # Check if the largest value is larger than any value in the top 10
  if (max_value_TCx > min(top_10_max_values)) {
    # Find the position of the smallest value in the top 10
    min_index <- which.min(top_10_max_values)
    
    # Update the largest value and the corresponding row name
    top_10_max_values[min_index] <- max_value_TCx
    top_10_max_rows[min_index] <- rownames(expression_data_1)[i]
  }
}

# Sort top_10_max_values in descending order
sorted_indices <- order(top_10_max_values, decreasing = TRUE)
top_10_max_values <- top_10_max_values[sorted_indices]
top_10_max_rows <- top_10_max_rows[sorted_indices]

# Print the 10 largest values and the corresponding row names
result <- data.frame(Row = top_10_max_rows, Max_Value = top_10_max_values)
cat("Vị trí hàng và giá trị của 10 gen biểu hiện mạnh nhất của vùng não TCx là")
print(result)


# Find the gene name of the 10 max TCx and save it to top_10_gen_name_TCx
top_10_gen_name_TCx <- character(10)

for (i in 1:10) {
  gen_max_TCx <- expression_data_1[top_10_max_rows[i], 1]
  top_10_gen_name_TCx[i] <- gen_max_TCx
  cat("Tên gen hàng", top_10_max_rows[i], "là", gen_max_TCx, "\n")
}

TCx <- data.frame(Gene_name = top_10_gen_name_TCx, Value = top_10_max_values)
TCx

# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
barplot(top_10_max_values, 
        col = colors,
        ylim = c(0, 12),
        ylab = "log2 RPKM", main ="Temporal neocortex", xaxt = "n",
        names.arg = top_10_gen_name_TCx, las=2)

# Add gene name
text(x = (1:length(top_10_max_values) * 1.2) - 0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_TCx,
     xpd = NA,
     ## Rotate the labels by 45 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 13. Primary auditory cortex
# Get all values from columns starting with "A1C"
values_A1C <- sapply(expression_data_1[, grepl("^A1C", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_A1C <- rowMeans(values_A1C, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_A1C <- head(sort(means_per_row_A1C, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_A1C <- head(order(-means_per_row_A1C), 10)

print(top_10_means_value_A1C)
print(top_10_means_rows_A1C)


# Find the gene name of 10 max A1C and save it to top_10_gen_name_A1C
top_10_gen_name_A1C <- character(10)

for (i in 1:10) {
  gen_max_A1C <- expression_data_1[top_10_means_rows_A1C[i], 1]
  top_10_gen_name_A1C[i] <- gen_max_A1C
  cat("Tên gen hàng", top_10_means_rows_A1C[i], "là", gen_max_A1C, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_A1C
sd_top_10_means_rows_A1C <- apply(values_A1C[top_10_means_rows_A1C, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_A1C)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_A1C, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Primary auditory cortex", xaxt = "n",
               names.arg = top_10_gen_name_A1C, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_A1C - sd_top_10_means_rows_A1C, x1 = mid, y1 = top_10_means_value_A1C + sd_top_10_means_rows_A1C, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_A1C) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_A1C,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 14. Primary visual cortex
# Get all values from columns starting with "V1C"
values_V1C <- sapply(expression_data_1[, grepl("^V1C", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_V1C <- rowMeans(values_V1C, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_V1C <- head(sort(means_per_row_V1C, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_V1C <- head(order(-means_per_row_V1C), 10)

print(top_10_means_value_V1C)
print(top_10_means_rows_V1C)


# Find the gene name of 10 max V1C and save it to top_10_gen_name_V1C
top_10_gen_name_V1C <- character(10)

for (i in 1:10) {
  gen_max_V1C <- expression_data_1[top_10_means_rows_V1C[i], 1]
  top_10_gen_name_V1C[i] <- gen_max_V1C
  cat("Tên gen hàng", top_10_means_rows_V1C[i], "là", gen_max_V1C, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_V1C
sd_top_10_means_rows_V1C <- apply(values_V1C[top_10_means_rows_V1C, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_V1C)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_V1C, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Primary visual cortex", xaxt = "n",
               names.arg = top_10_gen_name_V1C, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_V1C - sd_top_10_means_rows_V1C, x1 = mid, y1 = top_10_means_value_V1C + sd_top_10_means_rows_V1C, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_V1C) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_V1C,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 15. Primary motor-sensory cortex
# Get all values from columns starting with "M1C.S1C"
values_M1C.S1C <- sapply(expression_data_1[, grepl("^M1C.S1C", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_M1C.S1C <- rowMeans(values_M1C.S1C, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_M1C.S1C <- head(sort(means_per_row_M1C.S1C, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_M1C.S1C <- head(order(-means_per_row_M1C.S1C), 10)

print(top_10_means_value_M1C.S1C)
print(top_10_means_rows_M1C.S1C)


# Find the gene name of 10 max M1C.S1C and save it to top_10_gen_name_M1C.S1C
top_10_gen_name_M1C.S1C <- character(10)

for (i in 1:10) {
  gen_max_M1C.S1C <- expression_data_1[top_10_means_rows_M1C.S1C[i], 1]
  top_10_gen_name_M1C.S1C[i] <- gen_max_M1C.S1C
  cat("Tên gen hàng", top_10_means_rows_M1C.S1C[i], "là", gen_max_M1C.S1C, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_M1C.S1C
sd_top_10_means_rows_M1C.S1C <- apply(values_M1C.S1C[top_10_means_rows_M1C.S1C, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_M1C.S1C)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_M1C.S1C, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Primary motor-sensory cortex", xaxt = "n",
               names.arg = top_10_gen_name_M1C.S1C, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_M1C.S1C - sd_top_10_means_rows_M1C.S1C, x1 = mid, y1 = top_10_means_value_M1C.S1C + sd_top_10_means_rows_M1C.S1C, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_M1C.S1C) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_M1C.S1C,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 16. Lateral ganglionic eminence
# Get all values from columns starting with "LGE"
values_LGE <- sapply(expression_data_1[, grepl("^LGE", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_LGE <- rowMeans(values_LGE, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_LGE <- head(sort(means_per_row_LGE, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_LGE <- head(order(-means_per_row_LGE), 10)

print(top_10_means_value_LGE)
print(top_10_means_rows_LGE)


# Find the gene name of 10 max LGE and save it to top_10_gen_name_LGE
top_10_gen_name_LGE <- character(10)

for (i in 1:10) {
  gen_max_LGE <- expression_data_1[top_10_means_rows_LGE[i], 1]
  top_10_gen_name_LGE[i] <- gen_max_LGE
  cat("Tên gen hàng", top_10_means_rows_LGE[i], "là", gen_max_LGE, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_LGE
sd_top_10_means_rows_LGE <- apply(values_LGE[top_10_means_rows_LGE, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_LGE)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_LGE, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Lateral ganglionic eminence", xaxt = "n",
               names.arg = top_10_gen_name_LGE, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_LGE - sd_top_10_means_rows_LGE, x1 = mid, y1 = top_10_means_value_LGE + sd_top_10_means_rows_LGE, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_LGE) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_LGE,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 17. Medial ganglionic eminence
# Get all values from columns starting with "MGE"
values_MGE <- sapply(expression_data_1[, grepl("^MGE", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_MGE <- rowMeans(values_MGE, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_MGE <- head(sort(means_per_row_MGE, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_MGE <- head(order(-means_per_row_MGE), 10)

print(top_10_means_value_MGE)
print(top_10_means_rows_MGE)


# Find the gene name of 10 max MGE and save it to top_10_gen_name_MGE
top_10_gen_name_MGE <- character(10)

for (i in 1:10) {
  gen_max_MGE <- expression_data_1[top_10_means_rows_MGE[i], 1]
  top_10_gen_name_MGE[i] <- gen_max_MGE
  cat("Tên gen hàng", top_10_means_rows_MGE[i], "là", gen_max_MGE, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_MGE
sd_top_10_means_rows_MGE <- apply(values_MGE[top_10_means_rows_MGE, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_MGE)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_MGE, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Medial ganglionic eminence", xaxt = "n",
               names.arg = top_10_gen_name_MGE, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_MGE - sd_top_10_means_rows_MGE, x1 = mid, y1 = top_10_means_value_MGE + sd_top_10_means_rows_MGE, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_MGE) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_MGE,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 18. Caudal ganglionic eminence
# Get all values from columns starting with "CGE"
values_CGE <- sapply(expression_data_1[, grepl("^CGE", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_CGE <- rowMeans(values_CGE, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_CGE <- head(sort(means_per_row_CGE, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest valuet
top_10_means_rows_CGE <- head(order(-means_per_row_CGE), 10)

print(top_10_means_value_CGE)
print(top_10_means_rows_CGE)


# Find the gene name of 10 max CGE and save it to top_10_gen_name_CGE
top_10_gen_name_CGE <- character(10)

for (i in 1:10) {
  gen_max_CGE <- expression_data_1[top_10_means_rows_CGE[i], 1]
  top_10_gen_name_CGE[i] <- gen_max_CGE
  cat("Tên gen hàng", top_10_means_rows_CGE[i], "là", gen_max_CGE, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_CGE
sd_top_10_means_rows_CGE <- apply(values_CGE[top_10_means_rows_CGE, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_CGE)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_CGE, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Caudal ganglionic eminence", xaxt = "n",
               names.arg = top_10_gen_name_CGE, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_CGE - sd_top_10_means_rows_CGE, x1 = mid, y1 = top_10_means_value_CGE + sd_top_10_means_rows_CGE, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_CGE) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_CGE,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 19. Upper rhombic lip
# Get all values from columns starting with "URL"
values_URL <- sapply(expression_data_1[, grepl("^URL", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_URL <- rowMeans(values_URL, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_URL <- head(sort(means_per_row_URL, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_URL <- head(order(-means_per_row_URL), 10)

print(top_10_means_value_URL)
print(top_10_means_rows_URL)


# Find the gene name of 10 max URL and save it to top_10_gen_name_URL
top_10_gen_name_URL <- character(10)

for (i in 1:10) {
  gen_max_URL <- expression_data_1[top_10_means_rows_URL[i], 1]
  top_10_gen_name_URL[i] <- gen_max_URL
  cat("Tên gen hàng", top_10_means_rows_URL[i], "là", gen_max_URL, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_URL
sd_top_10_means_rows_URL <- apply(values_URL[top_10_means_rows_URL, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_URL)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_URL, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Upper rhombic lip", xaxt = "n",
               names.arg = top_10_gen_name_URL, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_URL - sd_top_10_means_rows_URL, x1 = mid, y1 = top_10_means_value_URL + sd_top_10_means_rows_URL, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_URL) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_URL,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 20. Mediodorsal nucleus of thalamus
# Get all values from columns starting with "MD"
values_MD <- sapply(expression_data_1[, grepl("^MD", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_MD <- rowMeans(values_MD, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_MD <- head(sort(means_per_row_MD, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_MD <- head(order(-means_per_row_MD), 10)

print(top_10_means_value_MD)
print(top_10_means_rows_MD)


# Find the gene name of 10 max MD and save it to top_10_gen_name_MD
top_10_gen_name_MD <- character(10)

for (i in 1:10) {
  gen_max_MD <- expression_data_1[top_10_means_rows_MD[i], 1]
  top_10_gen_name_MD[i] <- gen_max_MD
  cat("Tên gen hàng", top_10_means_rows_MD[i], "là", gen_max_MD, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_MD
sd_top_10_means_rows_MD <- apply(values_MD[top_10_means_rows_MD, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_MD)


# Bar chart ---------------------------------------------------------------
gray.scale <- seq(0, 1, length = 10)

colors <- gray(gray.scale)

old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_MD, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Mediodorsal nucleus of thalamus", xaxt = "n",
               names.arg = top_10_gen_name_MD, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_MD - sd_top_10_means_rows_MD, x1 = mid, y1 = top_10_means_value_MD + sd_top_10_means_rows_MD, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_MD) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_MD,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 21. Dorsal thalamus
# Get all values from columns starting with "DTH"
values_DTH <- sapply(expression_data_1[, grepl("^DTH", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_DTH <- rowMeans(values_DTH, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_DTH <- head(sort(means_per_row_DTH, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_DTH <- head(order(-means_per_row_DTH), 10)

print(top_10_means_value_DTH)
print(top_10_means_rows_DTH)


# Find the gene name of 10 max DTH and save it to top_10_gen_name_DTH
top_10_gen_name_DTH <- character(10)

for (i in 1:10) {
  gen_max_DTH <- expression_data_1[top_10_means_rows_DTH[i], 1]
  top_10_gen_name_DTH[i] <- gen_max_DTH
  cat("Tên gen hàng", top_10_means_rows_DTH[i], "là", gen_max_DTH, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_DTH
sd_top_10_means_rows_DTH <- apply(values_DTH[top_10_means_rows_DTH, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_DTH)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_DTH, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Dorsal thalamus", xaxt = "n",
               names.arg = top_10_gen_name_DTH, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_DTH - sd_top_10_means_rows_DTH, x1 = mid, y1 = top_10_means_value_DTH + sd_top_10_means_rows_DTH, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_DTH) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_DTH,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 22. Amygdaloid complex
# Get all values from columns starting with "AMY"
values_AMY <- sapply(expression_data_1[, grepl("^AMY", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_AMY <- rowMeans(values_AMY, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_AMY <- head(sort(means_per_row_AMY, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_AMY <- head(order(-means_per_row_AMY), 10)

print(top_10_means_value_AMY)
print(top_10_means_rows_AMY)


# Find the gene name of 10 max AMY and save it to top_10_gen_name_AMY
top_10_gen_name_AMY <- character(10)

for (i in 1:10) {
  gen_max_AMY <- expression_data_1[top_10_means_rows_AMY[i], 1]
  top_10_gen_name_AMY[i] <- gen_max_AMY
  cat("Tên gen hàng", top_10_means_rows_AMY[i], "là", gen_max_AMY, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_AMY
sd_top_10_means_rows_AMY <- apply(values_AMY[top_10_means_rows_AMY, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_AMY)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_AMY, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Amygdaloid complex", xaxt = "n",
               names.arg = top_10_gen_name_AMY, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_AMY - sd_top_10_means_rows_AMY, x1 = mid, y1 = top_10_means_value_AMY + sd_top_10_means_rows_AMY, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_AMY) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_AMY,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 23. Cerebellum
# Get all values from columns starting with "CB"
values_CB <- sapply(expression_data_1[, grepl("^CB", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_CB <- rowMeans(values_CB, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_CB <- head(sort(means_per_row_CB, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_CB <- head(order(-means_per_row_CB), 10)

print(top_10_means_value_CB)
print(top_10_means_rows_CB)


# # Find the gene name of 10 max CB and save it to top_10_gen_name_CB
top_10_gen_name_CB <- character(10)

for (i in 1:10) {
  gen_max_CB <- expression_data_1[top_10_means_rows_CB[i], 1]
  top_10_gen_name_CB[i] <- gen_max_CB
  cat("Tên gen hàng", top_10_means_rows_CB[i], "là", gen_max_CB, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_CB
sd_top_10_means_rows_CB <- apply(values_CB[top_10_means_rows_CB, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_CB)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_CB, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Cerebellum", xaxt = "n",
               names.arg = top_10_gen_name_CB, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_CB - sd_top_10_means_rows_CB, x1 = mid, y1 = top_10_means_value_CB + sd_top_10_means_rows_CB, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_CB) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_CB,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 24. Cerebellar cortex
# Get all values from columns starting with "CBC"
values_CBC <- sapply(expression_data_1[, grepl("^CBC", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_CBC <- rowMeans(values_CBC, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_CBC <- head(sort(means_per_row_CBC, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_CBC <- head(order(-means_per_row_CBC), 10)

print(top_10_means_value_CBC)
print(top_10_means_rows_CBC)


# Find the gene name of 10 max CBC and save it to top_10_gen_name_CBC
top_10_gen_name_CBC <- character(10)

for (i in 1:10) {
  gen_max_CBC <- expression_data_1[top_10_means_rows_CBC[i], 1]
  top_10_gen_name_CBC[i] <- gen_max_CBC
  cat("Tên gen hàng", top_10_means_rows_CBC[i], "là", gen_max_CBC, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_CBC
sd_top_10_means_rows_CBC <- apply(values_CBC[top_10_means_rows_CBC, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_CBC)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_CBC, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Cerebellar cortex", xaxt = "n",
               names.arg = top_10_gen_name_CBC, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_CBC - sd_top_10_means_rows_CBC, x1 = mid, y1 = top_10_means_value_CBC + sd_top_10_means_rows_CBC, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_CBC) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_CBC,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 25. Striatum
# Get all values from columns starting with "STR"
values_STR <- sapply(expression_data_1[, grepl("^STR", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_STR <- rowMeans(values_STR, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_STR <- head(sort(means_per_row_STR, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_STR <- head(order(-means_per_row_STR), 10)

print(top_10_means_value_STR)
print(top_10_means_rows_STR)


# Find the gene name of 10 max STR and save it to top_10_gen_name_STR
top_10_gen_name_STR <- character(10)

for (i in 1:10) {
  gen_max_STR <- expression_data_1[top_10_means_rows_STR[i], 1]
  top_10_gen_name_STR[i] <- gen_max_STR
  cat("Tên gen hàng", top_10_means_rows_STR[i], "là", gen_max_STR, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_STR
sd_top_10_means_rows_STR <- apply(values_STR[top_10_means_rows_STR, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_STR)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_STR, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Striatum", xaxt = "n",
               names.arg = top_10_gen_name_STR, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_STR - sd_top_10_means_rows_STR, x1 = mid, y1 = top_10_means_value_STR + sd_top_10_means_rows_STR, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_STR) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_STR,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# 26. Hippocampus
# Get all values from columns starting with "HIP"
values_HIP <- sapply(expression_data_1[, grepl("^HIP", names(expression_data_1))], function(col) as.numeric(gsub("[^0-9.-]", "", col)))


# Calculate the mean of each row
means_per_row_HIP <- rowMeans(values_HIP, na.rm = TRUE)


# Take the 10 means with the largest value
top_10_means_value_HIP <- head(sort(means_per_row_HIP, decreasing = TRUE), 10)

# Get the row position of the 10 means with the largest value
top_10_means_rows_HIP <- head(order(-means_per_row_HIP), 10)

print(top_10_means_value_HIP)
print(top_10_means_rows_HIP)


# Find the gene name of 10 max HIP and save it to top_10_gen_name_HIP
top_10_gen_name_HIP <- character(10)

for (i in 1:10) {
  gen_max_HIP <- expression_data_1[top_10_means_rows_HIP[i], 1]
  top_10_gen_name_HIP[i] <- gen_max_HIP
  cat("Tên gen hàng", top_10_means_rows_HIP[i], "là", gen_max_HIP, "\n")
}

# Calculate the standard deviation for each row in top_10_means_rows_HIP
sd_top_10_means_rows_HIP <- apply(values_HIP[top_10_means_rows_HIP, , drop = FALSE], 1, sd, na.rm = TRUE)
# Show standard deviation for each row
print(sd_top_10_means_rows_HIP)


# Bar chart ---------------------------------------------------------------
old_par <- par(mar = c(15, 8, 2, 0), xpd = NA)
mid <- barplot(top_10_means_value_HIP, 
               col = colors,
               ylim = c(0, 12),
               ylab = "log2 RPKM", main ="Hippocampus", xaxt = "n",
               names.arg = top_10_gen_name_HIP, las=2)


# Add error bars for standard deviation
arrows(x0 = mid, y0 = top_10_means_value_HIP - sd_top_10_means_rows_HIP, x1 = mid, y1 = top_10_means_value_HIP + sd_top_10_means_rows_HIP, angle = 90, code = 3, length = 0.1)


# Add gene name
text(x = (1:length(top_10_means_value_HIP) * 1.2) -0.45,
     y = par("usr")[3],             
     labels = top_10_gen_name_HIP,
     xpd = NA,
     ## Rotate the labels by 40 degrees.
     srt = 40, adj=1.1,
     cex = 0.95)
# Add x-axis name
text(6.1, -12, "Gen name", cex=1)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
