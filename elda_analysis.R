library(statmod)

# Open PDF device
pdf("Limiting_Dilution_Analysis_Results.pdf", width = 8.5, height = 11)

################################################################################
# ELDA (Extreme Limiting Dilution Analysis) Data Import Script
################################################################################
# Description: This script reads ELDA experimental data from a text file
#              provided as a command line argument and loads it into a dataframe
#
# Usage: Rscript script_name.R <path_to_data_file.txt>
#
# Input file format: Tab-delimited or comma-delimited text file with columns:
#   - cells: Number of cells plated per well
#   - positive: Number of wells with positive growth/response
#   - tested: Total number of wells tested
#   - group: Experimental group/condition identifier
#
# Author: Khoi Huynh
# Date: 1/14/2026
# Version: 1.0
################################################################################

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if file path argument was provided
if (length(args) == 0) {
  stop("Error: No input file specified.\nUsage: Rscript script_name.R <path_to_data_file.txt>")
}

# Store the file path from the first argument
input_file <- args[1]

# Check if file exists
if (!file.exists(input_file)) {
  stop(paste("Error: File not found -", input_file))
}

# Read the data file into a dataframe
# Assumes tab-delimited format; change sep parameter for other delimiters
elda_data <- read.table(
  input_file,
  header = TRUE,           # First row contains column names
  sep = "\t",              # Tab-delimited (use "," for CSV)
  stringsAsFactors = FALSE # Keep strings as character vectors
)

# Display summary of loaded data
cat("\nData successfully loaded from:", input_file, "\n")
cat("Number of rows:", nrow(elda_data), "\n")
cat("Number of columns:", ncol(elda_data), "\n\n")
cat("Column names:", paste(names(elda_data), collapse = ", "), "\n\n")
cat("First few rows:\n")
print(head(elda_data))

elda_data$group <- factor(elda_data$group)

# RUN MAIN ANALYSIS
elda_result <- elda(
  response = elda_data$positive,
  dose = elda_data$cells,
  tested = elda_data$tested,
  group = elda_data$group,
  observed = TRUE,
  test.unit.slope = TRUE
)

# Helper function to draw a nice table
draw_table <- function(table_x, table_y, headers, data_rows, col_widths, cell_height = 0.03) {
  n_rows <- length(data_rows) + 1  # +1 for header
  n_cols <- length(headers)
  total_width <- sum(col_widths)
  
  # Outer border
  rect(table_x, table_y - n_rows*cell_height, table_x + total_width, table_y, lwd = 2)
  
  # Horizontal lines
  for (i in 0:n_rows) {
    segments(table_x, table_y - i*cell_height, table_x + total_width, table_y - i*cell_height, lwd = 1)
  }
  
  # Vertical lines
  x_pos <- table_x
  for (i in 0:n_cols) {
    segments(x_pos, table_y, x_pos, table_y - n_rows*cell_height, lwd = 1)
    if (i < n_cols) {
      x_pos <- x_pos + col_widths[i + 1]
    }
  }
  
  # Headers
  x_pos <- table_x
  for (i in 1:n_cols) {
    text(x_pos + col_widths[i]/2, table_y - cell_height/2, headers[i], cex = 0.75, font = 2)
    x_pos <- x_pos + col_widths[i]
  }
  
  # Data rows
  for (i in 1:length(data_rows)) {
    row_y <- table_y - (i+1)*cell_height + cell_height/2
    x_pos <- table_x
    for (j in 1:length(data_rows[[i]])) {
      text(x_pos + col_widths[j]/2, row_y, data_rows[[i]][j], cex = 0.65)
      x_pos <- x_pos + col_widths[j]
    }
  }
  
  # Return the bottom y position
  return(table_y - n_rows*cell_height)
}

# CREATE TEXT OUTPUT FOR FIRST PAGE
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

y_pos <- 0.98
line_height <- 0.020

# Helper function to add text
add_text <- function(text, size = 0.75, font = 1, x = 0.05) {
  text(x, y_pos, text, adj = c(0, 1), cex = size, family = "mono", font = font)
  y_pos <<- y_pos - line_height
}

# Header information
add_text('The value of the confidence choice entered was "0.95"')
add_text('The value of the observed choice is "TRUE"')
add_text('The value of the test_unit_slope choice is "TRUE"')
add_text('The value of the test_difference choice is "TRUE"')
add_text('Limiting dilution data entered.')
y_pos <- y_pos - line_height

# Data table - using nice table format
table_y <- y_pos
headers <- c("Counter", "Dose", "Tested", "Response", "Group")
col_widths <- c(0.08, 0.06, 0.07, 0.09, 0.35)

data_rows <- list()
for (i in 1:nrow(elda_data)) {
  data_rows[[i]] <- c(
    as.character(i),
    as.character(elda_data$cells[i]),
    as.character(elda_data$tested[i]),
    as.character(elda_data$positive[i]),
    as.character(elda_data$group[i])
  )
}

y_pos <- draw_table(0.05, table_y, headers, data_rows, col_widths, cell_height = 0.025)
y_pos <- y_pos - line_height * 1.5

add_text(sprintf("The number of lines of data entered = %d", nrow(elda_data)))
add_text("Plot results has been checked")
y_pos <- y_pos - line_height

# Confidence intervals table
add_text("Confidence intervals for")
add_text("1/(stem cell frequency)")
y_pos <- y_pos - line_height

table_y <- y_pos
headers <- c("Group", "Lower", "Estimate", "Upper")
col_widths <- c(0.30, 0.10, 0.10, 0.10)

ci_table <- elda_result$CI
data_rows <- list()
for (i in 1:nrow(ci_table)) {
  data_rows[[i]] <- c(
    rownames(ci_table)[i],
    sprintf("%.2f", ci_table[i, "Lower"]),
    sprintf("%.2f", ci_table[i, "Estimate"]),
    sprintf("%.2f", ci_table[i, "Upper"])
  )
}

y_pos <- draw_table(0.05, table_y, headers, data_rows, col_widths, cell_height = 0.028)
y_pos <- y_pos - line_height * 1.5

# Overall test
text(0.05, y_pos, "Overall test for", adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "differences in stem", adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "cell frequencies", adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "between any of the", adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "groups", adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height * 1.2

# Overall test table
table_y <- y_pos
headers <- c("Chisq", "DF", "P.value")
col_widths <- c(0.08, 0.05, 0.12)

overall_chisq <- round(as.numeric(elda_result$test.difference[1]), 1)
overall_df <- as.numeric(elda_result$test.difference[3])
overall_p <- formatC(as.numeric(elda_result$test.difference[2]), format="e", digits=2)

data_rows <- list(
  c(as.character(overall_chisq), as.character(overall_df), overall_p)
)

y_pos <- draw_table(0.05, table_y, headers, data_rows, col_widths, cell_height = 0.025)
y_pos <- y_pos - line_height * 1.5

# Pairwise header
text(0.05, y_pos, "Pairwise tests for differences in stem cell frequencies", 
     adj = c(0, 1), cex = 0.75, family = "mono")

# Start second page for pairwise comparisons
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
par(mar = c(0, 0, 0, 0))

y_pos <- 0.95

# Pairwise comparisons table
groups <- levels(elda_data$group)
pw_data <- list()
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    sub_data <- subset(elda_data, group %in% c(groups[i], groups[j]))
    pw_result <- limdil(sub_data$positive, sub_data$cells, sub_data$tested, 
                        sub_data$group, observed = TRUE)
    pw_data[[length(pw_data) + 1]] <- c(
      groups[i],
      groups[j],
      sprintf("%.2f", as.numeric(pw_result$test.difference[1])),
      as.character(as.numeric(pw_result$test.difference[3])),
      formatC(as.numeric(pw_result$test.difference[2]), format="e", digits=2)
    )
  }
}

table_y <- y_pos
headers <- c("Group 1", "Group 2", "Chisq", "DF", "Pr(>Chisq)")
col_widths <- c(0.25, 0.25, 0.08, 0.05, 0.12)

y_pos <- draw_table(0.05, table_y, headers, pw_data, col_widths, cell_height = 0.03)
y_pos <- y_pos - line_height * 2

# Goodness of fit
text(0.05, y_pos, "Goodness of fit tests. These test whether the log-dose slope", 
     adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "equals 1. Rejection of the tests may be due either to batch effects", 
     adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "(heterogeneity in the stem cell frequencies or assay success rate)", 
     adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height
text(0.05, y_pos, "or to a failure of the stem cell hypothesis.", 
     adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height * 1.5

slope_text <- sprintf("Estimated slope is %.2f", 
                      elda_result$test.slope.wald["Estimate"])
text(0.05, y_pos, slope_text, adj = c(0, 1), cex = 0.75, family = "mono")
y_pos <- y_pos - line_height * 1.5

# Goodness of fit table
lr_chisq <- elda_result$test.slope.lr["z value"]^2
lr_p <- elda_result$test.slope.lr["Pr(>|z|)"]
score_chisq <- elda_result$test.slope.score.dose["z value"]^2
score_p <- elda_result$test.slope.score.dose["Pr(>|z|)"]

table_y <- y_pos
headers <- c("Test", "Chisq", "DF", "P Value")
col_widths <- c(0.42, 0.08, 0.05, 0.10)

data_rows <- list(
  c("Likelihood ratio test of single-hit model", 
    sprintf("%.1f", lr_chisq), "1", sprintf("%.4f", lr_p)),
  c("Score test of heterogeneity", 
    sprintf("%.2f", score_chisq), "1", sprintf("%.4f", score_p))
)

y_pos <- draw_table(0.05, table_y, headers, data_rows, col_widths, cell_height = 0.028)

# New page for plot
plot.new()
par(mar = c(5, 4, 4, 2) + 0.1)
plot(elda_result, main = "")

# Close PDF device
dev.off()

cat("PDF created successfully: Limiting_Dilution_Analysis_Results.pdf\n")
