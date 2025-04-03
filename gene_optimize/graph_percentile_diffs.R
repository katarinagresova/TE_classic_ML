# Load necessary libraries
library(ggplot2)

# Function to create and save a graph
create_graph <- function(data_file, output_pdf, graph_title) { # Added graph_title parameter
  # Read the CSV file
  data <- read.csv(data_file)
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Cell_Line = factor(rep(c("Minimizing", "Maximizing"), each = nrow(data)), levels = c("Minimizing", "Maximizing")), # Set factor levels
    Percentile = c(data$percentile_in_min.ing_cell_line, data$percentile_in_max.ing_cell_line),
    Group = rep(1:nrow(data), 2)
  )

    # Calculate the average difference
    data$percentile_diff <- data$percentile_in_max.ing_cell_line - data$percentile_in_min.ing_cell_line
    average_diff <- mean(data$percentile_diff) # Calculate the average difference
  
  # Generate the plot
  p <- ggplot(plot_data, aes(x = Cell_Line, y = Percentile, group = Group, color = as.factor(Group))) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Cell Line", 
      y = "Percentile", 
      title = paste0(graph_title, "\nAvg. Difference = ", round(average_diff, 2))
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10), # Title text size
      axis.title = element_text(size = 10), # Axis title text size
      axis.text = element_text(size = 10)   # Axis tick text size
    )
  
  # Save the plot to a PDF
  ggsave(output_pdf, plot = p, width = 3, height = 3)
}

# Generate graphs for each CSV file
create_graph("results/utr_dgd_codon_op.csv", "utr_dgd_codon_op_graph.pdf", "UTR + DGD Codon Optimization")
create_graph("results/utr_naive_codon_op.csv", "utr_naive_codon_op_graph.pdf", "UTR + Naive Codon Optimization")
create_graph("results/utr_only.csv", "utr_only_graph.pdf", "UTR Optimization")
