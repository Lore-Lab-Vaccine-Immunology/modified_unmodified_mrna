# deseq_helper_functions.R

##################################DESIGN FORMULA############################################

# Function to create the DESeq2 design formula
create_deseq_formula <- function(main, additional = NULL) {
  if (missing(main) || is.null(main)) stop("Error: 'main' effect must be specified.")
  formula_terms <- c(additional, main)  # Additional terms come before main effect
  return(as.formula(paste("~", paste(formula_terms, collapse = " + "))))
}

############################################################################################
##################################VOlCANO PLOT############################################

generate_volcano_plot <- function(results_list, 
                                  comparison_name = NULL,                
                                  label_list = c(),                      
                                  TEST_TEXT = NULL,                 
                                  DESIGN_FORMULA_TEXT = NULL,  
                                  REFERENCE_TIME = NULL,       
                                  figures_dir = NULL,             
                                  x_limits = c(-8, 12.5),                
                                  y_limits = c(0, 160)) {                
  
  # Use global values if the arguments are not provided explicitly
  TEST_TEXT <- ifelse(is.null(TEST_TEXT), get("TEST_TEXT", envir = .GlobalEnv), TEST_TEXT)
  DESIGN_FORMULA_TEXT <- ifelse(is.null(DESIGN_FORMULA_TEXT), get("DESIGN_FORMULA_TEXT", envir = .GlobalEnv), DESIGN_FORMULA_TEXT)
  REFERENCE_TIME <- ifelse(is.null(REFERENCE_TIME), get("REFERENCE_TIME", envir = .GlobalEnv), REFERENCE_TIME)
  figures_dir <- ifelse(is.null(figures_dir), get("figures.dir", envir = .GlobalEnv), figures_dir)
  
  # If no specific comparison is provided, loop through all comparisons in results_list
  comparison_names <- if (is.null(comparison_name)) names(results_list) else comparison_name
  
  for (comparison_name in comparison_names) {
    
    # Extract the result from the results_list
    res <- results_list[[comparison_name]]
    
    # Convert the result to a data frame
    res_df <- as.data.frame(res)
    
    # # Create a new column to classify significance and direction, including biologically insignificant points
    # res_df <- res_df %>%
    #   mutate(Significance = case_when(
    #     padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
    #     padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    #     padj < 0.05 & abs(log2FoldChange) <= 1 ~ "Low FoldChange",
    #     TRUE ~ "Not Significant"
    #   ))
    
    # Combine Gene_Status and Significance into a new column for coloring
    res_df <- res_df %>%
      mutate(Color_Status = case_when(
        # Gene_Status == "Unique to Group1" ~ "Unique to Group1",
        # Gene_Status == "Unique to Group2" ~ "Unique to Group2",
        # Gene_Status == "Unique to Group3" ~ "Unique to Group3",
        padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
        padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        padj < 0.05 & abs(log2FoldChange) <= 1 ~ "Low FoldChange",
        TRUE ~ "Not Significant"
      ))
    
    # Filter for significant genes that are in the predefined label list (if provided)
    label <- if (length(label_list) > 0) {
      subset(res_df, gene_name %in% label_list & abs(log2FoldChange) >= 1 & padj < 0.05)
    } else {
      data.frame()  # If no label_list, return an empty data frame
    }
    
    # Count the number of upregulated and downregulated genes
    deg_number_up <- nrow(subset(res_df, log2FoldChange >= 1 & padj < 0.05))
    deg_number_down <- nrow(subset(res_df, log2FoldChange <= -1 & padj < 0.05))

    # Create the volcano plot with custom colors (but no labels yet)
    # volcano_plot <- ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    #   geom_point(alpha = 0.7) +
    #   scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", 
    #                                 "Not Significant" = "grey80", "Low FoldChange" = "grey40")) +
    #   geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 1) + 
    #   geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
    #   
    #   # Set axis limits
    #   scale_x_continuous(breaks = seq(x_limits[1], x_limits[2], by = 2), limits = x_limits) +
    #   scale_y_continuous(limits = y_limits) +
    #   
    #   # Add the counts of upregulated and downregulated genes
    #   annotate("text", x = 1.5, y = y_limits[2], label = paste0("Upregulated: ", deg_number_up), color = "red", size = 4, hjust = 0) +
    #   annotate("text", x = x_limits[1], y = y_limits[2], label = paste0("Downregulated: ", deg_number_down), color = "blue", size = 4, hjust = 0) +
    #   
    #   labs(title = paste("Comparing:", comparison_name, "| ", TEST_TEXT),
    #        subtitle = paste("Design Formula:", DESIGN_FORMULA_TEXT, " | Reference Time:", REFERENCE_TIME),
    #        x = "Log2 Fold Change",
    #        y = "-Log10 Adjusted P-value") +
    #   theme_prism() +
    #   theme(panel.grid.minor = element_blank())
    
    # Create the volcano plot with updated colors for Gene_Status
    volcano_plot <- ggplot(data = res_df#, 
                           #aes(x = log2FoldChange, y = -log10(padj), color = Color_Status)
                           ) +
      # geom_point(alpha = 0.7) +
      # geom_jitter(width = 0.3, height = 5, alpha = 0.5) +
      
      # Points that need jitter
      geom_jitter(data = subset(res_df, Color_Status %in% c(#"Unique to Group1", 
                                                            #"Unique to Group2", 
                                                            #"Unique to Group3", 
                                                            "Upregulated", 
                                                            "Downregulated")),
                  aes(x = log2FoldChange, y = -log10(padj), color = Color_Status),
                  width = 0.3, height = 5, alpha = 0.7) +
      
      # Points that remain static
      geom_point(data = subset(res_df, !(Color_Status %in% c(#"Unique to Group1", 
                                                             #"Unique to Group2", 
                                                             #"Unique to Group3", 
                                                             "Upregulated", 
                                                             "Downregulated"))),
                 aes(x = log2FoldChange, y = -log10(padj), color = Color_Status),
                 alpha = 0.5) +
      # Custom colors for unique genes and defaults
      scale_color_manual(values = c(
        "Unique to Group1" = "#7392C0",  # Custom color for Group1
        "Unique to Group2" = "#997FAD",  # Custom color for Group2
        "Unique to Group3" = "#B6D47C",  # Custom color for Group3
        "Upregulated" = "red",           # Default colors
        "Downregulated" = "blue",
        "Not Significant" = "grey80",
        "Low FoldChange" = "grey40"
      )) +
      
      # Vertical and horizontal threshold lines
      geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 1) + 
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
      
      # Set axis limits
      scale_x_continuous(breaks = seq(x_limits[1], x_limits[2], by = 2), limits = x_limits) +
      scale_y_continuous(limits = y_limits) +
      
      # Add counts for upregulated and downregulated genes
      annotate("text", x = 1.5, y = y_limits[2], 
               label = paste0("Upregulated: ", deg_number_up), 
               color = "red", size = 4, hjust = 0) +
      annotate("text", x = x_limits[1], y = y_limits[2], 
               label = paste0("Downregulated: ", deg_number_down), 
               color = "blue", size = 4, hjust = 0) +
      
      # Add plot labels and theme
      labs(title = paste("Comparing:", comparison_name, "| ", TEST_TEXT),
           subtitle = paste("Design Formula:", DESIGN_FORMULA_TEXT, " | Reference Time:", REFERENCE_TIME),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value",
           color = "Gene Status") +   # Legend title
      theme_prism() +
      theme(panel.grid.minor = element_blank())
    
    # Add labels only if the label_list is provided and has valid entries
    if (nrow(label) > 0) {
      volcano_plot <- volcano_plot + 
        geom_label_repel(data = label, aes(x = log2FoldChange, y = -log10(padj), label = gene_name))
    }
    
    # Dynamically generate the filename based on DESIGN_FORMULA_TEXT, reference time, and comparison name
    filename <- paste0(figures_dir, 
                       "volcano_plot_", 
                       comparison_name, 
                       "_", TEST_TEXT,
                       ".pdf")
    
    # Save the plot to a file
    ggsave(filename = filename,
           plot = volcano_plot,
           width = 8, height = 6, dpi = 300)
    
    print(volcano_plot)
  }
}