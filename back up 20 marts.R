# backup 20 marts
#without ploty 



####################################################################################################################################
# ui.r
####################################################################################################################################
# Heatmap 
tabPanel("Heatmap",
         useShinyjs(), # Enable Shiny JavaScript operations for enhanced UI interactions.
         tabsetPanel(
           # Heatmap Visualization tab: Provides UI elements for processing and visualizing heatmap data.
           tabPanel("Heatmap Visualization",
                    actionButton("run_process", "Run Data Processing"), # Button to trigger data processing.
                    actionButton("show_help", "Show User Guide"), # Button to display a user guide or help information.
                    radioButtons("selected_dataset", "Select data frame:", # Radio buttons for selecting the dataset type.
                                 choices = c("Original Data" = "original", "Merged Data" = "merged"),
                                 selected = "original"),
                    
                    
                    textOutput("upload_message"), # Message area for displaying upload status or instructions.
                    column(width = 4,  # Adjust width as necessary for your layout
                           # Place the components that should go in the second column here
                           uiOutput("select_group_ui_heatmap"),
                           tableOutput("numerator_table"),
                           tableOutput("denominator_table"),
                    ),
                    column(width = 4,
                           uiOutput("select_lipid_ui"),
                    ),
                    column(width = 4,  # Adjust width as necessary for your layout
                           # Place the components that should go in the third column here
                           uiOutput("p_value_max_ui"),
                           uiOutput("logFC_input_ui"),
                    ),
                    uiOutput("lipid_display_message"),
                    
                    column(width = 12,
                           plotOutput("heatmapPlot", width = "100%", height = "650px") # Adjust height as necessary
                    ),
                    
                    # Table in its own row, below the plot
                    column(width = 12,
                           dataTableOutput("pValueTable")
                    ),
                    uiOutput("select_group_ui"), # Dynamic UI for selecting groups, populated server-side.
                    tableOutput("selected_group_table"), # Displays selected groups for verification.
           ),
           
           # Data of groups in Heatmap tab: Displays the data behind the groups used in the heatmap.
           tabPanel("Data of groups in Heatmap",
                    uiOutput("table_message"), # Dynamic message about the data table, populated server-side.
                    conditionalPanel(
                      condition = "input.run_process > 0", # Only display the following if data processing has been triggered.
                      tableOutput("groups_table"), # Shows table of groups.
                      textOutput("limitation_notice"), # Notice about any limitations or considerations.
                      textOutput("rows_removed_text"), # Information on data rows removed during processing.
                      tableOutput("raw_data_table") # Displays the raw data table.
                    ),
                    verbatimTextOutput("error_message") # Area to display any error messages.
           )
         )
),

####################################################################################################################################
# functions.r
####################################################################################################################################




# Heatmap
#Grouping the data
process_data <- function(sequence, data) {
  # Initialize an empty list to store the sample names by group and other information
  results <- list()
  
  sample_identifiers <- rownames(sequence)[sequence[, "labels"] == "Sample"]
  groups <- sequence[sample_identifiers, "class"]
  names <- sample_identifiers
  
  # Extract the 'Sample' labels and corresponding 'class' from 'sequence'
  sample_rows <- sequence[sequence$labels == "Sample", ]
  results$sample_rows <- sample_rows  # Store the sample rows in the results list
  
  unique_groups <- unique(sample_rows$class)
  results$unique_groups <- unique_groups  # Store the unique groups in the results list
  
  # Initialize an empty list within results for grouped_samples
  results$grouped_samples <- list()
  
  # Iterate over each group to get the corresponding sample names
  for (group in unique_groups) {
    if (!is.na(group)) {
      # Get the sample names for the current group
      samples_in_group <- sample_rows[sample_rows$class == group, 1]
      
      # Add the sample names to the list, named by their group
      results$grouped_samples[[paste("group", group)]] <- samples_in_group
      results[[paste("Group", group)]] <- samples_in_group  # Store each group separately in the results list
    }
  }
  
  # Check if any of the components are NULL or empty and handle accordingly
  if (length(results$grouped_samples) == 0) {
    results$grouped_samples <- list(group1 = NA)  # Placeholder if no groups found
  }
  
  
  
  # Return the list of results that now includes grouped samples, unique groups, sample rows, and each group
  return(results)
}



# Function to create grouped data frames based on sequence and data
create_grouped_data_frames <- function(sequence, data) {
  # Initialize an empty list to store data frames for each group
  grouped_data_frames <- list()
  
  # Extract the 'Sample' labels and corresponding 'class' from 'sequence'
  sample_rows <- sequence[sequence$labels == "Sample", ]
  unique_groups <- unique(sample_rows$class)
  
  # Iterate over each unique group to create data frames
  for (group in unique_groups) {
    if (!is.na(group)) {
      # Get the sample identifiers for the current group
      sample_identifiers <- rownames(sample_rows)[sample_rows$class == group]
      
      # Find the matching column indices, excluding NA values
      matching_indices <- match(sample_identifiers, colnames(data))
      matching_indices <- matching_indices[!is.na(matching_indices)]
      
      # Check if we have any matching columns at all
      if (length(matching_indices) > 0) {
        # Select only the columns for the current group
        group_data <- data[, matching_indices, drop = FALSE]
        
        # Store the filtered data frame in the list, named by the group
        grouped_data_frames[[paste("group", group)]] <- group_data
      } else {
        warning(paste("Group", group, "contains column names that are not in the data. Skipping this group."))
      }
    }
  }
  
  return(grouped_data_frames)
}




combine_data <- function(sequence, data) {
  
  # Initialize an empty data frame to store the combined data
  combined_data <- data.frame()
  sample_rows <- sequence[sequence$labels == "Sample", ] #358
  unique_groups <- unique(sample_rows$class)
  
  for (group in unique_groups) {
    if (!is.na(group) && group %in% sequence$class) {
      samples_in_group <- sequence[sequence$class == group & sequence$labels == "Sample", 1]
      # Find the matching column indices, excluding NA values
      matching_indices <- match(samples_in_group, colnames(data))
      matching_indices <- matching_indices[!is.na(matching_indices)]  # Exclude NA values
      
      # Check if we have any matching columns at all
      if (length(matching_indices) == 0) {
        next  # Skip this group if no matching columns are found
      }
      
      # Extract the data for the current group's samples
      group_data <- data[, matching_indices, drop = FALSE]
      
      # Combine the extracted data for the current group with the main combined data frame
      # This assumes that all data frames have the same number of rows and row order
      if (ncol(combined_data) == 0) {
        combined_data <- group_data
      } else {
        combined_data <- cbind(combined_data, group_data)
      }
    }
  }
  
  # Now 'combined_data' contains all the sample data from each group without the group prefix
  return(combined_data)
}



calculate_means_for_grouped_data <- function(grouped_data_frames) {
  # Initialize a new list to store the modified data frames
  new_grouped_data_frames <- list()
  
  # Iterate over each group's data frame in the list
  for (group_name in names(grouped_data_frames)) {
    # Clone the current group's data frame to avoid modifying the original
    group_data <- grouped_data_frames[[group_name]]
    
    # Assuming the first column is not numeric and should be excluded from the mean calculation
    # Calculate the mean for each row across all other columns
    means <- rowMeans(group_data[, drop = FALSE], na.rm = TRUE)
    
    # Append the calculated means as a new column to the cloned data frame
    group_data$Mean <- means
    
    # Add the modified data frame to the new list
    new_grouped_data_frames[[group_name]] <- group_data
  }
  
  # Return the new list of grouped data frames with means calculated
  return(new_grouped_data_frames)
}


####      ####      ####      ####      ####      ####      ####      ####      ####      ####
# Anita 
####      ####      ####      ####      ####      ####      ####      ####      ####      ####

# Function to group lipids by their class prefix (e.g., "CAR", "LP", etc.)
group_lipids_by_class <- function(data) {
  # Assuming the first column of 'data' contains the lipid names like "CAR(18:1)"
  lipid_names <- data[[1]]  # Replace 1 with the actual column name or index if different
  
  # Use a regular expression to extract the class prefix from lipid names
  # This matches any consecutive alphabetic characters at the beginning of the string
  lipid_classes <- sub("\\(([0-9]+:[0-9]+)\\).*", "", lipid_names)
  
  # Create a data frame that maps lipid names to their class
  class_mapping <- data.frame(Lipid_Name = lipid_names, Class = lipid_classes, stringsAsFactors = FALSE)
  
  # Optionally, if you want to return a list that names each group by its class
  # names(grouped_data) <- unique(lipid_classes)
  
  return(class_mapping)
}


# Data cleaning

# Function to extract patterns from compound names
# Removes all noise from compound name, so name and length is the only left: eg. going from "CAR 14:1'CAR'[M+H]+" to "CAR 14:1"
extract_pattern <- function(name) {
  
  # Pattern to find first part consisting of letters and numbers with a colon or a letter before the numbers
  pattern <- "([A-Za-z]+\\s[0-9]+:[0-9]+)|([A-Za-z]+\\s[[:alpha:]]?-?[0-9]+:[0-9]+)"
  matches <- regmatches(name, gregexpr(pattern, name))
  
  # Returns the first match, or the hole name if no match
  if (length(matches[[1]]) > 0) {
    return(matches[[1]][1])
  } else {
    return(name)
  }
}

# Function to format strings
# Puts the length and double bonds numbers into a (). Eg "CAR 14:1" to "CAR(14:1)"
format_strings <- function(input_strings) {
  # Use gsub with regular expression to remove all whitespace characters
  formatted_strings <- gsub("\\s+", "", input_strings)
  # Add parentheses around the numbers
  formatted_strings <- gsub("([A-Za-z]*)(\\d+):(\\d+)", "\\1(\\2:\\3)", formatted_strings)
  return(formatted_strings)
}

### pattern_column is that corret? 
# Function to filter rows based on the specified pattern, meaning removes any data that are not on X(C:D) format.
filter_data_by_pattern <- function(data) {
  # Define the regular expression pattern
  pattern <- "^.+\\(\\d+:\\d+\\)$"
  
  # Check if the first column of data matches the pattern
  filtered_data <- data[grepl(pattern, data[[1]]), ]
  
  return(filtered_data)
}


#merge duplicated names of the data
merge_duplicates <- function(data) {
  # Ensure the first column is treated as the Compound Name
  compound_name_col <- names(data)[1]
  
  # Group by the first column and then summarise all other columns by summing
  data_merged <- data %>%
    group_by(.data[[compound_name_col]]) %>%
    summarise(across(everything(), sum, na.rm = TRUE), .groups = 'drop')
  
  return(data_merged)
}

#duplicated names have add _1, _2 and _3 depending on how many duplicates. 
unique_compound_names <- function(data) {
  # Ensure that 'data' is a data frame and has at least one column
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("The input must be a data frame with at least one column.")
  }
  
  # Apply the processing to the first column of 'data'
  data[[1]] <- ave(data[[1]], data[[1]], FUN = function(x) {
    if (length(x) > 1) {
      # Extract the base name without parentheses
      base_name <- sub("\\(.*\\)", "", x)
      # Extract the part within parentheses
      suffix <- sub(".*\\(", "(", x)
      # Combine base name with sequence number and the part within parentheses
      paste0(base_name, "_", seq_along(x), suffix)
    } else {
      x
    }
  })
  
  # Return the modified data
  return(data)
}




####################################################################################################################################
# server.r
####################################################################################################################################

# Heatmap
# Values is used for HTML message display before and after data process
values <- reactiveValues(runProcessClicked = FALSE)

# When bottum clicked in interface, all the following will be processed
observeEvent(input$run_process, {
  values$runProcessClicked <- TRUE
  
  # Accessing sequence and data from active files
  sequence <- rv$sequence[[rv$activeFile]]
  data <- rv$data[[rv$activeFile]]
  
  #HVis before data clean, så det muligt at oprette den her
  
  # Removes anything that are not part of the data of the samples and name. 
  data <- data[, sequence[ , 'labels'] %in% c("Name","Sample")]
  sequence <- sequence[sequence[ , 'labels'] %in% c("Name","Sample"), ]
  
  # Check if the sequence is uploaded before proceeding
  if (is.null(sequence)) {
    return(NULL)  # Stop the observeEvent if no file is uploaded
  }
  
  # Capture the number of rows before filtering
  number_of_rows_before <- nrow(data)
  
  
  # Cleaning the noise of the data, by calling functions
  # Apply the `extract_pattern` function to the first column of the dataset
  # Removes all noise from compound name, so name and length is the only left: eg. going from "CAR 14:1'CAR'[M+H]+" to "CAR 14:1"
  data[, 1] <- sapply(data[, 1], extract_pattern)
  
  # Puts the length and double bonds numbers into a (). Eg "CAR 14:1" to "CAR(14:1)"
  data[, 1] <- sapply(data[, 1], format_strings)
  
  # Function to filter rows based on the specified pattern, meaning removes any data that are not on X(C:D) format.
  data <- filter_data_by_pattern(data)
  
  # This will make it possible to switch between original data and merged data. OG data: using _1 _2 ... _n. Merged will sum the values of the "duplicated" data. 
  if (input$selected_dataset == "original") {
    # Call a function to process the original data
    data <- unique_compound_names(data)
  } else if (input$selected_dataset == "merged") {
    # Call a function to process the merged data
    data <- merge_duplicates(data)
  }
  
  # Some data are in NA, calculations cannot read this, therefore NA values are set low. 
  data[is.na(data)] <- 0.01
  
  
  # For data counting, used in display of how many rows are removed. 
  # Capture the number of rows after filtering
  number_of_rows_after <- nrow(data)
  
  # Calculate the number of rows removed
  rows_removed <- number_of_rows_before - number_of_rows_after
  
  # Output the number of rows removed to the console
  print(paste("Rows removed:", rows_removed))
  
  # Alternatively, you can display this information in the UI using a textOutput
  output$rows_removed_text <- renderText({
    paste("Rows removed after data cleaning are:", rows_removed, ". The removal is be due to the names in the first column of the data file not being in the X(C:D) format. Keep in mind, that a merged data will also count as a removed row.")
  })
  
  # Notification text for the UI
  output$limitation_notice <- renderText({
    data <- rv$data[[rv$activeFile]]  # or however you access your raw data
    paste("Displaying first 25 observations and 5 variables out of", 
          nrow(data), "observations and", ncol(data), "variables in the dataset before data cleaning.")
  })
  
  
  # The following is used in the tab: 'Data of groups in heatmap'.
  # Call the process_data function and use the result directly within observeEvent
  processed_data <- process_data(sequence, data)
  
  output$groups_table <- renderTable({
    if (is.null(processed_data)) {
      return()
    }
    
    # Extract the 'Sample' labels and corresponding 'class' from 'sequence'
    sample_rows <- sequence[sequence$labels == "Sample", ]
    unique_groups <- unique(sample_rows$class)
    
    # Create the dataframe to be displayed as a table
    df_processed_data <- data.frame(
      Group = unique_groups,
      Samples = sapply(unique_groups, function(group) {
        sample_identifiers <- rownames(sample_rows)[sample_rows$class == group]
        paste(sample_identifiers, collapse = ", ")
      })
    )
    # Return the data frame to be rendered as a table
    df_processed_data
  })
  
  # Add this to render the raw data table, being displayed in "Data of groups in Heatmap"
  output$raw_data_table <- renderTable({
    
    # Notify the user of the display limitation and the total size of the data
    total_obs <- nrow(data)
    total_vars <- ncol(data)
    
    #limteing the displayed data
    limited_data <- head(data, 25)[, 1:5]
    
    return(limited_data)
  })
  
  
  # Table output of the table in the tab: 'Heatmap' used for testing. 
  observeEvent(input$run_process, {
    # Process your data here
    processed_results <- process_data(sequence, data)
    grouped_data_frames <- create_grouped_data_frames(sequence, data)
    
    # Add the first column of "data" to each grouped data frame
    compound_names <- data[[1]]  # Extract the first column which contains compound names
    
    # Assuming that each grouped data frame has rows in the same order as "data"
    for (i in seq_along(grouped_data_frames)) {
      grouped_data_frames[[i]] <- cbind(Compound_Name = compound_names, grouped_data_frames[[i]])
    }
    
    # Update the names of the grouped_data_frames if they're not already set
    names(grouped_data_frames) <- paste("Group", seq_along(grouped_data_frames))
    
    # Dynamically generate selectInput for group selection
    output$select_group_ui <- renderUI({
      selectInput("selected_group", "Select Group:",
                  choices = names(grouped_data_frames))  # Use group names as choices
    })
    
    # Dynamically generate table output for the selected group
    output$selected_group_table <- renderTable({
      req(input$selected_group)  # Ensure a group is selected
      grouped_data_frames <- grouped_data_frames[[input$selected_group]]
      if (is.null(grouped_data_frames)) {
        return(data.frame())  # Return an empty data frame if group data is not available
      }
      grouped_data_frames
    })
  })
  
  
  # Heatmap input selection  
  observeEvent(input$run_process, {
    
    # Process your data here
    processed_results <- process_data(sequence, data)
    grouped_data_frames <- create_grouped_data_frames(sequence, data)
    grouped_data_frames_with_means <- calculate_means_for_grouped_data(grouped_data_frames)
    
    # Add the first column of "data" to each grouped data frame
    compound_names <- data[[1]]  # Extract the first column which contains compound names
    
    # Assuming that each grouped data frame has rows in the same order as "data"
    for (i in seq_along(grouped_data_frames_with_means)) {
      grouped_data_frames_with_means[[i]] <- cbind(Compound_Name = compound_names, grouped_data_frames_with_means[[i]])
    }
    
    # Update the names of the grouped_data_frames if they're not already set
    names(grouped_data_frames_with_means) <- paste("Group", seq_along(grouped_data_frames_with_means))
    
    # Create UI for group selection
    output$select_group_ui_heatmap <- renderUI({
      tagList(
        selectInput("selected_group_for_numerator", "Select Group for numerator:",
                    choices = names(grouped_data_frames_with_means)),
        selectInput("selected_group_for_denominator", "Select Group for denominator:",
                    choices = names(grouped_data_frames_with_means))  # Use group names as choices
      )
    })
    
    # Message shown when hovering over Original data and merged data. # Remember to change this outside of the observe event, Search for addTooltip
    observe({
      addTooltip(session, "selected_dataset", 
                 "Choose 'Original Data' to work with the data as it was initially collected. Select 'Merged Data' for a combined and cleaned dataset.", 
                 placement = "bottem", 
                 trigger = "hover")
    })
    
    # Render UI for maximum p-value input
    output$p_value_max_ui <- renderUI({
      numericInput("p_value_max", 
                   "Maximum p-value:", 
                   value = 0.05, 
                   min = 0, 
                   step = 0.01)
    })
    
    # Render UI for logFC input
    output$logFC_input_ui <- renderUI({
      numericInput("logFC_input", 
                   "Enter logFC value:", 
                   value = 2)
    })
    
    # Dynamic p-values depended on interface
    reactiveFilteredData <- reactive({
      # Get the maximum p-value threshold from the input
      p_value_max <- input$p_value_max
      
      # Filter the data based on the maximum p-value
      filtered_data <- numerator_data %>%
        filter(P_Value <= p_value_max)
      
      # Now return the filtered data
      filtered_data
    })
    
    
    
    # Dynamic logFC depended on interface
    reactive_max_logFC <- reactive({
      input$logFC_input  # This will be the positive value
    })
    
    reactive_min_logFC <- reactive({
      -input$logFC_input  # This will be the negative value
    })
    
    
    
    
    # logFC calculation
    reactiveLogFC <- reactive({
      # The required data input for the data handling. 
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      req(reactive_max_logFC(), reactive_min_logFC())
      
      # Define data input, makes it more readable 
      numerator_data <- grouped_data_frames_with_means[[input$selected_group_for_numerator]]
      denominator_data <- grouped_data_frames_with_means[[input$selected_group_for_denominator]]
      
      
      # Ensure there is data to work with
      if(nrow(numerator_data) == 0 || nrow(denominator_data) == 0) {
        return(NULL)
      }
      
      # Calculate logFC
      logFC <- log2((numerator_data$Mean + 1e-6) / (denominator_data$Mean + 1e-6))
      numerator_data$logFC <- logFC
      
      
      
      
      # Filter based on the selected lipid(s), if not 'All'
      if(!"All" %in% input$selected_lipid) {
        numerator_data <- numerator_data[lipid_names$Class %in% input$selected_lipid, ]
      }
      
      # Filter based on the input logFC range
      filtered_data <- numerator_data[numerator_data$logFC >= reactive_min_logFC() & numerator_data$logFC <= reactive_max_logFC(), ]
      return(filtered_data)
    })
    
    
    
    
    # Calculation of p-value, using t-test
    reactiveP_value <- reactive({
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      
      numerator_data <- grouped_data_frames_with_means[[input$selected_group_for_numerator]]
      denominator_data <- grouped_data_frames_with_means[[input$selected_group_for_denominator]]
      
      if(nrow(numerator_data) == 0 || nrow(denominator_data) == 0) {
        return(NULL)
      }
      
      # Initialize a vector to store the p-values
      p_values <- numeric(nrow(numerator_data))
      
      # Loop through each row to perform the t-test
      for (i in 1:nrow(numerator_data)) {
        t_test_result <- t.test(numerator_data[i, -1], denominator_data[i, -1])
        p_values[i] <- t_test_result$p.value
      }
      
      # store the data
      numerator_data$p_values <- p_values
      
      # Filtrate/store the data, so it is ready to display in table in 'Heatmap'.
      filtered_data <- numerator_data %>%
        mutate(p_value = p_values) %>%
        filter(p_value <= input$p_value_max)
      
      return(filtered_data[, c("Compound_Name", "p_values")])
    })
    
    
    # Combine both logFC and p-values into one reactive expression
    reactiveFilteredData <- reactive({
      # Retrieve the filtered datasets based on logFC and p-values
      logFCData <- reactiveLogFC()
      pValuesData <- reactiveP_value()
      
      # Ensure both datasets are not NULL before proceeding
      req(logFCData, pValuesData)
      
      # Combine the datasets to have both logFC and p-value information
      # Assuming both datasets have a 'Compound_Name' column to join on
      combinedData <- merge(logFCData, pValuesData, by = "Compound_Name")
      
      # Now return the combined dataset
      combinedData
    })
    
    
    # Interface of selectinos of lipids to display
    output$select_lipid_ui <- renderUI({
      # Extract the lipid names from first column of the file 'data'
      lipid_names <<- group_lipids_by_class(data)
      
      ###### Above this in the definition of the lipid_names I instead tried to call the function. 
      
      selectizeInput("selected_lipid", "Select lipid(s) to display:",
                     choices = c("All", unique(lipid_names$Class)),
                     multiple = TRUE,
                     options = list(placeholder = 'Choose lipids...',
                                    onInitialize = I('function() { this.setValue(""); }')))
    })
    
    # Reactive expression to track the number of selected lipids or the "All" selection
    selected_lipid_count <- reactive({
      # If "All" is selected, we could set this to a value that causes the default text size to be used
      if ("All" %in% input$selected_lipid) {
        return(Inf)  # 'Inf' is used here as a flag for "All"
      } else {
        return(length(input$selected_lipid))
      }
    })
    
    
    
    
    
    # Heatmap display
    # Use the reactive expression in renderPlotly
    output$heatmapPlot <- renderPlot({
      filtered_data <- reactiveLogFC()
      
      # Take the input from user in interface and change p-value and logFC
      filtered_data <- reactiveFilteredData()
      
      
      # Ensure the data is not NULL and has rows to plot
      req(nrow(filtered_data) > 0)
      
      # Apply any necessary filtering based on logFC
      filtered_data <- filtered_data[filtered_data$logFC >= -2 & filtered_data$logFC <= 2, ]
      
      # Get the count of selected lipids
      num_of_lipids <- selected_lipid_count()
      
      
      # Adjust text size based on the count of selected lipids
      lipid_text_size <- if (num_of_lipids < 5) {
        10  # Smaller number of lipids, larger text
      } else if (num_of_lipids == Inf) {
        4  # "All" is selected, use default text size
      } else {
        5   # More than 5 lipids, smaller text
      }        
      # Ensure compound_names are available. If compound_names were defined earlier,
      # make sure they are accessible here, either as reactive values or as global variables.
      names.mapping <- map_lipid_names(x = filtered_data$Compound_Name)
      
      heatmap_plot <- heatmap_lipidome(
        x = filtered_data[ , c("Compound_Name", "logFC")],
        names.mapping = names.mapping,
        class.facet = "wrap",
        x.names = "Compound_Name",
        fill.limits = c(-2.5, 2.5), # Set limits to include -2.5 to 2.5
        fill.midpoint = 2.5, # Set midpoint explicitly to 2.5
        melt.value.name = "logFC",
        scales = "free"
      ) +
        scale_fill_gradient2(
          low = "blue",       # Color for low values
          mid = "white",     # Color for midpoint values
          high = "red",     # Color for high values
          midpoint = 0,      # Midpoint value, adjust as needed (e.g., the neutral point in your data)
          limit = c(-2.5, 2.5),  # Limits for the scale
          space = "Lab",     # Color space in which to calculate gradient
          name = "logFC"     # Legend title
        ) +
        scale_x_continuous(
          breaks = scales::pretty_breaks(n = lipid_text_size)  # Adjust n to control label frequency
        ) +
        scale_y_continuous(
          breaks = scales::pretty_breaks(n = lipid_text_size)  # Adjust n as needed for the y-axis
        ) 
      # Return the heatmap plot
      heatmap_plot
    })
    
    
    
    #Display af p_Values and logFC values under Heatmap
    output$pValueTable <- renderDataTable({
      # Access the logFC and p_values data
      logFCData <- reactiveLogFC()
      pValuesData <- reactiveP_value()
      
      
      # Ensure both are not NULL before attempting to merge
      req(logFCData, pValuesData)
      
      # Merge the dataframes based on the common "Compound_Name" column
      combinedData <- merge(logFCData, pValuesData, by = "Compound_Name")
      
      # Select only the columns you want to display
      dataTableToShow <<- combinedData[, c("Compound_Name", "logFC", "p_values")]
      
      
      # Round 'logFC' and 'p_values' to the desired number of decimal places
      dataTableToShow$logFC <- round(dataTableToShow$logFC, 5) # 2 decimal places for logFC
      dataTableToShow$p_values <- round(dataTableToShow$p_values, 5) # 4 decimal places for p-values
      
      # Render the selected data in a DataTable
      datatable(dataTableToShow, options = list(pageLength = 10, scrollX = TRUE))
    })
  })
  
  
  # Update the UI message 
  output$table_message <- renderUI({
    if (values$runProcessClicked) {
      HTML('<p>Data processing complete, see tables below.</p>')
    }
  })
  
}) # This finishes the first 'observeEvent' when 'Run data processing' is clicked

# Outside of the observeEvent, based on whether runProcessClicked is TRUE or FALSE, the message display will be placed on this: 
output$table_message <- renderUI({
  if (!values$runProcessClicked) {
    HTML('<p>Make sure sequences file is uploaded, when uploaded: Press "Run Data Processing" to get a display of data</p>')
  }
})

# Outside of the observeEvent, so the message both are shown before and after runProcessClicked is clicked. 
observe({
  addTooltip(session, "selected_dataset", 
             "Choose 'Original Data' to work with the data as it was initially collected. Select 'Merged Data' for a combined and cleaned dataset.", 
             placement = "bottom", 
             trigger = "hover")
})

# User guide inside 'Heatmap'
observeEvent(input$show_help, {
  showModal(modalDialog(
    title = "User Guide for the Heatmap",
    tags$ul(tags$li(tags$b("This Heatmap is designed for comparative lipidomic analysis."))),
    tags$ul("The user guide will explain how to use the Heatmap visualization.",
            tags$li("Upon clicking 'Run Data Processing', the app cleans the data and prepares it for plotting. This process may take a few seconds due to extensive data handling. The 'Data of groups in Heatmap' tab allows you to view the distribution of samples across groups."),
            tags$li("Original Data vs. Merged Data: Switch between 'Original Data' and 'Merged Data' by selecting the desired option and then clicking 'Run Data Processing'. This action will automatically reset the heatmap plot and any selected scales to their default settings."),
            tags$li(tags$b("When using the Heatmap:")),
            tags$li("Select groups for comparison from the dropdown menus."),
            tags$li("The heatmap will automatically update to display the data for the chosen groups."),
            tags$li("To display specific lipids, use 'Select lipid(s) to display:' to type or click your selections. To remove selections, use 'Backspace' on your keyboard. If there are no display of the selected lipids, if may be due to the thresholds of logFC and P-value"),
            tags$li("Adjustments to logFC and P-values are available. Click on 'Enter max logFC value' or 'Maximum P-value:' to enter new thresholds. The application accepts both comma and period as decimal separators. Lipids will be displayed within the specified logFC and P-value ranges. For instance, entering a logFC of 2 will automatically consider a range from -2 to 2."),
            tags$li("For further details on lipid data, scroll down to the table beneath the heatmap."),
            tags$li(tags$b("Data Calculations:")),
            tags$li("The color scale of the heatmap represents log-fold change (logFC) values. LogFC is computed as log2((numerator_data$Mean + 1e-6) / (denominator_data$Mean + 1e-6)), where numerator_data and denominator_data correspond to the groups selected. An offset of 1e-6 is included to avoid division by zero. The logFC scale is set to span from -2.5 to 2.5."),
            tags$li("P-values are determined using the 't.test' function in R, which conducts a statistical comparison between corresponding rows of lipid data from the selected groups."),
            
            
    ),
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})