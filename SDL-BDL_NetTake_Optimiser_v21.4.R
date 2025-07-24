# STREAMLINED SDL POPULATION OPTIMIZATION SOLVER
# Cleaned version with only essential components

# Required Libraries
library(SWTools)
library(lubridate)
library(tidyverse) 
library(zoo)
library(gridExtra)

# ============================================================================
# CONFIGURATION
# ============================================================================

scenVars <- list( 
  target = 37.80,                                    # Target net take (GL)
  pFolder = "E:/ACT_Uplift/Models/PopulationV3",     # Main working folder
  scenFile = "Population_Function_Setup.csv",       # Scenario setup file
  sourceProjectName = "UMB_ACT_POPV3.rsproj",      # Source model file
  outputFolder = "outputs",                          # Output directory
  inputSetFile = "ACT_PopulationV3_InputSet.txt",   # Input set file
  summaryFilename = "SDL_Opt_Crops_V3_Veneer_AnnualMeans.csv",  # Results file
  minPop = 800000,                                   # Population lower bound
  maxPop = 1200000,                                  # Population upper bound
  startDate = "01/07/1890",                          # Model start date
  endDate = "30/06/2024",                            # Model end date
  startWindow = "01/07/1895",                        # Analysis window start
  endWindow = "30/06/2024",                          # Analysis window end
  baseURL = "http://localhost:9876"                  # Veneer server URL
)

# Load existing functions
source("E:/R/shSDL/SDL_Opt_Functions.R")

# Setup files and folders
fileSetup(scenVars)

# Load scenario file
scenario_file_path <- paste0(scenVars$pFolder, "/", scenVars$scenFile)
if (file.exists(scenario_file_path)) {
  sf <- read.csv(scenario_file_path)
  cat("Scenario file loaded successfully.\n")
} else {
  stop("Scenario file not found: ", scenario_file_path)
}

# ============================================================================
# SDL SOLVER CLASS
# ============================================================================

create_sdl_solver <- function(scenVars, sf) {
  
  # Solver configuration
  config <- list(
    target = scenVars$target,
    tolerance = 0.01,
    max_iterations = 12,
    bounds = c(scenVars$minPop, scenVars$maxPop)
  )
  
  # Solver state
  state <- list(
    results = data.frame(
      Run = integer(),
      Population = numeric(),
      NetTake = numeric(),
      Residual = numeric()
    ),
    converged = FALSE,
    method_used = "secant"
  )
  
  # Evaluate population using Veneer
  evaluate_population <- function(population, iteration) {
    cat("Iteration", iteration, ": Population =", format(population, big.mark = ","))
    
    # Write input file and run model
    writeInputSetFile3(Population = population, scenVars, sf)
    VeneerRunSource(
      StartDate = scenVars$startDate,
      EndDate = scenVars$endDate,
      baseURL = scenVars$baseURL
    )
    
    # Get and process results
    tryCatch({
      # Get time series data
      df <- fortify.zoo(VeneerGetTS(
        TSURL = "/runs/latest/location/Functions/element/Functions/variable/Functions@ACT_Net_Take@$f_Net_Take"
      ))
      
      # Process data
      colnames(df) <- c("Date", "NetTake")
      df$Date <- date(df$Date)
      
      # Apply analysis window
      df <- df %>% 
        filter(
          Date >= date(dmy(scenVars$startWindow)),
          Date <= date(dmy(scenVars$endWindow))
        ) %>%
        mutate(
          Mon = month(Date),
          Year = year(Date),
          WY = getWaterYearAnnual(.)
        )
      
      # Calculate annual mean net take
      annual_means <- df %>%
        group_by(WY) %>%
        summarise(annual_take = sum(NetTake, na.rm = TRUE), .groups = 'drop')
      
      net_take_gl <- mean(annual_means$annual_take) / 1000
      residual <- abs(config$target - net_take_gl)
      
      # Store result
      new_result <- data.frame(
        Run = iteration,
        Population = population,
        NetTake = net_take_gl,
        Residual = residual
      )
      
      state$results <<- rbind(state$results, new_result)
      
      cat(", Net Take =", round(net_take_gl, 2), "GL, Residual =", 
          round(residual, 4), "GL\n")
      
      # Check convergence
      if (residual < config$tolerance) {
        state$converged <<- TRUE
      }
      
      return(new_result)
      
    }, error = function(e) {
      stop("Failed to get results from Veneer: ", e$message)
    })
  }
  
  # Secant method for next population estimate
  get_next_population_secant <- function() {
    n <- nrow(state$results)
    if (n < 2) return(NULL)
    
    # Use last two points
    x1 <- state$results$Population[n-1]
    x2 <- state$results$Population[n]
    f1 <- state$results$NetTake[n-1] - config$target
    f2 <- state$results$NetTake[n] - config$target
    
    # Check for numerical issues
    if (abs(f2 - f1) < 1e-10) return(NULL)
    
    # Secant formula
    x_new <- x2 - f2 * (x2 - x1) / (f2 - f1)
    
    # Apply bounds
    x_new <- max(config$bounds[1], min(config$bounds[2], x_new))
    
    return(round(x_new, 0))
  }
  
  # Newton-Raphson method (alternative)
  get_next_population_newton <- function() {
    n <- nrow(state$results)
    if (n < 2) return(NULL)
    
    # Use last two points for numerical derivative
    x1 <- state$results$Population[n-1] 
    x2 <- state$results$Population[n]
    y1 <- state$results$NetTake[n-1]
    y2 <- state$results$NetTake[n]
    
    # Check for numerical issues
    if (abs(x2 - x1) < 1e-6) return(NULL)
    
    # Numerical derivative
    dydx <- (y2 - y1) / (x2 - x1)
    if (abs(dydx) < 1e-10) return(NULL)
    
    # Newton step
    current_error <- config$target - y2
    x_new <- x2 + current_error / dydx
    
    # Limit step size to prevent wild jumps
    max_step <- 100000
    step <- x_new - x2
    if (abs(step) > max_step) {
      step <- sign(step) * max_step
      x_new <- x2 + step
    }
    
    # Apply bounds
    x_new <- max(config$bounds[1], min(config$bounds[2], x_new))
    
    return(round(x_new, 0))
  }
  
  # Generate initial population estimates
  get_initial_populations <- function() {
    range_size <- config$bounds[2] - config$bounds[1]
    pop1 <- config$bounds[1] + round(range_size * 0.35, 0)
    pop2 <- config$bounds[1] + round(range_size * 0.75, 0)
    return(c(pop1, pop2))
  }
  
  # Main solver function
  solve <- function(method = "secant") {
    # Print header
    cat("\n", rep("=", 60), "\n")
    cat("SDL POPULATION OPTIMIZATION\n")
    cat(rep("=", 60), "\n")
    cat("Target Net Take:", config$target, "GL\n")
    cat("Tolerance:", config$tolerance, "GL\n")
    cat("Population Range:", format(config$bounds[1], big.mark = ","), 
        "to", format(config$bounds[2], big.mark = ","), "\n")
    cat("Method:", toupper(method), "\n")
    cat("Analysis Period:", scenVars$startWindow, "to", scenVars$endWindow, "\n")
    cat(rep("=", 60), "\n\n")
    
    start_time <- Sys.time()
    state$method_used <<- method
    
    # Phase 1: Initial sampling
    cat("Phase 1: Initial sampling...\n")
    initial_pops <- get_initial_populations()
    
    for (i in seq_along(initial_pops)) {
      evaluate_population(initial_pops[i], i)
      if (state$converged) {
        cat("Early convergence achieved!\n")
        break
      }
    }
    
    # Phase 2: Iterative refinement
    if (!state$converged) {
      cat("\nPhase 2: Iterative refinement...\n")
      
      iteration <- nrow(state$results) + 1
      
      while (!state$converged && iteration <= config$max_iterations) {
        
        # Get next population estimate
        next_pop <- switch(method,
                           "newton" = get_next_population_newton(),
                           "secant" = get_next_population_secant(),
                           get_next_population_secant()  # default to secant
        )
        
        # Fallback strategy if method fails
        if (is.null(next_pop)) {
          cat("Method failed, using fallback strategy...\n")
          last_result <- tail(state$results, 1)
          error <- config$target - last_result$NetTake
          
          if (error > 0) {
            next_pop <- min(last_result$Population * 1.1, config$bounds[2])
          } else {
            next_pop <- max(last_result$Population * 0.9, config$bounds[1])
          }
          next_pop <- round(next_pop, 0)
        }
        
        # Bounds check
        next_pop <- max(config$bounds[1], min(config$bounds[2], next_pop))
        
        # Evaluate new population
        evaluate_population(next_pop, iteration)
        iteration <- iteration + 1
      }
    }
    
    # Calculate execution time
    end_time <- Sys.time()
    execution_time <- round(end_time - start_time, 2)
    
    # Print results
    final_result <- tail(state$results, 1)
    
    cat("\n", rep("=", 60), "\n")
    cat("OPTIMIZATION COMPLETE\n")
    cat(rep("=", 60), "\n")
    cat("Optimum Population:", format(final_result$Population, big.mark = ","), "\n")
    cat("Final Net Take:", round(final_result$NetTake, 2), "GL\n")
    cat("Final Residual:", round(final_result$Residual, 4), "GL\n")
    cat("Iterations Required:", nrow(state$results), "\n")
    cat("Converged:", ifelse(state$converged, "YES", "NO"), "\n")
    cat("Execution Time:", format(execution_time), "\n")
    cat(rep("=", 60), "\n")
    
    return(state$results)
  }
  
  # Create visualization plots
  create_plots <- function() {
    if (nrow(state$results) == 0) return()
    
    # Fit quadratic for comparison (visualization only)
    if (nrow(state$results) >= 3) {
      d <- data.frame(
        Population = state$results$Population, 
        Residual = state$results$Residual
      )
      
      qm <- lm(Residual ~ poly(Population, 2), data = d)
      
      pop_range <- seq(min(d$Population), max(d$Population), length.out = 100)
      curve_data <- data.frame(
        Population = pop_range,
        Residual = predict(qm, newdata = data.frame(Population = pop_range))
      )
    }
    
    # Main optimization plot
    p1 <- state$results %>%
      ggplot(aes(x = Population, y = Residual)) +
      {if(exists("curve_data")) geom_line(data = curve_data, color = "blue", alpha = 0.6)} +
      geom_point(size = 4, alpha = 0.7, color = "#69b3a2") +
      geom_hline(yintercept = config$tolerance, color = "red", linetype = "dashed", alpha = 0.7) +
      geom_point(data = tail(state$results, 1), 
                 aes(x = Population, y = Residual),
                 size = 6, color = "orange", alpha = 0.8) +
      scale_x_continuous(labels = scales::comma_format()) +
      labs(
        title = "SDL Population Optimization", 
        subtitle = paste("Target:", config$target, "GL | Method:", toupper(state$method_used)),
        x = "Population", 
        y = "Residual (GL)",
        caption = paste("Converged in", nrow(state$results), "iterations")
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    
    # Convergence history
    p2 <- state$results %>%
      ggplot(aes(x = Run, y = Residual)) +
      geom_line(color = "blue", size = 1) +
      geom_point(size = 3, color = "#69b3a2") +
      geom_hline(yintercept = config$tolerance, color = "red", linetype = "dashed") +
      scale_y_log10() +
      labs(
        title = "Convergence History", 
        x = "Iteration", 
        y = "Residual (GL, log scale)"
      ) +
      theme_minimal()
    
    # Net Take relationship
    p3 <- state$results %>%
      ggplot(aes(x = Population, y = NetTake)) +
      geom_point(size = 4, alpha = 0.7, color = "#69b3a2") +
      geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 0.6) +
      geom_hline(yintercept = config$target, color = "red", linetype = "dashed") +
      geom_point(data = tail(state$results, 1), 
                 aes(x = Population, y = NetTake),
                 size = 6, color = "orange", alpha = 0.8) +
      scale_x_continuous(labels = scales::comma_format()) +
      labs(
        title = "Net Take vs Population", 
        x = "Population", 
        y = "Net Take (GL)"
      ) +
      theme_minimal()
    
    # Arrange and save plots
    combined_plot <- grid.arrange(p1, p2, p3, 
                                  layout_matrix = rbind(c(1, 1), c(2, 3)))
    
    output_path <- paste0(scenVars$pFolder, "/", scenVars$outputFolder, "/ACT_SDL_Population.png")
    ggsave(filename = output_path, plot = combined_plot, 
           width = 12, height = 8, dpi = 300)
    
    cat("Plots saved to: ACT_SDL_Population.png\n")
    return(combined_plot)
  }
  
  # Save results to CSV
  save_results <- function() {
    # Main results file
    output_file <- paste0(scenVars$pFolder, "/", scenVars$summaryFilename)
    write_csv(state$results, output_file)
    cat("Results saved to:", scenVars$summaryFilename, "\n")
    
    # Summary file
    summary_data <- data.frame(
      Parameter = c("Target", "Final_Population", "Final_NetTake", "Final_Residual", 
                    "Iterations", "Converged", "Method", "Tolerance"),
      Value = c(
        config$target, 
        tail(state$results$Population, 1),
        round(tail(state$results$NetTake, 1), 3),
        round(tail(state$results$Residual, 1), 4),
        nrow(state$results),
        ifelse(state$converged, "Yes", "No"),
        toupper(state$method_used),
        config$tolerance
      )
    )
    
    summary_file <- paste0(scenVars$pFolder, "/optimization_summary.csv")
    write_csv(summary_data, summary_file)
    cat("Summary saved to: optimization_summary.csv\n")
  }
  
  # Return solver interface
  list(
    solve = solve,
    create_plots = create_plots,
    save_results = save_results,
    get_results = function() state$results,
    set_tolerance = function(tol) { config$tolerance <<- tol },
    set_max_iterations = function(max_iter) { config$max_iterations <<- max_iter }
  )
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("Initializing SDL solver...\n")

# Create solver
solver <- create_sdl_solver(scenVars, sf)

# Optional: Adjust settings
# solver$set_tolerance(0.005)
# solver$set_max_iterations(15)

# Run optimization
results <- solver$solve(method = "secant")

# Generate plots and save results
solver$create_plots()
solver$save_results()

# Display final results
cat("\nFinal Results:\n")
print(results)

cat("\n✅ Optimization completed successfully!\n")