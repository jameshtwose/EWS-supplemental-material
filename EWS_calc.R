EWS_calc <- function(df, plot_pre=FALSE){
  # Get the id of the user as an integer
  user <- df$user_id[1]
  # remove the user_id column
  df <- subset(df, select = -c(user_id))
  cat("Calculating Early Warning Signals and Cumulative Complexity Peaks for user == ", user, "\n")
  # remove columns with too many NAs
  df <- df[, colSums(is.na(df)) < nrow(df)/3]
  vars <- colnames(df)[2:length(colnames(df))]
  if (plot_pre == TRUE) {
    # plot the time series with missing values highlighted
    op <- par(mfrow = c(ceiling(length(vars)/4),3),mar =c(2,2,2,2))
    l_ply(seq_along(vars), function(c){
      imputeTS::plotNA.distribution(x = as.numeric(df[,c]), 
                                    main=colnames(df)[c], 
                                    xlab = "", ylab = "")
    })
    par(op)
  }
  
  # Impute missing values with Classification And Regression Trees / Random Forests
  # RF and CART return (identical) discrete numbers
  imp.cart  <- mice::mice(df[vars], method = 'cart', remove.constant = TRUE, remove.collinear = TRUE, printFlag = FALSE)
  df_imp  <- mice::complete(imp.cart)
  
  if (plot_pre == TRUE) {
    # Plot the timeseries with imputed values, where NAs used to be, in red
    par(og_par)
    for(c in c(vars)){
      # cat("Classification And Regression Trees\n")
      imputeTS::plotNA.imputations(x.withNA = as.numeric(df[,c]),
                                   x.withImputations = as.numeric(df_imp[,c]),
                                   main = paste(c,"CART"), xlab = "", ylab = "")
    }
  }
  
  # put each column between 0 and 1 using elastic scaler
  elasc_df <- data.frame(apply(df_imp, 2, elascer))
  # dynamic complexity of the variables with the imputed data
  win = 28 # here the window is set to 28 due to the slow changing nature of MS
  dc <- dc_win(elasc_df, win = win, scale_min=0, scale_max=1, doPlot = FALSE, colOrder = NA)
  datesIMP <- df$timestamp
  ccp.caseIMP <- dc_ccp(df_win = dc, alpha_item = 0.001, alpha_time = 0.001)
  
  if (plot_pre == TRUE){
    # Plot the Complexity Resonance Diagram Plot
    plotDC_res(df_win = dc, win = win, colOrder = NA, 
               useTimeVector = datesIMP, timeStamp = "99-01-31",
               title = paste("Complexity Resonance Diagram (CART) user == ", user))
    # Plot the Cumulative Complexity Peak Plot
    plotDC_ccp(df_ccp = ccp.caseIMP, win = win, colOrder = NA, 
               useTimeVector = datesIMP, timeStamp = "99-01-31",
               title = paste("Cumulative Complexity Peak Plot (CART) user == ", user))
  }
  
  user_complexity_df <- data.frame("user_id" = user, 
                                   "timestamp" = datesIMP, 
                                   "dynamic_complexity_sum" = rowSums(dc),
                                   "complexity_peaks" = ccp.caseIMP$sig.peaks)
  
  return(user_complexity_df)
}


