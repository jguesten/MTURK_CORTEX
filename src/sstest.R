#### SSTEST ####

### RUN SIMPLE AGE SLOPES FOR EACH DOMAIN ON LME MODEL ###
### (MODIFIED FROM PACKAGE REGHELPER [Hughes, 2020]) ###

# input: trialdata, subjectdata
# output: subjectdata dataframe with added domain-specific performance measures


sstest=function(model){
  call <- model@call
  mdata <- model@frame
  
  int_term <- which.max(attr(terms(model), 'order'))
  # get location of highest interaction term
  #int_vars <- names(which(attr(terms(model), 'factors')[, int_term] == 1))
  int_vars=c('Age','Domain')
  # get location of variables in the interaction
  
  # figure out which variables are categorical
  factor_vars_log <- vapply(int_vars, function(v) {
    is.factor(mdata[, v])
  }, logical(1))
  factor_vars <- names(factor_vars_log)[which(factor_vars_log == 1)]
  
  original_contrasts <- list()
  if (length(factor_vars) > 0) {
    for (i in 1:length(factor_vars)) {
      original_contrasts[[factor_vars[i]]] <- contrasts(mdata[, factor_vars[i]])
    }
  }
  
  # get points at which to test each variable
  factors <- .set_factors(mdata, int_vars, levels)
  
  # create grid of all tests to be done
  grids <- .create_grids(mdata, factors)
  
  form <- format(formula(model))
  
  template <- grids[[1]]
  models <- grids[[2]]
  
  # distinguish between lmer and lmerTest models
  if ('lmerModLmerTest' %in% class(model)) {
    models[, c('Test Estimate', 'Std. Error', 'df', 't value', 'Pr(>|t|)')] <- NA
    coef_cols <- 5
  } else {
    models[, c('Test Estimate', 'Std. Error', 't value')] <- NA
    coef_cols <- 3
  }
  
  est_count <- 1
  
  for (i in 1:nrow(template)) {
    new_form <- form
    test_var_name <- names(template)[which(startsWith(as.character(template[i, ]), 'sstest'))]
    test_var <- mdata[[test_var_name]]
    
    for (j in 1:ncol(template)) {
      vname <- colnames(template)[j]
      if (vname != test_var_name) {
        if (is.factor(mdata[[vname]])) {
          # for factors, we set the contrast, with reference group as
          # the one for that test
          contrasts(mdata[[vname]]) <- contr.treatment(
            levels(mdata[[vname]]),
            base=which(levels(mdata[[vname]]) == template[i, j])
          )
        } else {
          # for continuous, we replace the name of the variable in the
          # formula to shift the 0 point
          new_var <- paste0('I(', vname, ' - ', template[i, j], ')')
          new_form <- gsub(vname, new_var, new_form)
        }
      } else {
        # when testing a factor effect, revert to original contrasts
        # that the user had set
        if (is.factor(mdata[[vname]])) {
          contrasts(mdata[[vname]]) <- original_contrasts[[vname]]
        }
      }
    }
    call[['formula']] <- as.formula(new_form)
    call[['data']] <- quote(mdata)
    new_model <- eval(call)
    
    if (is.factor(test_var)) {
      contr <- original_contrasts[[test_var_name]]
      dummy_names <- paste0(test_var_name, colnames(contr))
      
      estimates <- as.data.frame(
        coef(summary(new_model))[dummy_names, ])
      
      # when only one contrast, the coefficients will be a vector, not a
      # matrix, so estimates ends up transposed
      if (ncol(contr) < 2) {
        estimates <- as.data.frame(t(estimates))
        rownames(estimates) <- 1
      }
      
      for (est in 1:nrow(estimates)) {
        models[est_count,
               (ncol(models)-coef_cols+1):ncol(models)] <- estimates[est, ]
        est_count <- est_count + 1
      }
    } else {
      models[est_count,
             (ncol(models)-coef_cols+1):ncol(models)
      ] <- coef(summary(new_model))[test_var_name, ]
      est_count <- est_count + 1
    }
  }
  class(models) <- c('simple_slopes', 'data.frame')
  return(models)
}

#### SIMPLE SLOPES (REGHELPER) ####
print.simple_slopes <- function(
  x,
  digits=max(3L, getOption('digits') - 3L),
  signif.stars=getOption('show.signif.stars'),
  ...) {
  
  model <- x
  
  if (!is.logical(signif.stars) || is.na(signif.stars)) {
    warning("option \"show.signif.stars\" is invalid: assuming TRUE")
    signif.stars <- TRUE
  }
  
  index <- c('Test Estimate', 'Std. Error', 't value', 'df')
  for (i in index) {
    if (i %in% colnames(model)) {
      model[, i] <- round(model[, i], digits=digits)
    }
  }
  
  if ('Pr(>|t|)' %in% colnames(model)) {
    if (signif.stars) {
      stars <- symnum(as.numeric(model[, 'Pr(>|t|)']), corr=FALSE, na=FALSE,
                      numeric.x=TRUE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                      symbols=c("***", "**", "*", ".", " "))
    }
    
    model[, 'Pr(>|t|)'] <- format.pval(as.numeric(model[, 'Pr(>|t|)']), digits=digits)
    
    if (signif.stars) {
      model[, 'Sig.'] <- as.character(stars)
    }
  }
  
  print.data.frame(model, quote=FALSE, right=TRUE, na.print='NA')
  invisible(model)
}


#' Print simple slopes.
#' 
#' \code{print} method for class "\code{simple_slopes_lme4}".
#' 
#' @param x An object of class "\code{simple_slopes_lme4}", usually, a result
#'   of a call to \code{\link{simple_slopes}}.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' @seealso \code{\link{simple_slopes}}
#' @export
print.simple_slopes_lme4 <- function(
  x,
  digits=max(3L, getOption('digits') - 3L),
  ...) {
  
  model <- x
  
  index <- c('Test Estimate', 'Std. Error', 't value')
  for (i in index) {
    model[, i] <- round(model[, i], digits=digits)
  }
  
  print.data.frame(model, quote=FALSE, right=TRUE, na.print='NA')
  invisible(model)
}


#' Create grid of simple slope tests.
#' 
#' Helper function creates a data frame with a list of all simple slope tests to
#' be done.
#' 
#' @param data Data from the linear model being tested.
#' @param factors List of all points to test for each variable in the
#'   interaction.
#' @return A data frame with one line for each simple slope test, indicating
#'   what point each variable is set to in the test.
#' @noRd
.create_grids <- function(data, factors) {
  grid <- with(data, expand.grid(factors))
  
  # we only want to use the models that are testing a single variable
  find_tests <- apply(grid, 1, function(x) {
    length(which(suppressWarnings(x == 'sstest'))) == 1
  })
  grid <- grid[find_tests, ]
  
  rownames(grid) <- seq(to=nrow(grid))
  
  # remove factor levels from factor variables
  grid <- as.data.frame(lapply(grid, function(x) {
    as.character(x)
  }), stringsAsFactors=FALSE)
  
  new_grid <- grid
  
  # look for factor variables with more than 2 levels -- they need extra rows
  for (var in names(factors)) {
    variable <- data[[var]]
    if (is.factor(variable) && length(levels(variable)) > 2) {
      find_rows <- which(new_grid[, var] == 'sstest')
      contr <- contrasts(variable)
      
      # count up number of times we should be repeating each row
      num_rep <- ifelse(1:nrow(new_grid) %in% find_rows, ncol(contr), 1)
      rep_index <- rep(row.names(new_grid), num_rep)
      
      new_grid <- new_grid[rep_index, ]
      
      # change rownames to letter subscripts, e.g., '4a', '4b'
      dupe_values <- data.frame(table(rep_index))
      dupe_values <- dupe_values[dupe_values$Freq > 1, ]
      for (i in dupe_values$rep_index) {
        indices <- which(rep_index == as.character(i))
        new_grid[indices, var] <- paste0('sstest (', colnames(contr), ')')
      }
      rownames(new_grid) <- 1:nrow(new_grid)
    }
  }
  
  return(list(grid, new_grid))
}

#' Find points to test variables.
#' 
#' Helper function calculates points at which to test slopes of a variable. For
#' categorical variables, this includes all its levels. For continuous
#' variables, this is -1 SD, the mean, and +1 SD.
#' 
#' @param model The linear model being tested.
#' @param int_vars The variables involved in the interaction being tested.
#' @param user_levels Any user-specified levels to be used instead of the
#'   defaults.
#' @param sstest Logical. Whether or not to insert 'sstest' as a factor level.
#' @return A list with all the factor points for each variable in the
#'   interaction.
#' @noRd
.set_factors <- function(model, int_vars, user_levels=NULL, sstest=TRUE) {
  factors <- list()
  for (term in int_vars) {
    term_data <- model[[term]]
    if (!is.null(user_levels) && term %in% names(user_levels)) {
      # if user specified levels, use those
      factors[[term]] <- user_levels[[term]]
    } else {
      if (class(term_data) == 'factor') {
        # factors are plotted at all levels
        if (sstest == TRUE) {
          factors[[term]] <- c('sstest', levels(term_data))
        } else {
          factors[[term]] <- levels(term_data)
        }
      } else {
        # continuous vars are plotted at -1SD, mean, and +1 SD
        if (sstest == TRUE) {
          factors[[term]] <- c(
            'sstest',
            .offset_point(term_data, -1),
            .offset_point(term_data, 0),
            .offset_point(term_data, 1)
          )
        } else {
          factors[[term]] <- c(
            .offset_point(term_data, -1),
            .offset_point(term_data, 0),
            .offset_point(term_data, 1)
          )
        }
      }
    }
  }
  return(factors)
}


#' Calculate points for testing simple slopes
#' 
#' Helper function calculates point at which to set a variable, based on its
#' mean and standard deviation.
#' 
#' @param var A vector of the variable to test.
#' @param offset_sd The value, in standard deviations, at which you wish to
#'   offset the variable.
#' @param digits Number of decimal places to round value.
#' @return The value of the variable at the offset point. For example, for a
#'   standard normal distribution, an offset point of -1 would return a value of
#'   -1.
#' @noRd
.offset_point <- function(var, offset_sd, digits=6) {
  point <- mean(var, na.rm=TRUE) + offset_sd * sd(var, na.rm=TRUE)
  return(round(point, digits))
}
