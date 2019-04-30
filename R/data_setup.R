read_cdc_data <- function(data_to_read, region = FALSE){
  #' Read CDC Data
  #'
  #' Reads in CDC data from file and puts data into correct format.
  #'
  #' @param data_to_read A string, path to data .csv file.
  #' @param region \code{TRUE} if incorporating interaction between regions,
  #'   \code{FALSE} otherwise.
  #' @return Dataframe with cleaned up CDC data.
  #' @examples dat <- read_cdc_data("D:/School/TB_Research/data/cdc_cases.csv",
  #'   region = FALSE)
  #' @section Warning:
  #' Data must contain 'cases' column and have time as the first column.

  data <- read.csv(data_to_read)

  data$cases <- as.numeric(as.character(data$cases))
  data$cases[is.na(data$cases)] <- 2
  data$epiweek <- data[[1]]
  data$weeknum <- as.numeric(data[[1]])

  if (region == TRUE){
    data$region <- as.numeric(data$region)
  }

  return(data)

}
