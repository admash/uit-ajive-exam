# function to filter data 
# data is the initial complete dataset to be filtered
# crit can be variance, IQR or other quantities
# keep is the number of variables to keep in the filtered data
# returns filtered data with keep n of variables
# (Erica Ponzi's code)

filter_top_cols <- function(data, filter, keep, names.only=FALSE){
  filter.list <- data %>% 
    apply(2, filter, na.rm = TRUE ) %>%
    order(decreasing = TRUE) %>%
    extract(1:keep)
  if(names.only){ return(colnames(data)[filter.list]) }
  else          { return(data[,filter.list]) }
}
