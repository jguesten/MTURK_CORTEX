#### function NUMEXTRACT 
# greps numbers from a character variable
# input: x, character vector
# output: y, character vector only containing the numbers in x

library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}