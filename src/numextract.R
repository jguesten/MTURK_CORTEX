#### NUMEXTRACT (from http://stla.github.io/stlapblog/posts/Numextract.html) ####

### GREPS SCALARS FROM A CHARACTER VECTOR ###
# input: x, character vector
# output: y, character vector only containing the numbers in x

library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}