
## @knitr env, include=FALSE, echo=FALSE, cache=FALSE
library("knitr")
opts_chunk$set(fig.align = 'center', 
               fig.show = 'hold', 
               par = TRUE,
               prompt = TRUE,
               eval = TRUE,
               stop_on_error = 1L,
               comment = NA)
options(replace.assign = TRUE, 
        width = 55)
set.seed(1)

## load the data using standard R data structures
## @knitr makedata1, tidy = FALSE
n <- 10
m <- 6
marray <- matrix(rnorm(n * m, 10, 5), ncol = m)
pmeta <- data.frame(sampleId = 1:m, 
                    condition = rep(c("WT", "MUT"), each = 3))
rownames(pmeta) <- colnames(marray) <- LETTERS[1:m]
fmeta <- data.frame(geneId = 1:n, 
                    pathway = sample(LETTERS, n, replace = TRUE))
rownames(fmeta) <- 
  rownames(marray) <- paste0("probe", 1:n)

## Store the data in a single variable. A list with 3 elements.
## @knitr makedata2, tidy= FALSE
maexp <- list(marray = marray,
              fmeta = fmeta,
              pmeta = pmeta)
rm(marray, fmeta, pmeta) ## clean up
str(maexp)


## @knitr access
maexp$pmeta
summary(maexp$marray[, "A"])
wt <- maexp$pmeta[, "condition"] == "WT"
maexp$marray["probe8", wt]
maexp[["marray"]]["probe3", !wt] ## different syntax


## @knitr bw1, dev='pdf', echo=TRUE
boxplot(maexp$marray)


## @knitr subset
x <- 1:5
y <- 1:3
marray2 <- maexp$marray[x, y]
fmeta2 <- maexp$fmeta[x, ]
pmeta2 <- maexp$pmeta[y, ]
maexp2 <- list(marray = marray2,
               fmeta = fmeta2,
               pmeta = pmeta2)
rm(marray2, fmeta2, pmeta2) ## clean up
str(maexp2)

## set up the data as a class with methods
## @knitr makeclass, tidy = FALSE
MArray <- setClass("MArray",
                   slots = c(marray = "matrix",
                             fmeta = "data.frame",
                             pmeta = "data.frame"))


## @knitr makeobject, tidy = FALSE
MArray() ## an empty object
#MArray(marray = 1:2) ## not allowed
# ma <- MArray(marray = maexp[[1]],
#              pmeta = maexp[["pmeta"]],
#              fmeta = maexp[["fmeta"]])       
ma <- MArray(marray = marray,
             pmeta = pmeta,
             fmeta = fmeta)       
rm(marray, fmeta, pmeta) ## clean up

class(ma)
ma


## @knitr accesswithat
ma@pmeta


## @knitr showmeth
show
isGeneric("show")
hasMethod("show")


## @knitr showmethod, tidy = FALSE
setMethod("show", 
          signature = "MArray", 
          definition = function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" ", nrow(object@marray), " features by ", 
                ncol(object@marray), " samples.\n", sep = "")
            invisible(NULL)
          })
ma 


## @knitr makegen # name the new generic and args for defining class-specific method
setGeneric("marray", function(object) standardGeneric("marray"))


## @knitr makegen2 # elipsis in the function allows for other methods to use more args
setGeneric("marray", function(object, ...) standardGeneric("marray"))


## @knitr makemeth, tidy = FALSE # make the method named above
## function "marray" is now retrieve the marray matrix in the MArray instance named in the arg
setMethod("marray", "MArray", 
          function(object) object@marray)
marray(ma)


## @knitr otheraccess, echo = FALSE 
# create function(object) to retrieve pmeta and fmeta from MArray instance
# .silent is a variable, don't know meaning or use
.silent <- setGeneric("pmeta", function(object) standardGeneric("pmeta"))
.silent <- setGeneric("fmeta", function(object) standardGeneric("fmeta"))
.silent <- setMethod("pmeta", "MArray", function(object) object@pmeta)
.silent <- setMethod("fmeta", "MArray", function(object) object@fmeta)


## @knitr syntaticsugar, tidy = FALSE
letters[1:3]
`[`(letters, 1:3)


## @knitr subset ma, tidy = FALSE
# build a bracket function for subsets of a MArray instance
setMethod("[", "MArray",
          function(x,i,j,drop="missing") {              
            .marray <- x@marray[i, j]
            .pmeta <- x@pmeta[j, ]
            .fmeta <- x@fmeta[i, ]
            MArray(marray = .marray,
                   fmeta = .fmeta,
                   pmeta = .pmeta)
          })
ma[1:5, 1:3]


## @knitr setval, tidy = FALSE
## setValidity function looks to verify that the data matches the meta-data structure
setValidity("MArray", function(object) {
  msg <- NULL
  valid <- TRUE
  if (nrow(marray(object)) != nrow(fmeta(object))) {
    valid <- FALSE
    msg <- c(msg, 
             "Number of data and feature meta-data rows must be identical.")
  }
  if (ncol(marray(object)) != nrow(pmeta(object))) {
    valid <- FALSE
    msg <- c(msg, 
             "Number of data rows and sample meta-data columns must be identical.")
  }
  if (!identical(rownames(marray(object)), rownames(fmeta(object)))) {
    valid <- FALSE
    msg <- c(msg, 
             "Data and feature meta-data row names must be identical.")        
  }
  if (!identical(colnames(marray(object)), rownames(pmeta(object)))) {
    valid <- FALSE
    msg <- c(msg, 
             "Data row names and sample meta-data columns names must be identical.")        
  }
  if (valid) TRUE else msg 
})
validObject(ma)


## @knitr notvalid
# the setValidity function for MArray now rejects, with explanation, a bad instance 
x <- matrix(1:12, ncol = 3)
y <- fmeta(ma)
z <- pmeta(ma)
MArray(marray = x, fmeta = y, pmeta = z)


## @knitr replacedirect
ma@marray <- 1
(broken <- ma)
broken@marray <- matrix(1:9, 3)
broken
validObject(broken)


## @knitr genreplacement, tidy = FALSE
# create a replacement method "object<-" to replace accessable slot in the instance
setGeneric("marray<-", 
           function(object, value) standardGeneric("marray<-"))


## @knitr replacement, tidy = FALSE
# define the replacement method
setMethod("marray<-", "MArray", 
          function(object, value) {
            object@marray <- value
            if (validObject(object))
              return(object)
          })


## @knitr replacement2, tidy = FALSE
tmp <- matrix(rnorm(n*m, 10, 5), ncol = m) 
marray(ma) <- tmp # replace with invalid data, throws error
colnames(tmp) <- LETTERS[1:m]
rownames(tmp) <- paste0("probe", 1:n)
head(marray(ma), n = 2)  #the old data
marray(ma) <- tmp # replace with new valid data, no error
head(marray(ma), n = 2)


## @knitr replacementex, echo = FALSE
# set up replacement for the meta-data
.silent <- setGeneric("fmeta<-", function(object, value) standardGeneric("fmeta<-"))
.silent <- setMethod("fmeta<-", "MArray", 
                     function(object, value) {
                       object@fmeta <- value
                       if (validObject(object))
                         return(object)
                     })
.silent <- setGeneric("pmeta<-", function(object, value) standardGeneric("pmeta<-"))
.silent <- setMethod("pmeta<-", "MArray", 
                     function(object, value) {
                       object@pmeta <- value
                       if (validObject(object))
                         return(object)
                     })


## @knitr replace pmeta, this replace adds a column
pmeta(ma)$sex <- rep(c("M", "F"), 3)
pmeta(ma)


## @knitr introspec
slotNames(ma)
getClass("MArray")


## @knitr introspec2
showMethods("marray")
showMethods(classes = "MArray")


## @knitr introspec3
getMethod("marray", "MArray")


## @knitr bioenv, echo=FALSE
suppressPackageStartupMessages(library("Biobase"))


## @knitr biob
library("Biobase")
getClass("ExpressionSet")


## @knitr sessioninfo, results='asis', echo=FALSE
toLatex(sessionInfo())