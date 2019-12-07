Class7
================
Philip Dai Le
10/22/2019

Lecture 7 - 10/22/19 \#Key commads: Ctrl+Alt+i = new code chunk
Alt+Enter = run code chunk

How to write normal functions

``` r
rescale <- function(x, na.rm=TRUE, plot=FALSE, ...) {
 rng <-range(x, na.rm=na.rm)

 answer <- (x - rng[1]) / (rng[2] - rng[1])
 if(plot) {
 plot(answer, ...)
 }
 return(answer)
}
```

## R functions revisted

Source my function from last day

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
#rescale(1:10)
rescale(c(1,10,5,NA,6), na.rm=TRUE)
```

    ## [1] 0.0000000 1.0000000 0.4444444        NA 0.5555556

Warning and Stop functions -To control unexpected occurences/scenarios
in a function -Warning() gives a warning -Stop() stops/breaks the
function

``` r
rescale2 <- function(x, na.rm=TRUE, plot=FALSE, ...) {
 if( !is.numeric(x) ) {
 stop("Input x should be numeric", call.=FALSE)
 }
 rng <-range(x, na.rm=na.rm)

 answer <- (x - rng[1]) / (rng[2] - rng[1])
 if(plot) {
 plot(answer, ...)
 }
 return(answer)
}
```

Error message useful to understand desired input and source of error

``` r
#rescale2(1:10, "barry")
```

\#\#A new function called both NA Write a function to find where there
are NA in two input vetors. Firstly, make simple input where answer is
known

``` r
# Lets define an example x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

Looked online and found the **is.na() function**

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
#tells which elements are true in regards to the function is.na()
which (is.na(x))
```

    ## [1] 3 5

``` r
which (is.na(y))
```

    ## [1] 1 3

The AND function requires two input TRUE to give. Taking the **sum()**
of TRUE FALSE vectore will tell me how many TRUE elements I have. This
is my working snippet of code.

``` r
#practice with &
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
#TRUE treated as 1 and FALSE treated as 0
sum(c(TRUE,TRUE,FALSE,TRUE))
```

    ## [1] 3

Now turning the working snippet into a function

``` r
#x <- c( 1, 2, NA, 3, NA)
#y <- c(NA, 3, NA, 3, 4)

both_na<-function(x,y){
sum(is.na(x) & is.na(y))
}
both_na(x,y)
```

    ## [1] 1

Testing various outputs in a function

``` r
#eejit proofing
x <- c(NA, NA, NA)
x2<- c(NA, NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
both_na(x, y1)
```

    ## [1] 2

``` r
both_na(x2, y2)
```

    ## [1] 3

``` r
#longer object length is not a multiple of shorter object [1] 3. Error due to recycling in vectors. X is not as long as Y so it reuses the first variable of X for the 4th variable of Y
```

``` r
x<- c(1, NA, NA)
y3<- c (1, NA, NA, NA, NA, NA, NA)
both_na (x,y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
#Useful for understanding plots. The output [1] 4 is the ouput of x<-c(1, NA, NA, 1, NA, NA, 1) when using x<-c(1, NA, NA)
```

Revisiting recycling in vector

``` r
#Notice how the color variables are recycled continuously after green is used.
plot(1:10, col =c( "red", "blue", "green" ))
```

![](Class7_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
#use of function length()
length(x)
```

    ## [1] 3

``` r
length( y3)
```

    ## [1] 7

``` r
length(x) == length (y3)
```

    ## [1] FALSE

``` r
both_na2<-function(x,y){
        if(length(x)!=length(y)) {
                stop("STOP, input lengths are not equal")
        }
sum(is.na(x) & is.na(y))

}
```

Writing a function practice on student grade

``` r
# student 1, out of 800 this totals out to 790
Stu1<- c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2, out of 800 this totals out to 637 if NA = 0
Stu2<- c(100, NA, 90, 90, 90, 90, 97, 80)
#if statement has to be first, stop() will stop code at missing input
grade<-function(x){
        if(any(is.na(x))){
                warning("Student is missing homework")
        }

        finalgrade<- (sum(x, na.rm=TRUE)/800)*100
        print("Your final grade percent is")
        return(finalgrade)
}
grade(Stu1)
```

    ## [1] "Your final grade percent is"

    ## [1] 98.75

``` r
grade(Stu2)
```

    ## Warning in grade(Stu2): Student is missing homework

    ## [1] "Your final grade percent is"

    ## [1] 79.625

Dr. Barry Grant’s code solution

``` r
# student 1, out of 800 this totals out to 790
s1<- c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2, out of 800 this totals out to 637 if NA = 0
s2<- c(100, NA, 90, 90, 90, 90, 97, 80)

#basic function body: s2[-which.min(s2)]
#mean (s2[-which.min(s2)],na.rm=TRUE)
GrantGrade<-function(x){
        NiceGrade<-mean(x[-which.min(x)], na.rm=TRUE)
        #print("Your adjusted grade percent is")
        return(NiceGrade)
}
GrantGrade(s1)
```

    ## [1] 100

``` r
GrantGrade(s2)
```

    ## [1] 92.83333

``` r
#"-" cause to print out everything but minimal 
```

We have our working code now turn it into a first function

``` r
#GrantGrade<-function(x){
 #       if (any(is.na(x)) ) {
  #              warning ("Student is missing a homework")
   #     }
    #    NiceGrade<-mean(x[-which.min(x)], na.rm=TRUE)
        #print("Your adjusted grade percent is")
     #   return(NiceGrade)
#}
#GrantGrade(s1)
#GrantGrade(s2)

#s3<- c(100, NA, NA, NA, NA)
#grade2<
#grade2(s3)
```

URL example for function

``` r
url <-"https://tinyurl.com/gradeinput"
hw<- read.csv(url, row.names=1)
apply(hw,1, GrantGrade)
```

    ##  student-1  student-2  student-3  student-4  student-5  student-6  student-7 
    ##   91.75000   82.50000   84.25000   88.00000   88.25000   89.00000   94.00000 
    ##  student-8  student-9 student-10 student-11 student-12 student-13 student-14 
    ##   93.75000   87.75000   81.33333   86.00000   91.75000   92.25000   87.75000 
    ## student-15 student-16 student-17 student-18 student-19 student-20 
    ##   83.33333   89.50000   88.00000   97.00000   82.75000   82.75000

``` r
#install.packages("msa")
library (msa)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit
