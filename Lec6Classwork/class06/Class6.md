Class06 R functions
================
Philip Dai Le
10/17/2019

This is my work flow for class 06 in **BIMM 143**: “Intro to
bioinformatics”.

``` r
#this is to demo a code chunk. This is a code chunk in R
plot(1:10)
```

![](Class6_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

\#\#Practice reading files (again…) read.table() functions differ in
parameters: read.delim; uses sep=“ read.csv; uses sep=”," read.csv2;
uses sep=“;” read.table; uses sep=" "

``` r
read.table("test1.txt", sep=",", header=TRUE)
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
read.table("test2.txt")
```

    ##               V1
    ## 1 Col1$Col2$Col3
    ## 2          1$2$3
    ## 3          4$5$6
    ## 4          7$8$9
    ## 5          a$b$c

``` r
read.table("test3.txt", sep="\t", header = FALSE)
```

    ##           V1
    ## 1 1   6   a 
    ## 2  2   7   b
    ## 3 3   8   c 
    ## 4  4   9   d
    ## 5  5   10  e

``` r
read.csv 
```

    ## function (file, header = TRUE, sep = ",", quote = "\"", dec = ".", 
    ##     fill = TRUE, comment.char = "", ...) 
    ## read.table(file = file, header = header, sep = sep, quote = quote, 
    ##     dec = dec, fill = fill, comment.char = comment.char, ...)
    ## <bytecode: 0x0000000015211738>
    ## <environment: namespace:utils>

``` r
read.table("https://bioboot.github.io/bimm143_S19/class-material/test2.txt", sep="$",header= TRUE)
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
add <- function(x, y=1) {
 # Sum the input x and y
 x + y
}
```

Use ctrl+alt+I to insert code chunk

``` r
add(1)
```

    ## [1] 2

``` r
add(5,5) 
```

    ## [1] 10

``` r
#overides default y value
#add(5,"barry"), add(5,5,5) 
#can't use a third value for Z cause not initiated and can't use non-numerical value
```

“Our wee function is vectorized :-)” - Barry Grant

``` r
add(x=1, y=4)
```

    ## [1] 5

``` r
add(1, 4)
```

    ## [1] 5

``` r
add(1)
```

    ## [1] 2

``` r
add( c(1, 2, 3) )
```

    ## [1] 2 3 4

``` r
add( c(1, 2, 3), 4 )
```

    ## [1] 5 6 7

``` r
#add(1, 2, 2), error due to third value not initialized 
#add(x=1, y=“b”), error due to y="b"
```

What does this code do? df\(a <- (df\)a - min(df\(a)) / (max(df\)a) -
min(df\(a)) df\)b \<- (df\(b - min(df\)a)) / (max(df\(b) - min(df\)b))
df\(c <- (df\)c - min(df\(c)) / (max(df\)c) - min(df\(c)) df\)d \<-
(df\(d - min(df\)d)) / (max(df\(a) - min(df\)d))

Silent killers in code, trying to use “a” for df\(b and df\)d

Here the the intent is far more clear: df\(a <- rescale (df\)a) Purpose
of code is clearer, reduce mistakes from copy/paste, updating code is
easier, and reduce duplication and facilitate re-use

Start with code snippet to simplify and reduce calculatino duplication
Ex:x \<- (x - min(x)) / (max(x) - min(x))

``` r
rescale <- function(x) {
 rng <-range(x)
 (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
?range
```

    ## starting httpd help server ... done

``` r
range (1:10)
```

    ## [1]  1 10

``` r
range (1:20)
```

    ## [1]  1 20

``` r
y<-range(c(5,1,2,5,20))
y[1]
```

    ## [1] 1

``` r
y[2]
```

    ## [1] 20

``` r
c("a","b")
```

    ## [1] "a" "b"

Test some

``` r
rescale(c(1,2,NA,3,10))
```

    ## [1] NA NA NA NA NA

``` r
x<-c(1,2,NA,3,10)
rng<-range(x, na.rm = TRUE)
rng
```

    ## [1]  1 10

``` r
rescale2<- function(x) {
  rng<-range(x, na.rm = TRUE)
  (x-rng[1])/(rng[2]-rng[1])
}
```

``` r
rescale2 (c(1,2,NA,3,10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
 rng <-range(x, na.rm=na.rm)
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 return(answer)
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 return(answer)
}
```

``` r
rescale3(1:10, plot = FALSE)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale3(1:10, plot=TRUE)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Additional Info regarding header sizes \#This is H1 Creates small header

## This is H2

Creates large header

### This is H3

Creates medium header
