#Class 5 - Data visualization

#Section 1 X
x <- rnorm(1000)
mean(x)
sd(x)
summary(x)
boxplot(x)
hist(x)

#section 2 scatterplots
#read input file"weight_chart"
#parameters in read.csv() corresponds to adjustable traits
read.table("bimm143_05_rstats/weight_chart.txt")
read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

#2B
mouse<- read.table("bimm143_05_rstats/feature_counts.txt" ,sep="\t")

#barplot(features$Count,horiz=TRUE)


#barplot(features$Count,horiz=TRUE, names.arg=mouse$Feature,main="mouse", las=1)

par(mar=c(5,5,2,2))

?par

mouse<-read.delim("bimm143_05_rstats/feature_counts.txt",header = TRUE, sep="\t")
barplot(mouse$Count, names.arg = mouse$Feature, main= "Number of features in the mouse GRCm38 genome", las=1, horiz=TRUE)
par(mar=c(5,11,2,2))$mar


#section 3
maleFemale <- read.delim("bimm143_05_rstats/male_female_counts.txt")
barplot(maleFemale$Count,main="Male and Female Count", col=rainbow(10), xlab="Gender split between males and female",ylab="gender frequency")

