## example R code (intro to R part 2)
## author: alice yue


## TRUE/FALSE variables ####

x <- TRUE
y <- FALSE

if (x) {
    print("x is true")
}
if (y) {
    print("y is true")
}

# condition statements evaluate to TRUE/FALSE!
x1 <- 1
x <- x1<=10 & x1>-5 & x1!=0
x
x1 <= 10
if (x) {
    print("x1 is betwen -4 and 10 and is not 0")
}

## TRUE/FALLSE can be used for indexing vectors
xv1 <- c(1, 3, 5, 10, -1, 0, Inf, -Inf)
xv2 <- append(xv1, c(NA, NULL, log(-1))) # append() merges two vectors or lists
xv1
xv2

xv1i <- xv1 < 10
xv1i
xv1I <- which(xv1i)
xv1I
xv1[xv1i]
xv1[xv1I]

xv2i <- xv2 < 10
xv2i
which(xv2i)
xv2[xv2i] # by default, all NaN and NA are converted to NA

## TRUE/FALLSE can be used for indexing in matrices and data.frames()
m1 <- matrix(c(1:20), ncol=5, nrow=4)
m1
colnames(m1) <- c("a", "b", "c", "4", "5")
m1

m1i <- m1[,"a"] > 2
m1[m1i, "a"]
m1[which(m1i), "a"]
m1[m1i,] # it can be used here because: length(m1i) == nrow(m1) is TRUE

m1j <- m1[2,] > 2
m1[2, m1j]
m1[2, which(m1j)]
m1[ , m1j]


## plotting ####
m2 <- matrix(c(1:10), ncol=2, nrow=5)
m2
colnames(m2) <- c("x", "y")
plot(m2)
lines(m2)
abline(h=7)
abline(v=2)
abline(h=c(8,9))
abline(v=c(2,3), lty="dashed", col="red")
graphics.off() # close plot

# we can also change the size and colour of our points
plot(m2, cex=0.1)
plot(m2, cex=0.5)
plot(m2, cex=1)

plot(m2, col="red")
plot(m2, cex=2, col="red")

# we can also customize the size for each row
plot(m2, cex=c(0.1, 0.5, 1, 2, 3)) 
plot(m2, cex=c(0.1, 0.5, 1, 2, 3), col=c("red", "blue", "green", "black", "yellow")) 

# since size and colours are vectors, we can use condition indexing!
colours <- rep("red", 5)
colours
m2i <- m2[,"x"] < 2 & m2[,"y"] < 8
colours[m2i] <- "blue"
colours
plot(m2, cex=0.5, col=colours)

# add a title with main!
plot(m2, cex=0.5, col=colours, main="my first plot!")


## TRY: practice problems ####

# 4. write a function called "myfunction4" that outputs the following:
#    (you can assume the user gives 3 whole numbers as input arguments)

a <- myfunction4(2, -5, 6)
print(a)
# > a <- myfunction4(2, -5, 6)
# [1] 8
# [1] 7
# [1] 6
# [1] 5
# [1] 4
# [1] 3
# [1] 2
# [1] 1
# > print(a)
# [1] 1

a <- myfunction4(1, 1, 4)
print(a)
# > a <- myfunction4(1, 1, 4)
# [1] "first two arguments shouldn't be the same!"
# > print(a)
# [1] 4

a <- myfunction4(0, 10, 3)
print(a)
# > a <- myfunction4(0, 10, 3)
# [1] 3
# [1] 4
# [1] 5
# [1] 6
# [1] 7
# [1] 8
# [1] 9
# [1] 10
# [1] 11
# [1] 12
# [1] 13
# > print(a)
# [1] 13

# 5. write a function called "myfunction5" that takes in two numeric variables
#    a and b. it will return a vector of length 3. 
#    if a is greater than b, it will return c(TRUE, FALSE, FALSE)
#    if a is equal to b, it will return c(FALSE, TRUE, FALSE)
#    if a is less than b, it will return c(FALSE, FALSE, TRUE)
myfunction5(4,3)
# TRUE FALSE FALSE
myfunction5(5,5)
# FALSE TRUE FALSE
myfunction5(1,3)
# FALSE FALSE TRUE

# 6. write a function called "myfunction6" that print() out the following given
#    three numeric variables.
myfunction6(1,3,50)
# [1] "x is less than y"
# [1] "x is less than z"
# [1] "y is less than z"
myfunction6(-9,0,-10)
# [1] "x is less than y"
# [1] "x is greater than z"
# [1] "y is greater than z"
myfunction6(2,100,3)
# [1] "x is less than y"
# [1] "x is less than z"
# [1] "y is greater than z"


## practice problem solutions
# 4. 
myfunction4 <- function(x, y, z) {
    if (x==y) {
        print("first two arguments shouldn't be the same!")
        return(z)
    }
    for (i in c(x:y)) {
        a <- z + i
        print(a)
    }
    return(a)
}

# 5.
myfunction5 <- function(a, b) {
    bools <- c(a < b, a == b, a > b)
    return(bools)
}

# 6.
myfunction6 <- function(x, y, z) {
    xy <- myfunction5(x, y)
    xz <- myfunction5(x, z)
    yz <- myfunction5(y, z)
    answer <- c("less than", "equal to", "greater than")
    print(paste0("x is ", answer[xy], " y"))
    print(paste0("x is ", answer[xz], " z"))
    print(paste0("y is ", answer[yz], " z"))
}
