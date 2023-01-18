## example R code (intro to R)
## author: alice yue


# Adding "#" at the beginning of this line makes this a comment!
# Use comments to tell the reader what you are doing (documentation).
# Comments are not run.

## variables

# numerics, used to perform arithmetics
1 + 1
2 * 3
8 / 2
(2 + 3) * 5

x1 <- 5
x1

x1 <- 1 # replaces previously assignmed value, 5
x1

x2 <- 2
x3 <- x1 + x2


# strings, a sequence of characters
s1 <- "one plus one"
s2 <- "is"
s3 <- paste0(s1, " ", s2, ": ", x3)

# print a message to console
print(s3)


# vectors, lists; indexing, length
xv1 <- c(x1, x2, x3, 5)
xv2 <- c(1:5)
sv1 <- c(s1, s2, "last string")
sv2 <- rep("a", 4) # rep = repeat

xv1
xv1[1]
length(xv1)

xv2
xv2[2]
length(xv2)

sv1
sv1[c(2:3)]
length(sv1)

sv2
sv2[length(sv2)]
length(sv2)

# lists
l1 <- list(xv1, xv2, sv1, sv2)
l1
length(l1)

names(l1) <- c("a", "b", "c", "d")
l1[[1]]
l1[c(1:3)]
l1[["c"]]


## matrices; dimensions, column/row naming, indexing
sm <- matrix("a", nrow=3, ncol=2)
sm
ncol(sm)
nrow(sm)
dim(sm) # dim=dimensions

xm <- matrix(0, nrow=2, ncol=3)
xm

nm <- matrix(c(1:9), nrow=3, ncol=3)
nm
nmcolumns <- paste0("column", c(1:ncol(nm)))
nmrows <- paste0("row", c(1:nrow(nm)))
colnames(nm) <- nmcolumns
rownames(nm) <- nmrows
nm

nm[2,1] # single numeric
nm[c(1:2),2] # matrix
nm[1,] # vector
nm[,2]
nm["row1",]
nm[,"column2"]

nm2 <- nm[,2] # vector
length(nm2)

nm2 <- nm[,2, drop=FALSE] # matrix
dim(nm2)

newm <- rbind(xm, nm) # row bind two matrices together; column bind with cbind()
newm
dim(newm)

# data.frame() is like matrix but each column 
# can contain different types of variables
d <- data.frame(name=c("bob", "mike", "emily"), age=c(10, 12, 13))
d
dim(d)
nrow(d)
ncol(d)
d[["name"]] # columns can be accessed like a list element
d[["age"]]


## if condition
x1 <- 10

if (x1 < 10) {
    print("x1 is less than 10")
}

if (x1 <= 10) {
    print("x1 is less than or equal to 10")
}
if (x1 != 10) {
    print("x1 is not 10")
}
if (x1 == 10) {
    print("x1 is 10")
}
if (x1 >= 10) {
    print("x1 is greater than or equal to 10")
}
if (x1 > 10) {
    print("x1 is greater than 10")
}

if (x1 < 10) {
    print("x1 is less than 10")
} else {
    print("x1 is not less than 10")
}

if (x1 < 10) {
    print("x1 is less than 10")
} else if (x1 == 10) {
    print("x1 is 10")
} else {
    print("x1 is greater than 10")
}

if (x1 > 5 & x1 < 15) {
    print("x1 is greater than 5 and less than 15")
}
if (x1 > 25 | x1 < 15) {
    print("x1 is either greater than 25 or less than 15")
}


## loops: for, while
for (i in c(3:10)) {
    print(i)
}

i <- 3
while (i <= 10) {
    print(i)
    i <- i + 1
}


## functions, a set of pre-written code
## to see a description of the function, type "?" in front of it
?print
?length
?matrix
?dim
?ncol
?nrow
?colnames
?rownames

# why do we need functions?
# imagine we need to write a piece of code over and over again:
p1 <- 2
p2 <- 4
p3 <- p1 * p2
print(p3)

p1 <- 3
p2 <- 5
p3 <- p1 * p2
print(p3)

p1 <- 6
p2 <- 7
p3 <- p1 * p2
print(p3)

# you can create your own function to make this code shorter!
myfunction <- function(p1, p2) {
    p3 <- p1 * p2
    return(p3)
}

# same code but shorter
p3 <- myfunction(2, 4)
print(p3)

p3 <- myfunction(3, 5)
print(p3)

p3 <- myfunction(6, 7)
print(p3)

# functions are encapsulated: 
# everything that happens in a function stays in the function.
p1 <- 3
p2 <- 5
p3 <- 7
r1 <- myfunction(6, 8)

p1
p2
p3
r1

# if you write many functions, you can bundle these functions into a "package"
# and share it with the community :)


## TRY: practice problems ####

# 1. write a function called "function1" that does the following when called:

function1(2, 3, 6)
# output printed (print()) to console:
# [1] 2
# [1] 6
# [1] 18
# [1] 54
# [1] 162
# [1] 486

# 2. write a function called "function2" that, when given two vectors,
#    return() TRUE if the two vectors are identicle, and FALSE otherwise:

function2(c(1,2,3), c(2,3,4))
# [1] FALSE
function2(c(1,1,2), c(1,1))
# [1] FALSE
function2(c(2,3), c(2,3))
# [1] TRUE

# 3. write a function called "function3" that rbind() all the matrices
#    that have the same number of columns in a list given by the user.
#    return() a list that contains all the rbind()-ed matrices.
#    hint: try using the function unique(), (use ?unique to read about it!)

list3 <- list(
    matrix(0, nrow=2, ncol=4),
    matrix(1, nrow=3, ncol=5),
    matrix(2, nrow=4, ncol=4),
    matrix(3, nrow=2, ncol=5),
    matrix(4, nrow=5, ncol=3)
)
function3(list3)
# returns a list with 3 elements:
# 
# [[1]]
# [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    0
# [3,]    2    2    2    2
# [4,]    2    2    2    2
# [5,]    2    2    2    2
# [6,]    2    2    2    2
# 
# [[2]]
# [,1] [,2] [,3] [,4] [,5]
# [1,]    1    1    1    1    1
# [2,]    1    1    1    1    1
# [3,]    1    1    1    1    1
# [4,]    3    3    3    3    3
# [5,]    3    3    3    3    3
# 
# [[3]]
# [,1] [,2] [,3]
# [1,]    4    4    4
# [2,]    4    4    4
# [3,]    4    4    4
# [4,]    4    4    4
# [5,]    4    4    4


## practice problem solutions ####

# 1.
function1 <- function(a, b, c) {
    for (i in c(1:c)) {
        print(a)
        a <- a*b
    }
}


# 2.
function2 <- function(a, b) {
    if (length(a) != length(b)) {
        return(FALSE)
    }
    for (i in c(1:length(a))) {
        if (a[i] != b[i]) {
            return(FALSE)
        }
    }
    return(TRUE)
}

function2 <- function(x, y) {
    if (length(x) == length(y)) {
        # two vectors have same length
        i <- 1
        while (i <= length(x)) {
            if (x[i] != y[i]) {
                return(FALSE)
            }
            i <- i + 1
        }
        return(TRUE)
    } else {
        return(FALSE)
    }
}


# 3.
function3 <- function(list3) {
    ncol_n <- rep(0, length(list3))
    for (l3i in c(1:length(list3))) {
        ncol_n[l3i] <- ncol(list3[[l3i]])
    }
    ncol_u <- unique(ncol_n)
    listr <- list()
    for (nui in c(1:length(ncol_u))) {
        l3d <- list3[ncol_n==ncol_u[nui]]
        l3m <- l3d[[1]]
        if (length(l3d) > 1) {
            for (l3m_ in l3d[-1]) {
                l3m <- rbind(l3m, l3m_)
            }
        }
        listr[[nui]] <- l3m
    }
    return(listr)
}

funtion3 <- function(l3) {
    # rbind() all matrices in l3 with the same number of columns
    if (length(l3) <= 1) {
        return(l3)
    }
    # sapply(l3, ncol)
    
    # find the number of columns for each matrix
    ncolumns <- rep(0, length(l3))
    for (i in c(1:length(l3))) {
        ncolumns[i] <- ncol(l3[[i]])
    }
    uc <- unique(ncolumns)
    l3_new <- l3
    l3_new <- list()
    for (i in uc) {
        a <- matrix(0, nrow=0, ncol=0)
        for (j in c(1:l3)) {
            if (i==ncolumns[j]) {
                if (ncol(a)!=i) {
                    a <- l3[[j]]
                } else {
                    a <- rbind(a, l3[[j]])
                }
                
            }
        }
        l3_new <- append(l3_new, a)
    }
}