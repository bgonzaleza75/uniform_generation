library(microbenchmark)
library(ggplot2)
library(tidyverse)
library(rootSolve)

ntrials <- 50 # from 500
dim.min <- 3
dim.max <- 10 # from 100
dim.range <- dim.max-dim.min+1

#dev.off()
dev.new()

# Define the RealScalingNormalMethod function
RealScalingNormalMethod <- function(d) {
    BallNormal <- rnorm(d, 0, 1)
    norm <- sqrt(sum(BallNormal^2))
    return(BallNormal / norm)
}

########################
#                      #    
#     Gamma Method     #
#                      #
########################

#generate random Gammas from runif
Random.Gamma <- function(num, alpha = 1, beta = 1) {
    gammas <- rep(0, num)
    
    if (alpha > 1) {
        for (i in 1:num) {
            stop <- FALSE
            while (!stop) {   
                u1 <- runif(1, 0, 1)
                u2 <- runif(1, 0, 1)
                
                v <- (u1 * (alpha - (1 / (6 * alpha)))) / (u2 * (alpha - 1))
                
                if ((((2 * (u2 - 1)) / (alpha - 1)) + v + (1 / v)) <= 2) {
                    gammas[i] <- (alpha - 1)*v
                    stop <- TRUE
                } else if ((((2 * log(u2)) / (alpha - 1)) - log(v) + v) <=  1) {
                    gammas[i] <- (alpha - 1)*v
                    stop <- TRUE
                }
            }
        }
    } else if (alpha == 1) {
        for (i in 1:num) {    
            gammas[i] <- -1*log(runif(1, 0, 1))
        }
    } else {
        t <- 0.07 + 0.75(sqrt(1 - alpha))
        b <- 1 + ((exp(-t) * alpha) / t)
        
        for(i in 1:num) {
            stop <- FALSE
            while (!stop) {
                u1 <- runif(1, 0, 1)
                u2 <- runif(1, 0, 1)
                
                v  <- b * u1
                
                if (v <= 1) {
                    x <- t * v^(1 / alpha)
                    
                    if (u2 <= ((2 - x) / (2 + x))) {
                        gammas[i] <- x
                        stop <- TRUE
                    } else if (u2 <= exp(-x)) {
                        gammas[i] <- x
                        stop <- TRUE
                    }
                } else {
                    x <- -1 * log((t * (b - v)) / alpha)
                    y <- x / t
                    
                    if ((u2 * (alpha + y*(1 - alpha))) <= 1) {
                        gammas[i] <- x
                        stop <- TRUE
                    } else if (u2 <= y^(alpha - 1)) {
                        gammas[i] <- x
                        stop <- TRUE
                    }
                }
            }
        }
    }
    
    return(gammas / beta)
}

# Generate Bd(1) using Gamma Transformation Method
GammaMethod <- function(d) {
    r <- 1
    Ball <- rep(0, d)
    for (i in 1:(d - 2)) {
        dim <- d - i + 1
        XY <- Random.Gamma(2, (dim - 1)/2, 1)
        xi.star <- (XY[1] - XY[2])/(XY[1] + XY[2])
        Ball[i] = r * xi.star
        r <- sqrt(r^2 - Ball[i]^2)
    }
    theta <- 2 * pi * runif(1, 0, 1)
    Ball[d - 1] <- r * cos(theta)
    Ball[d] <- r * sin(theta)
    return(Ball)
}

#######################
#                     #
#     Beta Method     #
#                     #
#######################

rand.t <- function(df) {
    u1 <- runif(1, 0, 1)
    u2 <- runif(1, 0, 1)
    t <- sqrt(df*(u1^(-2/df) - 1)) * cos(2*pi*u2)
}

# generate random Betas from runif using a trick since the generation always uses alpha = beta
Random.Symmetric.Beta <- function(num, alpha = 1) {
    betas <- rep(0, num)
    
    for (i in 1:num) {
        t <- rand.t(2*alpha)
        betas[i] <- 0.5 * (1 + (t / sqrt(2*alpha + t^2)))
    }
    
    return(betas)
}

# Generate Bd(1) using the Beta Transformation Method
BetaMethod <- function(d) {
    r <- 1
    Ball <- rep(0, d)
    for (i in 1:(d - 2)) {
        dim <- d - i + 1
        xi.star <- 2*Random.Symmetric.Beta(1, (dim - 1)/2) - 1 
        Ball[i] = r * xi.star
        r <- sqrt(r^2 - Ball[i]^2)
    }
    theta <- 2 * pi * runif(1, 0, 1)
    Ball[d - 1] <- r * cos(theta)
    Ball[d] <- r * sin(theta)
    return(Ball)
}

##################################################
#                                                #
#     Scaling of Normal Vectors Method by        #
#      Generating Normal with Box-Muller         #
#                                                #
##################################################

# generate random Normals from runif
Random.Normal <- function(num, mu = 0, sigma = 1) {
    d <- ceiling(num/2)
    normals <- rep(0, 2*d)
    for (i in 1:d) {
        u1 <- runif(1, 0, 1)
        u2 <- runif(1, 0, 1)
        normals[i] <- sqrt(-2*log(u1))*cos(2*pi*u2)
        normals[d + i] <- sqrt(-2*log(u1))*sin(2*pi*u2)
    }
    normals <- sigma*normals + mu
    return(normals[1:num])
}

# Generate Bd(1) using the Scaling of Normal Vectors Method
ScalingNormalMethodBox <- function(d) {
    BallNormal <- Random.Normal(d,0,1)
    norm <- sqrt(sum(BallNormal^2))
    return(BallNormal / norm)
}

##################################################
#                                                #
#     Scaling of Normal Vectors Method by        #
#     Generating Normal with Ziggurat Method     #
#                                                #
##################################################

func <- function(x) exp(-x*x / 2)

func.inv <- function(y) sqrt(-2*log(y))

try_r_value <- function(r){
    v <- r*func(r) + exp(-0.5*r*r)/r;
    x <- rep(0, 256)
    x[256] <- r
    x[1] <- 0
    for (i in 256:3) {
        x[i - 1] <- func.inv((v / x[i]) + func(x[i]))
    }
    return(x[2]*(1 - func(x[2]))/v)
}

find.r <- function() {
    a <- 0
    b <- 10
    stop <- FALSE
    while (!stop) {
        aa <- a
        bb <- b
        r <- 0.5*(a + b)
        q = try_r_value(r)
        if ((q > 10^(-5)) & (!is.nan(q))) {
            b = r
        } else {
            a = r
        }
        if (!(aa<r && r<bb)) stop <- TRUE
    }
    return(r)
}

build.normal.tables <- function() {
    r <- find.r()
    v <- r*func(r) + (exp(-0.5*r*r) / sqrt(2*pi))/r;
    x <- rep(0, 256)
    x[256] <- r
    x[1] <- 0
    for (i in 256:3) {
        x[i - 1] <- func.inv((v / x[i]) + func(x[i]))
    }
    return(x)
}

normal.tables <- data.frame("x" = rev(build.normal.tables()))
normal.tables$y <- func(normal.tables$x)

RandNormZigg <- function(num, mu = 0, sigma = 1) {
    norms <- rep(0, num)
    iter <- 1
    
    while (iter <= num) {
        layer <- floor(255*runif(1, 0, 1)) + 1
        x <- (2*runif(1, 0, 1) - 1) * normal.tables$x[layer]
        if (abs(x) < normal.tables$x[layer + 1]) {
            norms[iter] <- x
            iter <- iter + 1
        } else if (layer == 1) {
            stop <- FALSE
            while (!stop) {
                x <- -log(runif(1, 0, 1)) / normal.tables$x[2]
                y <- -log(runif(1, 0, 1))
                if (2*y > x*x) {
                    norms[iter] <- x + normal.tables$x[2]
                    iter <- iter + 1
                    stop <- TRUE
                }
            }
        } else {
            y <- normal.tables$y[layer] + runif(1, 0, 1)*(normal.tables$y[layer + 1] - normal.tables$y[layer])
            fx <- exp(-0.5*x*x)
            if (y < fx) {
                norms[iter] <- x
                iter <- iter + 1
            }
        }
    }
    
    return(norms)
}

# Generate Bd(1) using the Scaling of Normal Vectors Method
ScalingNormalMethodZigg <- function(d) {
    BallNormal <- RandNormZigg(d,0,1)
    norm <- sqrt(sum(BallNormal^2))
    return(BallNormal / norm)
}

########################################
#                                      #
#     Spherical Coordinates Method     #
#                                      #
########################################

# Function to generate theta_i for the Spherical Coordinates function
randsin <- function(p) {
    if (round(p)!=p | p <=0){
        warning("power is not natural number.")
        return(0)
    }
    test <- FALSE
    while (!test) {
        u <- runif(1, 0, 1)
        v <- runif(1, 0, 1)
        t <- (sin(pi*v))^p
        if (u <= t) {
            theta <- pi*v
            test <- TRUE
        }
    }
    if (theta>0 && theta<pi) return(theta)
    else{
        warning("theta is not in (0,pi).")
        return(0)
    }
}

# Generate Bd(1) using the Spherical Coordinates Method
SphericalCoordinatesMethod <- function(d) {
    theta <- rep(0, (d - 1))
    for (i in 1:(d - 2)) {
        theta[i] <- randsin(d - i - 1)
    }
    theta[d - 1] <- 2*pi*runif(1, 0, 1)
    
    Ball <- rep(0, d)
    Ball[1] <- prod(sin(theta))
    for(i in 2:(d - 1)) {
        Ball[i] <- cos(theta[d - i + 1])*prod(sin(theta[1:(d - i)]))
    }
    Ball[d] <- cos(theta[1])
    return(Ball)
}

####################################
#                                  #
#     Ratio of Uniforms Method     #
#                                  #
####################################

findrects <- function() {
    rects <- rep(0, (dim.range - 1))
    for (d in 4:(dim.range+2)) {
        r <- (d - 3) / 2
        area.rect <- 2 * sqrt((r^r) / ((1 + r)^(1 + r)))
        rects[d - 3] <- area.rect
    }
    return(rects)
}

findomegas <- function() {
    omegas <- rep(0, (dim.range - 1))
    for (d in 4:(dim.range+2)) {
        r <- (d - 3) / 2
        area.omega <- (sqrt(pi) * gamma(r + 1)) / (2 * gamma(r + 1.5))
        omegas[d - 3] <- area.omega
    }
    return(omegas)
}

area.rect <- findrects()
area.omega <- findomegas()
M <- area.rect / area.omega

# density function g for Ratio of Uniforms
f.ratio <- function(x, r) {
    return((2 * x * sqrt(1 - x^(2 / r))) / area.omega[2*r])
}

# generate random xi value from the density
randRatio <- function(d) {
    if (d == 3) {
        return(2*runif(1, 0, 1) - 1)
    } else {
        r <- (d - 3) / 2
        stop <- FALSE
        while (!stop) {
            u <- runif(1, 0, 1)
            y <- runif(1, 0, 1)
            if (u <= (f.ratio(y, r) / M[d - 3])) {
                U <- y
                stop <- TRUE
            }
        }
        s <- U*sqrt(1 - U^(2 / r))
        V <- 2*s*runif(1, 0, 1) - s
        return(V / U)
    }
}

# Generate Bd(1) using the 1-Dimensional Ratio of Uniforms Method
ROUMethod <- function(d) {
    r <- 1
    Ball <- rep(0, d)
    for (i in 1:(d - 2)) {
        dim <- d - i + 1
        xi.star <- randRatio(dim)
        Ball[i] = r * xi.star
        r <- sqrt(r^2 - Ball[i]^2)
    }
    theta <- 2 * pi * runif(1, 0, 1)
    Ball[d - 1] <- r * cos(theta)
    Ball[d] <- r * sin(theta)
    return(Ball)
}



######################################
#                                    #
#     Accept/Reject Involving f*     #
#                                    #
######################################

# calculate the surface area of an n-sphere
NSphereSurfaceArea <- function(N, r = 1) (N * pi^(N / 2) * r^(N - 1)) / gamma((N / 2) + 1)

# f* density to generate from
f.star <- function(x, d) {2*(NSphereSurfaceArea(d - 1) / NSphereSurfaceArea(d)) * (1 - x^2)^((d - 3) / 2)}

f.star.inv <- function(y, d) {sqrt(1 - ((y*NSphereSurfaceArea(d)) / (2 * NSphereSurfaceArea(d - 1)))^(2 / (d - 3)))}

g.tail <- function(x, r, d) {(f.star(r, d) * (x - 1)) / (r - 1)}

try.r.value.star <- function(r, d){
    v <- r*f.star(r, d) + integrate(f.star, lower=r, upper=1, d=d)$value
    x <- rep(0, 256)
    x[256] <- r
    for (i in 256:2) {
        x[i - 1] <- f.star.inv(((v / x[i]) + f.star(x[i], d)), d)
    }
    return(x[2]*(1 - func(x[2]))/v)
}

find.r.star <- function(d) {
    a <- 0
    b <- 1
    stop <- FALSE
    while (!stop) {
        aa <- a
        bb <- b
        r <- 0.5*(a + b)
        q <- try.r.value.star(r, d)
        if (q >  10^(-5) & !is.nan(q)) {
            b = r
        } else {
            a = r
        }
        # print(paste("a =", a, "b =", b))
        if (!(aa<r && r<bb)) stop <- TRUE
    }
    return(r)
}

build.ancillary.tables <- function(d) {
    r <- find.r.star(d)
    v <- r*f.star(r, d) + integrate(f.star, lower=r, upper=1, d=d)$value
    x <- rep(0, 256)
    x[256] <- r
    x[1] <- 0
    for (i in 256:3) {
        x[i - 1] <- f.star.inv(((v / x[i]) + f.star(x[i], d)), d)
    }
    return(x)
}

x.table <- matrix(rep(0, 97*256), 256, 97)
for (i in 1:97) {
    x.table[, i] <- rev(build.ancillary.tables(i + 3))
}

y.table <- matrix(rep(0, 97*256), 256, 97)
for (i in 1:97) {
    y.table[, i] <-  f.star(x.table[, i], i + 3)
}

randZigg <- function(d) {
    if (d == 3) {
        return(2*runif(1, 0, 1) - 1)     
    } else { 
        stop <- FALSE
        sign <- sample(c(-1, 1), 1)
        
        while (!stop) {
            layer <- floor(255*runif(1, 0, 1)) + 1
            x <- runif(1, 0, 1) * x.table[layer, (d - 3)]
            if (abs(x) < x.table[layer + 1, (d - 3)]) {
                return(sign*x)
                stop <- TRUE
            } else if (layer == 1) {
                stopper <- FALSE
                while (!stopper) {
                    z <- runif(1, x.table[1], 1)
                    u <- runif(1, 0, 1)
                    test <- f.star(z, d) / g.tail(z, x.table[1], d)
                    if (u <= test) {
                        return(sign*x)
                        stop <- TRUE
                        stopper <- TRUE
                    }
                }
            } else {
                y <- y.table[layer, (d - 3)] + runif(1, 0, 1)*(y.table[layer + 1, (d - 3)] - y.table[layer, (d - 3)])
                fx <- f.star(x, d)
                if (y < fx) {
                    return(sign*x)
                    stop <- TRUE
                }
            }
        }
    }
}

ZigguratMethod <- function(d) {
    r <- 1
    Ball <- rep(0, d)
    for (i in 1:(d - 2)) {
        dim <- d - i + 1
        xi.star <- randZigg(dim) 
        Ball[i] = r * xi.star
        r <- sqrt(r^2 - Ball[i]^2)
    }
    theta <- 2 * pi * runif(1, 0, 1)
    Ball[d - 1] <- r * cos(theta)
    Ball[d] <- r * sin(theta)
    return(Ball)
}

###########################
#                         #
#     Test Uniformity     #
#                         #
###########################

# norms <- rep(0, 1000) 
# for (i in 1:1000) {
#     norms[i] <- norm(SphericalCoordinatesMethod(37), type = "2")
# }
# mean(norms)
# 
# dims <- 25
# 
# test.beta.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.beta.matrix[i,] <- BetaMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.beta.matrix[, dim]), main = paste("Beta Dim", dim))
# }
# 
# test.gamma.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.gamma.matrix[i,] <- GammaMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.gamma.matrix[, dim]), main = paste("Gamma Dim", dim))
# }
# 
# test.sphere.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.sphere.matrix[i,] <- SphericalCoordinatesMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.sphere.matrix[, dim]), main = paste("Sphere Dim", dim))
# }
# 
# test.scale.box.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.scale.box.matrix[i,] <- RealScalingNormalMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.scale.box.matrix[, dim]), main = paste("Box Normal Dim", dim))
# }
# 
# test.scale.zigg.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.scale.zigg.matrix[i,] <- RealScalingNormalMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.scale.zigg.matrix[, dim]), main = paste("Zigg Normal Dim", dim))
# }
# 
# test.rou.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.rou.matrix[i,] <- ROUMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.rou.matrix[, dim]), main = paste("RoU Dim", dim))
# }
# 
# test.zigg.matrix <- matrix(rep(0, dims*10000), 10000, dims)
# for (i in 1:10000) {
#     test.zigg.matrix[i,] <- ZigguratMethod(dims)
# }
# par(mfrow = c(5, 5))
# for (dim in 1:dims) {
#     plot(density(test.zigg.matrix[, dim]), main = paste("Ziggurat Dim", dim))
# }



###################################
#                                 #
#     Microbenchmark and Plot     #
#                                 #
###################################

# data frame to store median times of calculations
dimTime <- data.frame(dim = dim.min:dim.max,
                      beta = rep(0, dim.range),
                      gamma = rep(0, dim.range),
                      boxscaled = rep(0, dim.range),
                      spherical = rep(0, dim.range),
                      ratio = rep(0, dim.range),
                      arzigg = rep(0, dim.range))

# benchmark the calculation times and store the median time for 3 <= d <= 100
for (d in dim.min:dim.max) {
    dimTime[d - 2, -1:0] <- (summary(microbenchmark(
                            beta = BetaMethod(d),
                            gamma = GammaMethod(d),
                            boxscaled = ScalingNormalMethodBox(d),
                            spherical = SphericalCoordinatesMethod(d),
                            ratio = ROUMethod(d),
                            arzigg = ZigguratMethod(d),
                            times = ntrials))$median) / 1000
}

# convert dimTime to make plotting easier
dimTimePivoted <- pivot_longer(data = dimTime, 
                               cols = c("beta", "gamma", "boxscaled", "spherical", "ratio", "arzigg"),
                               values_to = "values",
                               names_to = "colored_var")

# Plot the different methods against each other
p <- ggplot(data = dimTimePivoted, aes(x = dim, y = values, color = colored_var)) +
    geom_line(size=1.25) + 
    geom_point() +
    theme(legend.position = "bottom") +
    labs(colour = "Method:") +
    scale_color_hue(labels = c("Beta", "Gamma", "Scaling of Normal Vectors", "Spherical Coordinates", "Ratio of Uniforms", "Accept-Reject Using f*"))
print(p)

##############################
#                            #
#     Miscellaneous Bits     #
#                            #
##############################

# summary(microbenchmark(
#     real = RealScalingNormalMethod(100),
#     #mine = RealScalingNormalMethod2(100),
#     times = 1000))
# 
# summary(microbenchmark(
#     new = ROUMethod(100),
#     real = ROUMethod.old(100),
#     times = 1000))

d <- 20
thetas <- matrix(rep(0, 10000*d), 10000, d)
for (j in 1:10000) {
    for (i in 1:(d - 1)) {
        thetas[j, i] <- randsin(d - i)
    }
    thetas[j, d] <- 2*pi*runif(1, 0, 1)
}





# QUESTIONS
# 1. how do generate a vector of independent gamma randome variable