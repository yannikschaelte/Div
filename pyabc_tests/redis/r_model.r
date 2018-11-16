model <- function(p){
    x <- p$mean + rnorm(10000)
}

sumstat <- function(x){
    list(y0=x)
}

distance <- function(x, y){
    sum((x$y0-y$y0)^2)
}

y_obs <- list(y0=rep(0, 10000))
