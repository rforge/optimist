##
##  r a s t r i g i n . R  Rastrigin Test Function
##


##  Rastrigin function
# No. of Vars.:  n >= 1
# Bounds:  -5.12 <= xi <= 5.12
# Local minima:  many
# Minimum:  0.0
# Solution:  xi = 0, i=1:n
rastrigin <- function(x) {
    n <- length(x)
    10*n + sum(x^2 - 10*cos(2*pi*x))
}
