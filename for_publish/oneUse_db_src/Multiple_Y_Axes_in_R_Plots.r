
#par(mar = c(6, 6 , 2, 1))

# par(mar = rep(1,4));# down left up right


layout(matrix(c(1, 2)), heights=c(3, 1))

par(mar = c(0, 2, 3, 2) )

plot(table(rpois(100,5)), type = "h", col = "red", lwd=10, xaxt="n")

#plot(sin, -pi, 2*pi)

par(mar = c(2, 2, 0, 2) )

plot(table(rpois(100,5)), type = "h", col = "blue", lwd=10, las =2)
par(new=T)

plot(x <- sort(rnorm(47)), type = "s", yaxt = "n", xaxt = "n")
axis(4, col='red', ylim = c(0,5) )

