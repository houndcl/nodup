library(imager)
library(imagerExtra)
library(ggplot2)
library(ggthemes)
library(parallel)
library(dplyr)
library(viridis)

file.name <- 'Figure 1a.png'
img <- load.image(file.name)

# hand-picked coordinates
# already optimized
x1 <- 684
y1 <- 18
x2 <- 952 
y2 <- 288
width <- 224 
height <- 338

# optimize the coordinates to maximize the overlap
optFunc <- function(x1, y1, x2, y2, width, height, dx1, dy1, dx2, dy2) {
    x1 <- x1 + dx1
    y1 <- y1 + dy1
    x2 <- x2 + dx2
    y2 <- y2 + dy2
    img.1 <- as.cimg(img[x1 : (x1 + width), y1 : (y1 + height), 1, 1:3])
    img.2 <- as.cimg(img[x2 : (x2 + width), y2 : (y2 + height), 1, 1:3])

    RR <- cor(R(img.1), R(img.2))^2
    RG <- cor(G(img.1), G(img.2))^2
    RB <- cor(B(img.1), B(img.2))^2

    return(mean(c(RR, RG, RB)))
}

gridSearchWindow <- -5:5

grids <- expand.grid(dx1 = gridSearchWindow, dy1 = gridSearchWindow, dx2 = gridSearchWindow, dy2 = gridSearchWindow)
grids$R2 <- unlist(mclapply(1:nrow(grids), function(i) {
           optFunc(x1, y1, x2, y2, width, height, grids[i,1], grids[i,2], grids[i,3], grids[i,4])
}, mc.cores = detectCores()))

grids <- arrange(grids, desc(R2))
x1 <- x1 + grids[1, 1]
y1 <- y1 + grids[1, 2]
x2 <- x2 + grids[1, 3]
y2 <- y2 + grids[1, 4]
img.1 <- as.cimg(img[x1 : (x1 + width), y1 : (y1 + height), 1, 1:3])
img.2 <- as.cimg(img[x2 : (x2 + width), y2 : (y2 + height), 1, 1:3])


# generate the moving focal plot
## the assumption is that we should see a hug spike of R2 if we explore the suboptimal 
## conditions by shifting coordinates by several pixels
focal <- expand.grid(dx1 = 0, dy1 = 0, dx2 = -10:10, dy2 = -10:10)
focal$R2 <- unlist(mclapply(1:nrow(focal), function(i) {
           optFunc(x1, y1, x2, y2, width, height, focal[i,1], focal[i,2], focal[i,3], focal[i,4])
}, mc.cores = detectCores()))


focalPlot <- 
    ggplot(focal, aes(x=dx2, y=dy2)) + 
    geom_tile(aes(fill=R2), color="white", size=0.1) + 
    scale_fill_viridis(name="R2") + 
    theme_tufte(base_family="Helvetica") + 
    labs(x = "delta x (pixel)", y = "delta y (pixel)", title = "Focal plot") +
    theme(plot.title = element_text(hjust = 0))
ggsave("moving focal plot.pdf", plot = focalPlot, width = 5, height = 4)


# generate image & edge comparison
pdf("Image and Edge comparison.pdf", width = 6, height = 12)
layout(matrix(1:6, nrow=3, byrow = FALSE))
cor <- cor(R(img.1), R(img.2))^2
img.1.edge <- cannyEdges((img.1), alpha = 0.5, sigma = 2)
img.2.edge <- cannyEdges((img.2), alpha = 0.5, sigma = 2)
plot(img.1, main = 'Image1')
plot(img.2, main = 'Image2')
plot(img.1 - img.2, main = 'delta RGB')

plot(-img.1.edge, main = 'Image1 Edges')
plot(-img.2.edge, main = 'Image2 Edges')
plot(img.1.edge - img.2.edge, main = 'delta Edge')
dev.off()


# generate channel correlation plot
pdf("Channel correlation plot.pdf", width = 9, height = 4)
layout(t(1:3))

## R channel
plot(c(R(img.1)), c(R(img.2)), cex = 0.01, 
     xlim = c(0,1), ylim = c(0,1), xlab='Image1', ylab = 'Image2', 
     main = sprintf("R channel R2 = %.2f", cor(R(img.1), R(img.2))^2))
abline(0,1, col='red')

## G channel
plot(c(G(img.1)), c(G(img.2)), cex = 0.01,
     xlim = c(0,1), ylim = c(0,1), xlab='Image1', ylab = 'Image2',
     main = sprintf("G channel R2 = %.2f", cor(G(img.1), G(img.2))^2))
abline(0,1, col='red')

## B channel
plot(c(B(img.1)), c(B(img.2)), cex = 0.01,
     xlim = c(0,1), ylim = c(0,1), xlab='Image1', ylab = 'Image2',
     main = sprintf("B channel R2 = %.2f", cor(B(img.1), B(img.2))^2))
abline(0,1, col='red')
dev.off()


# exploratory SPE correction. not used
#    # for simplicity, grayscale the image and find the best lambda for SPE
#    optFunc.SPE <- function(x, img1, img2, func) {
#        l1 <- x[1]
#        l2 <- x[2]
#        img1 <- SPE(func(img1), l1)
#        img2 <- SPE(func(img2), l2)
#        return(1 -  cor(img1, img2)^2)
#    }
#    
#    Rlambda <- optim(par = c(0.001, 0.001), fn = optFunc.SPE, method = "L-BFGS-B", 
#                     lower = 0.0001, upper = 0.01, 
#                     img1 = img.1, img2 = img.2, func = R)
#    
#    Glambda <- optim(par = c(0.001, 0.001), fn = optFunc.SPE, method = "L-BFGS-B",
#                     lower = 0.0001, upper = 0.01,
#                     img1 = img.1, img2 = img.2, func = G)
#    
#    Blambda <- optim(par = c(0.001, 0.001), fn = optFunc.SPE, method = "L-BFGS-B",
#                     lower = 0.00005, upper = 0.05,
#                     img1 = img.1, img2 = img.2, func = B)
#    
#    Rl1 <- Rlambda$par[1]
#    Rl2 <- Rlambda$par[2]
#    Gl1 <- Glambda$par[1]
#    Gl2 <- Glambda$par[2]
#    Bl1 <- Blambda$par[1]  # blue channel doesn't converge
#    Bl2 <- Blambda$par[2]  # blue channel doesn't converge
#    
#    im1 <- img.1
#    im2 <- img.2
#    R(im1) <- SPE(R(img.1), Rl1); G(im1) <- SPE(G(img.1), Gl1); B(im1) <- 0
#    R(im2) <- SPE(R(img.2), Rl2); G(im2) <- SPE(G(img.2), Gl2); B(im2) <- 0
#    gimg.2 <- SPE(grayscale(img.2), l2)
