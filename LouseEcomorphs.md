Phthiraptera Ecomorphs: R Markdown Tutorial
================
Stanislav Kolencik
1/21/2023

## *Note (0): Have all the files in the same directory, and set up this directory before running the scripts.*

## *Note (1): Some of the code lines are commented with “#” due to the inability to run them in R markdown. Remove “#” in your script to be able to run these lines as well. Doubled # (“##”) shows a comment inside the script.*

## **Step 1: Landmark Data (set up your working directory where are all files located).**

### 1.A: Read ply files

### 1.B Create landmark data file (.nts)

### Repeat steps 1.A-B for each .ply file by changing input and output names

### Example:

``` r
#setwd(".")
#read.ply("Clayiella_F.ply")->Clayiella_F
#digit.fixed(spec = Clayiella_F, fixed = 14, index = FALSE, ptsize = 1, center = TRUE)->Clayiella_F
```

## **Step 2: Read the files.**

### 2.A: Read in classifier file

### 2.B: Read in NTS files into mydata

``` r
setwd(".")

classifier  <- read.csv("Classifier.csv", header=T, row.names = 1)
filelist <- list.files(path = ".", pattern = ".nts")
mydata <- readmulti.nts(filelist)
dim(mydata)
```

    ## [1] 14  3 88

``` r
data(mydata)
```

## **Step 3: Procrustes Analysis of landmark data and generating of pc scores.**

``` r
Y.gpa<-gpagen(mydata)

##PCA of procrustes landmark data (use PC1 and PC2 axes values for Figure 3)
plot(gm.prcomp(Y.gpa$coords))
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
PCA <- gm.prcomp(Y.gpa$coords)

##Plot PCA results and choose axis - standard is 1&2; it also shows the variation. "main" = name of the plot

# plot(PCA, main = "PCA")
# plot(PCA, main = "PCA",  axis1 = 1, axis2 =2, pch = classifier$Ecomorph)

##You can write PCA scores to csv file and then copy them to you classiefier file. After that you have to save it and read classifier again
PCscores <- PCA$x
write.csv(PCscores, file = "PCscores.csv")
pc1 <- read.csv("PCscores.csv", header=TRUE, row.names=1)
pc1 <- as.matrix(pc1)[,1]
# pc1
pc2 <- read.csv("PCscores.csv", header=TRUE, row.names=1)
pc2 <- as.matrix(pc2)[,2]
# pc2
pc3 <- read.csv("PCscores.csv", header=TRUE, row.names=1)
pc3 <- as.matrix(pc3)[,3]
# pc3

###Recover percent variation for PC1 and 2)
summary(PCA)->PCAsum
PCAsum$PC.summary[2,1]*100->PC1pct
PCAsum$PC.summary[2,2]*100->PC2pct
trunc(PC1pct,10^4)
paste("PC1:",round(PC1pct,digits = 2),"%")->PC1label
paste("PC2:",round(PC2pct,digits = 2),"%")->PC2label
##add pc scores to the columns of your classifier.csv file and saved it, read it here again!(in the step 4)
```

## **Step 4: PCA plot for head shape morphospace with graphics and add average points for each group (calculated in excel as average of their pc values –> X marks in plot).**

``` r
#change classifier to have corresponding label numbers
classifier  <- read.csv("ClassifierLabels.csv", header=T, row.names = 1)
par(mfrow=c(1,1))
classifier$Color->Color
unique(classifier$Sex)->Sexs
unique(classifier$Ecomorph)->Ecomorphs
unique(classifier$Color)->EcomorphColor
plot(classifier$pc1, classifier$pc2, pch = c(21, 22)[as.factor(classifier$Sex)], bg = as.character(classifier$Color),  main= "Phthiraptera morphospace", xlab = PC1label, ylab = PC2label)
points(0.121354352, 0.027218973, pch = 4, col = "darkgreen", cex = 1)
points(-0.014160578, -0.003089477, pch = 4, col = "goldenrod4", cex = 1)
points(0.041069224, -0.05550678, pch = 4, col = "dodgerblue3", cex = 1)
points(-0.12675011, 0.028453611, pch = 4, col = "firebrick4", cex = 1)
text(classifier$pc1, classifier$pc2, labels = row.names(classifier), cex = 0.3, pos = 1, col = "black")
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## **Step 5: Then you can use a polygon function.**

``` r
#  Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
#    hpts <- chull(x = xcoord, y = ycoord)
#    hpts <- c(hpts, hpts[1])
#    lines(xcoord[hpts], ycoord[hpts], col = lcolor)
#  }
# 
# Plot_ConvexHull(classifier[1:21,"pc1"],classifier[1:21,"pc2"],lcolor = "seagreen3")
# Plot_ConvexHull(classifier[64:88,"pc1"],classifier[64:88,"pc2"],lcolor ="firebrick3")
# Plot_ConvexHull(classifier[42:63,"pc1"],classifier[42:63,"pc2"],lcolor = "deepskyblue")
# Plot_ConvexHull(classifier[22:41,"pc1"],classifier[22:41,"pc2"],lcolor = "goldenrod2")
```

## **Step 6: Adding tree data**

``` r
tree <- read.tree("LouseTree_Binary.tre")

##prune tips
classifier <-read.csv("ClassifierCohab.csv", header=TRUE, row.names = 1)
species.to.keep <-read.csv("Classifier.csv", header=TRUE, row.names = 1)

##treedata(tree, Y.gpa$coords)

pruned.tree<-drop.tip(tree,tree$tip.label[-na.omit(match(species.to.keep[,1],tree$tip.label))])
write.tree(pruned.tree, file = "pruned_tree.phy")
mytree <- read.tree("pruned_tree.phy")
name.check(tree, species.to.keep)
name.check(tree, Y.gpa$coords)
##name change - ".nts_1" can make a problem with compatibility, so you can remove it.
dimnames(Y.gpa$coords)[[3]] = gsub(pattern = ".nts_1", replacement = "", x = dimnames(Y.gpa$coords)[[3]])
match(mytree$tip.label,dimnames(Y.gpa$coords)[[3]])

classifier[match(x = as.vector(mytree$tip.label),table = rownames(classifier)),]->classifier2
```

## **Step 7: Adding phylogenetic tree data to geomorph data frame and comparing net rates of morphological evolution for groups.**

``` r
###geomorph data frame including phy data:
#gdf <- geomorph.data.frame(Y.gpa, Ecomorph = classifier2$Ecomorph,Sex = classifier2$Sex, Alone = classifier2$Alone, phy = mytree)

gdf <- geomorph.data.frame(Y.gpa, Ecomorph = classifier2$Ecomorph, Sex = classifier2$Sex, Alone = classifier2$Alone, EcoAlone = classifier2$EcoAlone, phy = mytree)

##Compare evolutionary rates; "method = simulation" and/or "method = permutation"
##Setting up groups
#for ecomorph
gp.end <- factor(classifier2$Ecomorph)
names(gp.end)<-gdf$phy$tip.label
#SIMULATION
ER.E<-compare.evol.rates(Y.gpa$coords, phy=gdf$phy, gp=gp.end, iter = 10000, method = c("simulation"), print.progress = TRUE)
summary(ER.E)
plot(ER.E)
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
##for cohabiting status
gp.end <- factor(classifier2$Alone)
names(gp.end)<-gdf$phy$tip.label
#SIMULATION
ER<-compare.evol.rates(Y.gpa$coords, phy=gdf$phy, gp=gp.end, iter = 10000, method = c("simulation"), print.progress = TRUE)
summary(ER)
plot(ER)
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
##or for both ecomorph and cohabiting status
#gp.end <- factor(classifier2$EcoAlone)
#names(gp.end)<-gdf$phy$tip.label

##or for sex
#gp.end <- factor(classifier2$Sex)
#names(gp.end)<-gdf$phy$tip.label

#ER<-compare.evol.rates(Y.gpa$coords, phy=gdf$phy, gp=gp.end, iter = 10000, method = c("simulation"), print.progress = TRUE)
#summary(ER)
#plot(ER)
```

## **Step 8: Run a test for phylogenetic signal of the head shape.**

``` r
PS.shape<-physignal(gdf$coords, gdf$phy, iter=9999)
summary(PS.shape)
```

    ## 
    ## Call:
    ## physignal(A = gdf$coords, phy = gdf$phy, iter = 9999) 
    ## 
    ## 
    ## 
    ## Observed Phylogenetic Signal (K): 0.37486
    ## 
    ## P-value: 1e-04
    ## 
    ## Effect Size: 9.47115
    ## 
    ## Based on 10000 random permutations

## **Step 9: Phylogenetic and Phylogenetically-aligned PCA**

``` r
###phyPCA: Phylogenetic PCA with ordinary least-squares (OLS)
par(mfrow=c(1,1))
phyPCA <- gm.prcomp(Y.gpa$coords, phy = mytree)
#summary(phyPCA)

phyPCAScores <- phyPCA$x
write.csv(phyPCAScores, file = "phyPCAScores.csv")
phypc1 <- read.csv("phyPCAScores.csv", header=TRUE, row.names=1)
phypc1 <- as.matrix(phypc1)[,1]
#phypc1
phypc2 <- read.csv("phyPCAScores.csv", header=TRUE, row.names=1)
phypc2 <- as.matrix(phypc2)[,2]
#phypc2

plot(phyPCA, phylo = TRUE, phylo.par = list(tip.labels = FALSE, tip.txt.cex = 0.4, tip.txt.adj = c(-0.2, -0.2), edge.color = "grey", edge.width = 0.8, node.labels = FALSE, node.pch = 21, node.cex = 0.4, node.bg = "black", anc.states = TRUE), pch = 21, bg = as.character(classifier$Color))
text(phypc1, phypc2, labels = classifier2$number, cex = 0.3, pos = 1, col = "black")
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
###phyPCA: Phylogenetic PCA with generalized least-squares (GLS)
par(mfrow=c(1,1))
phyPCAG <- gm.prcomp(Y.gpa$coords, phy = mytree)
#summary(phyPCAG)

phyPCAGScores <- phyPCAG$x
write.csv(phyPCAGScores, file = "phyPCAGScores.csv")
phypcG1 <- read.csv("phyPCAGScores.csv", header=TRUE, row.names=1)
phypcG1 <- as.matrix(phypcG1)[,1]
#phypcG1
phypcG2 <- read.csv("phyPCAGScores.csv", header=TRUE, row.names=1)
phypcG2 <- as.matrix(phypcG2)[,2]
#phypcG2

plot(phyPCAG, phylo = TRUE, phylo.par = list(tip.labels = FALSE, tip.txt.cex = 0.4, tip.txt.adj = c(-0.2, -0.2), edge.color = "grey", edge.width = 0.8, node.labels = FALSE, node.pch = 21, node.cex = 0.4, node.bg = "black", anc.states = TRUE), pch = 21, bg = as.character(classifier$Color))
text(phypc1, phypc2, labels = classifier2$number, cex = 0.3, pos = 1, col = "black")
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
###PaCA: Phylogenetically-aligned PCA with ordinary least-squares (OLS)
par(mfrow=c(1,1))
PaCA <- gm.prcomp(Y.gpa$coords, phy = mytree, align.to.phy = TRUE)
#summary(PaCA)

PaCAScores <- PaCA$x
write.csv(PaCAScores, file = "PaCAScores.csv")
pac1 <- read.csv("PaCAscores.csv", header=TRUE, row.names=1)
pac1 <- as.matrix(pac1)[,1]
#pac1
pac2 <- read.csv("PaCAScores.csv", header=TRUE, row.names=1)
pac2 <- as.matrix(pac2)[,2]
#pac2

plot(PaCA, phylo = TRUE, phylo.par = list(tip.labels = FALSE, tip.txt.cex = 0.4, tip.txt.adj = c(-0.2, -0.2), edge.color = "grey", edge.width = 0.8, node.labels = FALSE, node.pch = 21, node.cex = 0.4, node.bg = "black", anc.states = TRUE), pch = 21, bg = as.character(classifier$Color))
text(pac1, pac2, labels = classifier2$number, cex = 0.3, pos = 1, col = "black")
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
###PaCA with generalized least-squares (GLS)
PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = mytree, GLS = TRUE, align.to.phy = TRUE)
#summary(PaCA.gls)
plot(PaCA.gls, phylo = TRUE, phylo.par = list(tip.labels = FALSE, tip.txt.cex = 0.4, tip.txt.adj = c(-0.2, -0.2), edge.color = "grey", edge.width = 0.8, node.labels = FALSE, node.pch = 21, node.cex = 0.4, node.bg = "black", anc.states = TRUE), pch = 21, bg = as.character(classifier$Color))

PaCA.glsScores <- PaCA.gls$x
write.csv(PaCA.glsScores, file = "PaCA.glsscores.csv")
pacg1 <- read.csv("PaCA.glsscores.csv", header=TRUE, row.names=1)
pacg1 <- as.matrix(pacg1)[,1]
#pacg1
pacg2 <- read.csv("PaCA.glsscores.csv", header=TRUE, row.names=1)
pacg2 <- as.matrix(pacg2)[,2]
#pacg2
text(pacg1, pacg2, labels = classifier2$number, cex = 0.3, pos = 1, col = "black")
```

![](LouseEcomorphs_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## **Optional add-on -> Step 10: 3D plot Phylomorphospace.**

``` r
# plot(phyPCA, time.plot = TRUE, phylo.par = list(edge.color = "grey", edge.width = 1.5, tip.labels = FALSE, node.labels = F, anc.states = F), pch = 21, bg = as.character(classifier$Color))

# plot(PaCA.gls, time.plot = TRUE, phylo.par = list(tip.labels = FALSE, edge.color = "grey", edge.width = 0.8, node.labels = FALSE, node.pch = 21, node.cex = 0.4, node.bg = "black", anc.states = TRUE), pch = 21, bg = as.character(classifier$Color))
```

## **Additional substep: HTML widget for 3D plot:**

``` r
# `movie3d(spin3d(axis=c(0,0,1), rpm=4), duration = 15, dir = "./")
# widget <- (spin3d)
# rglwidget(elementId = "PaCA.gls")
# widget <- rglwidget()
# if (interactive() || in_pkgdown_example())
#   widget
# 
# 
# if (interactive() && !in_pkgdown_example()) {
#   # Save it to a file.  This requires pandoc
#   filename <- tempfile(fileext = ".html")
#   htmlwidgets::saveWidget(rglwidget(), filename)
#   browseURL(filename)
# }
```
