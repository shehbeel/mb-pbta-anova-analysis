# Load MB subset from PBTA dataset
dat <- read.csv("mb_pbta_subtypes.csv", header=T, row.names=1)

# The values in dataset are of "character" class. The dataset will be converted to "numeric" class
# library(dplyr)
# dat <- mutate_all(dat, function(x) as.numeric(x))

# Remove missing values
# dat <- na.omit(dat)

### Principal Component Analysis###
# Find the principal components
dat.pca <- prcomp(t(dat))

# Look at the first four components
dat.loadings <- dat.pca$x[,1:4]

# Plot the principal component 1 and 2
plot(range(dat.loadings[,1]),range(dat.loadings[,2]),type="n",xlab='p1',ylab='p2',main='PCA plot of MB PBTA miRNA Dataset\np2 vs. p1')
points(dat.loadings[,1], dat.loadings[,2],col=1,bg='red',pch=21,cex=1.5)
text(dat.loadings,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# Plot the principal component 1 and 3
plot(range(dat.loadings[,1]),range(dat.loadings[,3]),type="n",xlab='p1',ylab='p3',main='PCA plot of MB PBTA miRNA Dataset\np3 vs. p1')
points(dat.loadings[,1], dat.loadings[,3],col=1,bg='red',pch=21,cex=1.5)
text(dat.loadings,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# Plot the principal component 1 and 4
plot(range(dat.loadings[,1]),range(dat.loadings[,4]),type="n",xlab='p1',ylab='p4',main='PCA plot of MB PBTA miRNA Dataset\np4 vs. p1')
points(dat.loadings[,1], dat.loadings[,4],col=1,bg='red',pch=21,cex=1.5)
text(dat.loadings,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# Plot the principal component 2 and 3
plot(range(dat.loadings[,2]),range(dat.loadings[,3]),type="n",xlab='p2',ylab='p3',main='PCA plot of MB PBTA miRNA Dataset\np3 vs. p1')
points(dat.loadings[,2], dat.loadings[,3],col=1,bg='red',pch=21,cex=1.5)
text(dat.loadings,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# Plot the principal component 2 and 4
plot(range(dat.loadings[,2]),range(dat.loadings[,4]),type="n",xlab='p2',ylab='p4',main='PCA plot of MB PBTA miRNA Dataset\np4 vs. p2')
points(dat.loadings[,2], dat.loadings[,4],col=1,bg='red',pch=21,cex=1.5)
text(dat.loadings,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# Plot the principal component 3 and 4
plot(range(dat.loadings[,3]),range(dat.loadings[,4]),type="n",xlab='p3',ylab='p4',main='PCA plot of MB PBTA miRNA Dataset\np4 vs. p3')
points(dat.loadings[,3], dat.loadings[,4],col=1,bg='red',pch=21,cex=1.5)
text(dat.loadings,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# Another way to look at the plots
biplot(dat.pca, main='PCA plot of MB PBTA miRNA Dataset\np2 vs. p1') # Does not load
biplot(dat.pca)

# Cluster dendrogram to further see if samples can be distinguished
dat.dist <- dist(t(dat), method = "euclidean") # Calculate distance on transposed data
dat.clust <- hclust(dat.dist, method="single") # Calculate clusters
plot(dat.clust,labels=names(dat), main="Cluster Dendrogram of PBTA miRNA Dataset", xlab="Samples", cex=0.75)	# Plot cluster tree


#### Result: The samples aren't very clustered and distinguished after dimensionality reduction :(

### Gene Filtering to better distinguish the samples
# First log2 transform the dataset
dat.log <- log2(dat)

# Remove any NA's
dat.log[is.na(dat.log) | dat.log == "Inf" | dat.log == "-Inf"] <- NA
dat.log.clean <- na.omit(dat.log)

# Store sample names
shh <- colnames(dat.log[,c(2,3,4,8,18,19,22,25)])
wnt <- colnames(dat.log[,c(14,20)])
group3 <- colnames(dat.log[,c(6,7,9,11,26)])
group4 <- colnames(dat.log[,c(1,5,10,12,13,15,16,17,21,23,24)])

# FUNCTION: 1-factor ANOVA with ten levels
aov.all.genes <- function(x,s1,s2,s3,s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])

  fac <- c(rep("S",length(x1)), rep("W",length(x2)), 
           rep("G3",length(x3)), rep("G4",length(x4))
          )
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4))
  names(a.dat) <- c("factor","express")
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}

# Run the anova function
aov.run <- apply(dat.log.clean,1, aov.all.genes,
                 s1=shh,s2=wnt,s3=group3,s4=group4
                )

# Reorder the log-transformed dataset by probe names and add the p-values from ANOVA test as column
dat.log.clean[order(row.names(dat.log.clean)), ]
# Add ANOVA p-values as column
dat.anova <- cbind(dat.log.clean, aov.run)
# Rename last column
names(dat.anova)[names(dat.anova) == "aov.run"] <- "p.val"

# Filter out for genes with p-value < 0.10
top10.indexes <- which(dat.anova["p.val"] < 0.10) # Get the index of the top 5% DE genes
dat.anova.top10 <- dat.anova[top10.indexes,]

# Filter out for genes with p-value < 0.05
top5.indexes <- which(dat.anova["p.val"] < 0.05) # Get the index of the top 5% DE genes
dat.anova.top5 <- dat.anova[top5.indexes,]


# Filter out for genes with p-value < 0.01
top1.indexes <- which(dat.anova["p.val"] < 0.01) # Get the index of the top 1% DE genes
dat.anova.top1 <- dat.anova[top1.indexes,]


# Filter out for genes based on bonferroni correction
bc.indexes <- which(dat.anova["p.val"] < (0.05/2083)) # Get the index of the top DE genes
dat.anova.bc <- dat.anova[bc.indexes,]


# Number of genes
nrow(dat.anova.top10) # 3 genes
nrow(dat.anova.top5) # 1 genes
nrow(dat.anova.top1) # 1 genes
nrow(dat.anova.bc) # 0 genes

# Perform dimensionality reduction again for p<0.10 subset dataset
dat.anova.top10.pca <- prcomp(t(dat.anova.top10[,1:26]))
# Look at the first two components
dat.anova.top10.loadings <- dat.anova.top10.pca$x[,1:2]
# PCA Plot
plot(range(dat.anova.top10.loadings[,1]),range(dat.anova.top10.loadings[,2]),type="n",xlab='p1',ylab='p2',main='PCA Plot After Gene Filtering based on \np<0.01 PC2 vs. PC1')
points(dat.anova.top10.loadings[,1], dat.anova.top10.loadings[,2],col=1,bg='red',pch=21,cex=1.5)
text(dat.anova.top10.loadings,label=dimnames(dat.anova.top10)[[2]],pos=1,cex=0.5)

# Export the dataset
write.csv(dat.anova, "/Users/arifs2/Documents/mb_pbta_anova_analysis.csv")
