#first download bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library(caret)
#exercise 1.
#erwthma 1
y = NULL
p_values= NULL
true = NULL
r=NULL
for (i in 1:1000) {
n=500
p=100
x = matrix(rnorm(n*p),nrow = n, ncol = p)


b = numeric(p)
if( runif(1) < 0.3) {b[1] = rnorm(1)
                    true = rbind(true,1)
  }       else {
  true = rbind(true,0)
}
table(true)

y = x%*%b + rnorm(n)
y

model = lm(y~x)
m =summary(model)
p_values[i] = pf(m$fstatistic[1],m$fstatistic[2],m$fstatistic[3],lower.tail=FALSE)
}

methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")

adjusted_p_values <- p.adjust(p_values, method = 'hochberg')
r =cbind(adjusted_p_values <=0.05,true)
test= r[,1]
truth= r[,2]
conf = confusionMatrix(factor(test),factor(truth))
power = as.numeric(conf$byClass[3])
power

adjusted_p_values <- p.adjust(p_values, method = 'BY')
r =cbind(adjusted_p_values <=0.05,true)
test= r[,1]
truth= r[,2]
conf = confusionMatrix(factor(test),factor(truth))
power = as.numeric(conf$byClass[3])
power




adjusted_p_values <- p.adjust(p_values, method = 'BH')
r =cbind(adjusted_p_values <=0.05,true)
test= r[,1]
truth= r[,2]
conf = confusionMatrix(factor(test),factor(truth))
power = as.numeric(conf$byClass[3])
power




adjusted_p_values <- p.adjust(p_values, method = 'bonferroni')
r =cbind(adjusted_p_values <=0.05,true)
test= r[,1]
truth= r[,2]
conf = confusionMatrix(factor(test),factor(truth))
power = as.numeric(conf$byClass[3])
power


adjusted_p_values <- p.adjust(p_values, method = 'hommel')
r =cbind(adjusted_p_values <=0.05,true)
test= r[,1]
truth= r[,2]
conf = confusionMatrix(factor(test),factor(truth))
power = as.numeric(conf$byClass[3])
power



adjusted_p_values <- p.adjust(p_values, method = 'holm')
r =cbind(adjusted_p_values <=0.05,true)
test= r[,1]
truth= r[,2]
conf = confusionMatrix(factor(test),factor(truth))
power = as.numeric(conf$byClass[3])
power








#erwthma 2
power = NULL
alphas <- seq(0.01, 1, by = 0.01)
for (i in alphas){
  adjusted_p_values <- p.adjust(p_values, method = 'hochberg')
  r =cbind(adjusted_p_values <=i,true)
  test= r[,1]
  truth= r[,2]
  dat =table(truth,test)
  conf = confusionMatrix(factor(test),factor(truth))
  power = c(power,as.numeric(conf$byClass[3]))

}


plot(alphas,power,type= 'l')
round(p_values,2)
r =cbind(p_values <=0.05,true)
sum(r[,1] == r[,2])

###erwthma 3
cbind(p_values<=0.05,true)
confusionMatrix((p_values<=0.05),factor(true))
sum(p_values>0.05)/1000
hist(p_values,prob = TRUE)
abline(h =0.721,col='red')



adjusted_p_values <- p.adjust(p_values, method = 'hochberg')
sum(adjusted_p_values <= 0.05)
a = r[,1]
b = r[,2]
sum(b)
table(b,a)
adjusted = cbind(adjusted_p_values <=0.05,true)
test= adjusted[,1]
truth= adjusted[,2]
table(truth,test)

#power
741/800
powers[i] <- tp / (tp + fn)

alphas <- seq(0.01, 1, by = 0.01)

#pvalue adjusted
methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")
qvalue



BiocManager::install("qvalue")
library(qvalue)
qvalue(p=p_values)

library('leukemiasEset')
library('limma')
library('qvalue')
install.packages('VennDiagram')
library('VennDiagram')
library('plotrix')

data(leukemiasEset)
leukemiasEset

table(leukemiasEset$LeukemiaType)

ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
ourData$LeukemiaType <- factor(ourData$LeukemiaType)


# names of variables
sampleNames(ourData)
# get expression data only
x <- exprs(ourData)

pairs(x[,c(1,2,13,14)], panel=function(x,y){smoothScatter(x,y,add=T)})

design <- model.matrix(~ ourData$LeukemiaType - 1)
colnames(design) <- c("ALL", "NoL")
contrast.matrix <- makeContrasts("ALL-NoL", levels = design)
fit1 <- lmFit(ourData, design)
fit2 <- contrasts.fit(fit1, contrast.matrix)
fit3 <- eBayes(fit2)
res <- topTable(fit3, number = Inf, sort = "none")
hist(res$P.Value)

qStar <- 0.01
limma_genes <- res$adj.P.Val < qStar



n <- dim(x)[1]
pval <- numeric(n)
for(i in 1:n){
  pval[i] <- t.test(x[i,1:12], x[i,13:24])$p.value
}


fdr_genes <- qvalue(p = pval, lambda = 0, fdr.level = qStar)$significant
qvalue_genes <- qvalue(p = pval, fdr.level = qStar)$significant

# venn it
n1 <- sum(limma_genes)
n2 <- sum(fdr_genes)

n12 <- sum(limma_genes * fdr_genes)

grid.newpage()
draw.pairwise.venn(area1 = n1, area2 = n2, cross.area = n12,
                   category = c('limma', 't-test'), lty = 3, cat.cex = 2,
                   fill = c("blue", "red"),
                   alpha = rep(0.5, 2),
                   cat.pos = c(-30,140), cat.dist = rep(0.025, 2), cex = 2, cat.col = c('blue', 'red'))

ordP <- sort(res$adj.P.Val, decreasing = F)
perm <- order(res$adj.P.Val, decreasing = F)
testcol <- color.scale(ordP, c(1,1,0),c(0,1,1),0)
plot(fit1$coefficients[perm,], col = testcol, pch = 16, cex = 0.5)
color.legend(13,3.0,14,6, c('0', '0.5', '1'),color.scale(seq(0,1,length = 10), c(1,1,0),c(0,1,1),0),gradient="y")


par(mfrow = c(1, 2))
myCol = c('palegreen4', 'firebrick')[as.numeric(limma_genes) + 1]
plot(fit1$coefficients, col = myCol, pch = 16, cex = 0.5)
points(c(0,100), c(0,100), type = 'l', lty = 3, col = 'black', lwd = 2)
plot(fit1$sigma, fit1$coefficients[,1]-fit1$coefficients[,2], pch = 16, cex = 0.5, col = myCol, xlab = 'estimated sigma', ylab = 'estimated difference')
points(c(0,10), c(0,0), type = 'l', lty = 3, col = 'black', lwd = 2)
legend('topright', c('DE', 'non-DE'), col = names(table(myCol)), pch = 16)


oasidvidsojfsjdf