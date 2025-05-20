setwd("/home/j/perGENOTYPE/PerNCONDITION/histogram_modules/")
modules <- read.table("out.1.8.number")
colnames(modules) <- "Modules"
svg("histogram_modules.svg")
bp <- barplot(modules$Modules,names.arg = rownames(modules), xlab = "Module", ylab = "Genes", cex.axis = 1.8, col = "black", density = 200)
lines(bp, rep(10, length(modules$Modules)), col = "black", lty = 5)
lines(rep(28, length(modules$Modules)), bp, col = "black", lty = 5)
dev.off()

