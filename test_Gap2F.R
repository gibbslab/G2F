# Easy
install.packages(c("sybil","sybilSBML","RbioRXN"))
library(sybil)
library(sybilSBML)
library(RbioRXN)
sink("Easy_KB.log")
sessionInfo()
sink()

# Medium
sink("Medium_KB.log")
Glycolysis<- as.vector(read.csv("Glycolysis.txt",header = FALSE)[,1])
grep("H2O",Glycolysis)
Glycolysis[grep("H2O",Glycolysis)]
sink()

# Hard
sink("Hard_KB.log")
Glycolysis<- as.vector(read.csv("Glycolysis.txt",header = FALSE)[,1])
grep("<?=>?[[:print:]]+?H2O",Glycolysis)
Glycolysis[grep("<?=>?[[:print:]]+?H2O",Glycolysis)]
sink()
