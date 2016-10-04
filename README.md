g2f : Find and Fill Gaps in Metabolic Reconstructions
======
The **g2f** package was designed as a tool to find and fill gaps in metabolic reconstructions.
Install:
--------
This package required R version 2.10 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

For install the latest stable version this package directly from GitHub:
```
# Install 'devtools' R package
install.packages("devtools")

# Install 'minval' package
devtools::install_github("gibbslab/g2f")
library("g2f")
```

Available functions:
-------------------
|Function | Description |
|:--------|:------------|
|getReference|Download all the set of gene-associated stoichiometric reactions for a specific organism from the KEGG database|

Citation
--------
Kelly Botero, Daniel Osorio, Janneth Gonzalez and Andres Pinzon-Velasco (2016). **g2f: Find and Fill Gaps in Metabolic Reconstructions**. R package version 0.1.
