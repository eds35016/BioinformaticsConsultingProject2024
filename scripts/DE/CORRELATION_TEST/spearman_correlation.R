# Set CRAN repository
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check and install readr if necessary
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
library (dplyr)



