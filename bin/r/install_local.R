library(devtools)


print('Installing')


devtools::install_github('pmandros/TCGAPurityFiltering')
devtools::install_github('immunogenomics/presto')
devtools::install_github('aet21/EpiSCORE')
devtools::install_github('https://github.com/QuackenbushLab/NetworkDataCompanion', force=TRUE, dependencies=T)