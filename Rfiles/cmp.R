library('genomescope')
library('argparse')

parser <- ArgumentParser()
parser$add_argument("-g", "--genomescope_folder")
parser$add_argument("-i", "--hist_file")
parser$add_argument("-k", "--kmer_length", type = "integer", default = 31, help = "kmer length used to calculate kmer spectra [default 31]")
parser$add_argument("-o", "--output")

arguments = parser$parse_args()

outFolder = arguments$genomescope_folder
histFile = arguments$hist_file
k = arguments$kmer_length
out = arguments$output

# value for k can be got from the summary file, can't be bothered but will add if is useful
plot_residuals_auto <- function(outFolder, histFile, k){
  d       = as.numeric(substr(system(paste('grep -G \"d  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  r1      = as.numeric(substr(system(paste('grep -G \"r1  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  kmercov = as.numeric(substr(system(paste('grep -G \"kmercov  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  bias    = as.numeric(substr(system(paste('grep -G \"bias  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  length  = as.numeric(substr(system(paste('grep -G \"length  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))

  plot_residuals(d, r1, kmercov, bias, length, k, histFile)
}

# d, r1, kmercov, bias, length, k, file
plot_residuals <- function(d, r1, kmercov, bias, length, k, histFile){
  correctValues = read.table(file = histFile, header = FALSE, sep = "")$"V2"

  errorMin = 20

  y_limit = max(correctValues[errorMin:length(correctValues)])*1.1
  x_limit = min(c(900, which(correctValues == max(correctValues[errorMin:length(correctValues)])) * 3))

  x = c(1:length(correctValues))

  modelValues = length*predict2_0(r1, k, d, kmercov, bias, x)

  plot(correctValues, type = 'h', xlim = c(0,x_limit), ylim=c(-y_limit/20,y_limit), main = histFile)
  # points(correctValues ~ x, cex = 0.3, pch = 20)
  lines(modelValues ~ x, lwd = 2, col = 'black')
  lines((correctValues - modelValues) ~ x, lwd = 2, lty = 1, col = 'purple')

  residuals = correctValues - modelValues

  #print(paste(histFile, ": ", length(residuals[which(residuals > 0)]), " ", length(residuals[which(residuals < 0)]),sep = ""))
#
  #print(paste(sum(residuals[which(residuals > 0)])/sum(abs(residuals)), sum(residuals[which(residuals < 0)])/ sum(abs(residuals))))
  
  return(residuals)
}


pdf(out)
  residual = plot_residuals_auto(outFolder, histFile, k)
  ResidualHistogram <- hist(residual, breaks = 100000, xlim = c(-5e5,5e5))
waste <- dev.off()



length  = as.numeric(substr(system(paste('grep -G \"length  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
print(paste(histFile, ": ", sum(abs(residual / length)), sep = ""))

#sum((ResidualHistogram$counts*ResidualHistogram$mids)[which(ResidualHistogram$breaks == 0)+20:length(ResidualHistogram$counts)], na.rm = TRUE)
#sum((ResidualHistogram$counts*ResidualHistogram$mids)[1:(which(ResidualHistogram$breaks == 0) - 1)])
#
#sum(ResidualHistogram$counts*ResidualHistogram$mids)
#sum((ResidualHistogram$counts*ResidualHistogram$mids)[which(ResidualHistogram$breaks == 0):which(ResidualHistogram$breaks == 0) + 20])
#
#ResidualHistogram$mids[1:39]

# If this value is less than ~400-ish then the output probably isn't very accurate

#if(sum(ResidualHistogram$counts*sign(ResidualHistogram$mids)) < 800){
#  print(paste(histFile, ": ", sum(ResidualHistogram$counts*sign(ResidualHistogram$mids)),sep = ""))
#}


#sum(ResidualHistogram$counts)
#sum(ResidualHistogram$counts[which(ResidualHistogram$breaks <0 )])

#plot_residuals(4.476e-02, 
#               8.053e-03, 
#               2.566e+01, 
#               9.000e-01,
#               4.486e+08,
#               31, 
#               "testingData/beetles/icChrHerb2uli")
#
#plot_residuals(6.261e-02, 
#               8.911e-03, 
#               2.245e+01, 
#               5.880e-01,
#               7.005e+08,
#               31, 
#               "ORIS/genomescope-convergence/testingData/beetles/icChrGram1")
#
#plot_residuals(4.866e-02, 
#               4.339e-03, 
#               2.152e+01, 
#               5.822e-01,
#               7.493e+08,
#               31, 
#               "ORIS/genomescope-convergence/testingData/beetles/icChrOric1")
#
#plot_residuals(3.067e-02, 
#               2.639e-02, 
#               3.897e+01, 
#               1.319e-01,
#               2.781e+08,
#               31, 
#               "testingData/good/ilNapInac1.hist")
#
#plot_residuals(3.905e-01, 
#               3.015e-02, 
#               3.628e+01, 
#               2.708e+00,
#               3.341e+07,
#               31, 
#               "testingData/lowCoverage/lsAllScho1.hist")
#
#plot_residuals(1.735e-01, 
#               5.912e-02, 
#               1.936e+01, 
#               2.647e-01,
#               2.705e+08,
#               31, 
#               "ORIS/genomescope-convergence/testingData/needsWork/drAcePsed1.hist")
#
#plot_residuals(8.159e-02, 
#               9.483e-03, 
#               2.465e+01, 
#               7.202e-01,
#               5.265e+08,
#               31, 
#               "ORIS/genomescope-convergence/testingData/good/wsLamColu1.hist")
#
#plot_residuals(8.159e-02, 
#               9.483e-03, 
#               2.465e+01, 
#               7.202e-01,
#               5.265e+08,
#               31, 
#               "ORIS/genomescope-convergence/testingData/good/ilNapInac1.hist")
#
#plot_residuals(8.052e-02, 
#               2.739e-02, 
#               1.649e+01, 
#               1.042e-01,
#               3.483e+08,
#               31, 
#               "ORIS/genomescope-convergence/testingData/good/ilPseFlor1.hist")
#
#
#plot_residuals_auto("genomescopeOutput/beetles/icChrOric1",
#                    "testingData/beetles/icChrOric1", 31)