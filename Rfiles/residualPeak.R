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

  #print(paste(sum(residuals[which(residuals > 0)])/sum(abs(residuals)), sum(residuals[which(residuals < 0)])/ sum(abs(residuals))))
  
  return(residuals)
}

findTurningPoints <- function(values){
    turningPoints = c()
    prevValue = Inf
    for(i in c(1:(120))){
        currValue = values[i] - values[i + 1]
        if(sign(prevValue) != (sign(currValue))){
            turningPoints <- append(turningPoints, i)
        }
        prevValue = currValue
    }

    return(turningPoints)
}

d       = as.numeric(substr(system(paste('grep -G \"d  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
r1      = as.numeric(substr(system(paste('grep -G \"r1  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
kmercov = as.numeric(substr(system(paste('grep -G \"kmercov  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
bias    = as.numeric(substr(system(paste('grep -G \"bias  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
length  = as.numeric(substr(system(paste('grep -G \"length  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))

residuals <- plot_residuals(d, r1, kmercov, bias, length, k, histFile)
turningPoints <- findTurningPoints(residuals)
smoothedResiduals <- smooth(residuals, twiceit = TRUE)
smoothedTurningPoints <- findTurningPoints(smoothedResiduals)

print(turningPoints)
print(smoothedTurningPoints)

lines(smoothedResiduals ~ c(1:length(smoothedResiduals)), lwd = 2, col = 'blue')

residuals[turningPoints[1]] / length