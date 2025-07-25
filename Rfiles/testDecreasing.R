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

get_gradient_values_model_auto <- function(outFolder, k){
  d       = as.numeric(substr(system(paste('grep -G \"d  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  r1      = as.numeric(substr(system(paste('grep -G \"r1  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  kmercov = as.numeric(substr(system(paste('grep -G \"kmercov  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  bias    = as.numeric(substr(system(paste('grep -G \"bias  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))
  length  = as.numeric(substr(system(paste('grep -G \"length  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))

  get_gradient_values_model(d, r1, kmercov, bias, length, k)
}

get_gradient_values_model <- function(d, r1, kmercov, bias, length, k){
    modelValues = length*predict2_0(r1, k, d, kmercov, bias, c(1:999))
    gradient = c(0:(length(modelValues) - 1))

    for(i in c(1:(length(modelValues) - 1))){
        gradient[i] = modelValues[i + 1] - modelValues[i]
    }

    return(gradient)
}

get_gradient_values <- function(histFile){
    rawData = read.table(file = histFile, header = FALSE, sep = "")$"V2"
    gradient = c(0:(length(rawData) - 1))

    for(i in c(1:(length(rawData) - 1))){
        gradient[i] = rawData[i + 1] - rawData[i]
    }

    return(gradient[1:999])
}


modelGradient <- get_gradient_values_model_auto(outFolder, k)
dataGradient <- get_gradient_values(histFile)
length  = as.numeric(substr(system(paste('grep -G \"length  *[0-9]\" ', outFolder, '/model.txt', sep=""),intern=TRUE), 9, 17))

errorTurningPoint <- which(dataGradient > 0)[1]

normalisedGradient <- modelGradient / length

gradientResiduals <- dataGradient - modelGradient
gradientResidualsFiltered <- gradientResiduals[errorTurningPoint: length(gradientResiduals)]

normalisedGradientResiduals <- gradientResiduals / length
normalisedGradientResidualsFiltered <- normalisedGradientResiduals[errorTurningPoint:length(normalisedGradientResiduals)]


print(out)
print(paste("Gradient sum: ", sum(modelGradient), sep=""))
print(paste("Normalised gradient sum: ", sum(normalisedGradient), sep=""))
print(paste("Residual sum: ", sum(abs(gradientResidualsFiltered)), sep=""))
print(paste("Normalised residual sum: ", sum(abs(normalisedGradientResidualsFiltered)), sep=""))
print(paste("Turning point: ", errorTurningPoint, sep=""))

pdf(paste(out, ".pdf", sep=""))
    plot(modelGradient ~ c(1:length(modelGradient)), type = "l", xlim = c(0,200))
    lines(dataGradient ~ c(1:length(dataGradient)), col = "red")
    lines(gradientResiduals ~ c(1:length(gradientResiduals)))
junk <- dev.off()