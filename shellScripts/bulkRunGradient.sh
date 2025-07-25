for cat in ../dataFolders/testingData/*
do
    for filename in $cat/*
    do 
        Rscript ../Rfiles/testDecreasing.R -i $filename -g ../dataFolders/genomescopeOutput/${filename#../dataFolders/testingData/} -o ../dataFolders/gradientOutput/${filename#../dataFolders/testingData/} -k 31
    done
done