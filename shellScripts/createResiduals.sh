for cat in genomescopeOutput/*
do
    for filename in $cat/*
    do 
        Rscript cmp.R -i testingData/${filename#genomescopeOutput/} -g genomescopeOutput/${filename#genomescopeOutput/} -o residualOutput/${filename#genomescopeOutput/}.pdf
    done
done