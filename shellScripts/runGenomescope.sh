for cat in testingData/*
do
    for filename in $cat/*
    do 
        genomescope.R -i $filename -o genomescopeOutput/${filename#testingData/} -k 31
    done
done