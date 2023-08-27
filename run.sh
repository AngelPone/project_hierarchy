
for i in {0..10}
do
    echo $i
    #Rscript mortality/mortality.R $i
    Rscript mortality/ts_error.R $i
    #Rscript mortality/naturalvsrandom.R $i
done
