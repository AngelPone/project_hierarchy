
for i in {3..11}
do
    Rscript mortality/mortality.R $i
    Rscript mortality/ts_error.R $i
    Rscript mortality/naturalvsrandom.R $i
done
