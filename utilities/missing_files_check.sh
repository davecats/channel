for ii in $(seq $1 $2)
do 
    if [[ -f Dati.cart.$ii.out ]] # if file exists
    then
        true
    else
        echo $ii" missing!"
    fi
done