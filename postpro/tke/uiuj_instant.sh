if [ $1 == -h ]
then
    echo "Calculates instantaneous MKE/TKE statistics using uiuj. Syntax:"
    echo
    echo "   uiuj_instant.sh nmin nmax"
    echo
    echo "The average field is calculated on fields from nmin to nmax; then, statistics are calculated for each field, without averaging. They are then stored in bin files in folder instant_profiles."
fi

if [[ -d instant_profiles ]]
then
    true
else
    mkdir instant_profiles
fi

dir=$(dirname ${BASH_SOURCE[0]})/

for ii in $(seq $1 $2)
do
    ${dir}uiuj 1 1 localhost $ii $ii 1 --custom_mean $1 $2 1
    mv cm_profiles/uiuj.bin instant_profiles/uiuj.${ii}.bin
    ${dir}uiuj2ascii 1 1 localhost -n $ii
done