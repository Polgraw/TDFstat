#! /bin/bash

if [ -z "$1" -o -z "$2" ]; then
    echo "Generate xdat & grid files for given band(s)"
    echo "Usage: make_xdat.sh <nod> <start band> [<end band>]"
    exit
fi

nod=$1

#---------------  EDIT HERE ---------------

# exec locations
genseg="./genseg-hdf"
gridgen="./gridgen"
dets="H1 L1"
# generate grid for this segment
refseg=20

#------------------------------------------

iseg3=$(printf %03d $refseg)
echo "grid segment = ${iseg3}"

# set ephemeris generation
for det in $dets; do
    genseg_infile="${det}_bbbb_${nod}d.g2d"
    echo "genseg_infile=${genseg_infile}"
    odir="$(awk -F'=' '/^DataDir/{print $2}'  $genseg_infile | tr -d '[:blank:]')/"
    echo -n ${odir}${iseg3}/${det}/DetSSB.bin
    if [[ -f ${odir}${iseg3}/${det}/DetSSB.bin ]] ; then
	echo " exists, disabling gen_eph"
	sed -i 's/^gen_eph.*/gen_eph = no/' $genseg_infile 
    else
	echo " is missing, enabling gen_eph"
	sed -i 's/^gen_eph.*/gen_eph = yes/' $genseg_infile
    fi
done

# set range of bands to generate
bstart=$2
bend=$2
[[ -v "$3" ]] && bend=$3
echo "Band range: $bstart $bend"

declare -A Nexp2
Nexp2["1"]=18
Nexp2["2"]=19
Nexp2["4"]=20
Nexp2["6"]=20
Nexp2["8"]=21
Nexp2["10"]=21
Nexp2["12"]=21
Nexp2["16"]=22
Nexp2["24"]=22

sup_nod=$( echo ${!Nexp2[@]} | xargs -n1 | sort -n | xargs )
if [[ $sup_nod =~ $nod ]]; then
    nexp=${Nexp2[$nod]}
    echo "N = 2^$nexp"
else
    echo "Unsupported nod... Supported values: $sup_nod "
    exit
fi

# create logs dir
mkdir -p logs

# DO IT
for band in $(seq -s' ' $bstart $bend); do
    B4=$(printf %04d $band)
    echo -n "Generating band $B4 : [genseg:"
    for det in $dets; do
	echo -n " $det"
	genseg_infile="${det}_bbbb_${nod}d.g2d"
	$genseg $genseg_infile $B4 >& logs/genseg-${nod}d-${det}-${B4}.log
	#echo "$genseg $genseg_infile $B4 >& logs/genseg-${det}-${B4}.log"
    done
    echo -n "][gridgen]"
    $gridgen -c 0.75 -a s1 -s 1 -n 20 -d xdat_O4_${nod}d/ $iseg3 ${dets} $B4 >& logs/gridgen-${nod}d-${B4}.log
    #echo "$gridgen -c 0.75 -a s1 -s 1 -n $nexp -d $odir $iseg3 $dets $B4"
    echo " done"
done
