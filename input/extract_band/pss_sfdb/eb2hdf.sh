#! /bin/bash

if [ -z "$1" -o -z "$2" ];
then
    echo "Extract narrow bands from SFDB"
    echo "Usage: eb2hdf.sh <det name> <start band> [<end band>]"
fi

[ "$1" != "L1" -a "$1" != "H1" -a "$1" != "V1" ] && (echo "wrong detector name" && exit )
det=$1

#---------------  EDIT HERE ---------------

# exec location
eb=/work/chuck/virgo/scripts/data-gen/extract_band/pss_sfdb/extract_band_hdf_td

# overlap and bandwidth (see definition of fpo)
#ov=5
# overlap in O3 didn't followed 2^-ov formula, it was just 0.1, 
# remember to uncomment proper definition in fpo() function
ov=0.1
B=0.25

# list of input SFDB files
ilist=$(pwd)/sfdb_${det}.list

# output directory
odir=sts_B${B}_ov${ov}/${det}

function fpo () {
    #awk -v ov=$ov -v B=$B -v OFMT=%.15g "BEGIN{print 10.+(1.-2^-ov)*B*$1}"
    # formula used for run O3
    awk -v ov=$ov -v B=$B -v OFMT=%.15g "BEGIN{print 10.+(1.-ov)*B*$1}"
}
#-----------------------------------------

mkdir -p $odir
[[ ! -d $odir ]] && (echo "$odir is missing..."; exit)
cd $odir

# max. number of sfts: 200000 * 1024 s = 2370 days , should be enough
bstart=$2
if [ -z "$3" ];
then
    bend=$2
else
    bend=$3
fi

for band in $( seq -s' ' $bstart $bend )
do
    B4=$(printf %04d $band)
    echo -n "Extracting band $B4 ... "
    infile=${B4}_${det}.in
    echo sts_${B4}_${det}.h5 > $infile
    echo $ilist >> $infile
    f=$(fpo $band)
    echo $f >> $infile
    echo $B >> $infile
    $eb < $infile &>> eb2hdf-${B4}_${det}.out
    echo "$eb < $infile &>> eb2hdf-${B4}_${det}.out"
    echo "done!"
done
