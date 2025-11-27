#! /bin/bash

if [ -z "$1" -o -z "$2" ];
then
    echo "Extract narrow bands from SFDB"
    echo "Usage: eb2out.sh <det name> <start band> [<end band>]"
fi

[ "$1" != "L1" -a "$1" != "H1" -a "$1" != "V1" ] && (echo "wrong detector name" && exit )
det=$1

#---------------  EDIT HERE ---------------
# exec location
eb=/work/chuck/virgo/scripts/data-gen/extract_band/pss_sfdb/extract_band_out
# overlap and bandwidth (see definition of fpo)
overlap=0.1
B=0.25
# output directory
odir=sft_ov=${overlap}_B=${B}/${det}

function fpo () {
    #awk -v on=$overlap_n -v B=$B -v OFMT=%.15g "BEGIN{print 10.+(1.-2^-ov)*B*$1}"
    awk -v ov=$overlap -v B=$B -v OFMT=%.15g "BEGIN{print 10.+(1.-ov)*B*$1}"
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
    echo 200000 > $infile
    echo sfft_${B4}_${det}.out >> $infile
    echo ../../sfdb_files_${det}.list >> $infile
    f=$(fpo $band)
    echo $f >> $infile
    echo $B >> $infile
    $eb < $infile >& eb2out-${B4}_${det}.out
    echo "done!"
done
