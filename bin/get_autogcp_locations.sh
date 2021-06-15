#!/bin/bash
programname=$0

function usage {
      echo "Find location of automatically-detected control points in individual images using mm3d XYZ2Im."
      echo "usage: $programname Dir-ORI GCP_MEASURES_FILE <Images>"
      exit 0
}

if [ $# -lt 3 ]; then
  usage
fi

ori=$1
meas_file=$2
shift 2

cp -rv $ori $ori-NoDist
tmp_autocal=$(ls $ori-NoDist/AutoCal*)
new_autocal=$(echo $tmp_autocal | sed 's_/_\\/_g')

sed -i 's/\(<CoeffDist.*>\)[^<>]*\(<\/CoeffDist.*\)/\10.0\2/' $tmp_autocal

for im in "$@"; do
    sed -i "s/\(<FileInterne.*>\)[^<>]*\(<\/FileInterne.*\)/\1$new_autocal\2/" $ori-NoDist/Orientation-$im.xml
    mm3d XYZ2Im $ori-NoDist/Orientation-$im.xml $meas_file NoDist-$im.txt
    mm3d XYZ2Im $ori/Orientation-$im.xml $meas_file Auto-$im.txt
done

# combine_auto_measures.py Auto-*.txt --no_distortion
# rm -r $ori-NoDist
# rm NoDist-*.txt
# rm Auto-*.txt
