#!/bin/bash
#
# set default resolution to 0.014 mm/pixel
res=0.014
while getopts "r:h" opt; do
  case $opt in
    h)
      echo "Resample KH-9 halves using mm3d ReSampFid."
      echo "usage: resamp_hexagon.sh -r 'RES'"
      echo "    -r RES      : Scan resolution to re-sample to (default: 0.014 mm/pixel)"
      echo "    -h          : displays this message and exits."
      echo " "
      exit 0
      ;;
    r)
      res=$OPTARG
      ;;
  esac
done

mm3d ReSampFid "DZB.*a.tif" $res BoxCh=[0,0,232,220]
mm3d ReSampFid "DZB.*b.tif" $res BoxCh=[-2,0,230,220]

