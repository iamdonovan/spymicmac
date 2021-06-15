#!/bin/bash
# post-process everything output from MicMac into nice files to use:
#	- masked DEM
#	- hillshade
#	- georeferenced orthophoto
# use: post_process_micmac.sh -z UTMZONE -n NAME -l LEVEL,
#	where utm_zone has the form "6 +north" for the projection used in processing.
utm_set=0
sub_set=0
level=9
name=HexagonOUT
dir=$(pwd)
orig_dir=$(pwd)
while getopts "z:n:l:d:h" opt; do
  case $opt in
    h)
      echo "Post-process outputs from MicMac into nice files to use."
      echo "usage: PostProcessMicMac.sh -z 'UTMZONE' -n 'NAME' -l 'LEVEL' -h"
      echo "    -z UTMZONE  : UTM Zone of area of interest. Takes form 'NN +north(south)'"
      echo "    -n NAME     : Output name to use for the DEM and related files."
      echo "    -l LEVEL    : MicMac processing level to use (default: 9)"
      echo "    -d DIR      : Directory to process (default: current directory)"
      echo "    -h          : displays this message and exits."
      echo " "
      exit 0
      ;;
    z)
      UTM=$OPTARG
      utm_set=1
      ;;
    n)
      name=$OPTARG
      ;;
    l)
      level=$OPTARG
      ;;
    d)
      dir=$OPTARG
      ;;
    \?)
      echo "PostProcessMicMac.sh: Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "PostProcessMicMac.sh: Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ $utm_set -eq 0 ]; then
      echo "Error: UTM Zone has not been set."
      echo "call PostProcessMicMac.sh -h for details on usage."
      echo " "
      exit 1
fi

resize_rasters () {
    image1=$1
    image2=$2
    # first, get the raster sizes and check that they're the same
    img1size=($(gdalinfo $image1 | grep 'Size is' | grep -o '[0-9]*'))
    img2size=($(gdalinfo $image2 | grep 'Size is' | grep -o '[0-9]*'))
    # if the rasters are the same size, we continue.
    if [[ "${img1size[0]}" -eq "${img2size[0]}" && "${img1size[1]}" -eq "${img2size[1]}" ]]; then
        echo "$image1 and $image2 are the same size. Exiting..."
        return 1
    fi
    # get the upper left and lower right corners of image1
    ul=$(gdalinfo $image1 | grep 'Upper Left' | grep -Eo '[+-]?[0-9]*\.[0-9]*\,\s*?[+-]?[0-9]*\.[0-9]*' )
    lr=$(gdalinfo $image1 | grep 'Lower Right' | grep -Eo '[+-]?[0-9]*\.[0-9]*\,\s*[+-]?[0-9]*\.[0-9]*' )
    # split into two arrays    
    ul_arr=($(echo $ul | tr , ' '))
    lr_arr=($(echo $lr | tr , ' '))
    echo "gdalwarp -te ${ul_arr[0]} ${lr_arr[1]} ${lr_arr[0]} ${ul_arr[1]} -ts ${img1size[@]} $image2 ${image2%.tif}_resize.tif"
    echo "Re-sizing $image2 to agree with $image1 size."
    gdalwarp -te ${ul_arr[0]} ${lr_arr[1]} ${lr_arr[0]} ${ul_arr[1]} -ts ${img1size[@]} $image2 ${image2%.tif}_resize.tif
    mv -v ${image2%.tif}_resize.tif $image2
}

cd $dir
mkdir -p ../post_processed

#cd MEC-Malt
clevel=$(($level-1))

cp Z_Num$level\_DeZoom1_STD-MALT.tfw Correl_STD-MALT_Num_$clevel.tfw
cp Z_Num$level\_DeZoom1_STD-MALT.tfw AutoMask_STD-MALT_Num_$clevel.tfw

gdal_translate -a_nodata 0 -a_srs "+proj=utm +zone=$UTM +ellps=WGS84 +datum=WGS84 +units=m +no_defs" Correl_STD-MALT_Num_$clevel.tif tmp_corr.tif
gdal_translate -a_srs "+proj=utm +zone=$UTM +ellps=WGS84 +datum=WGS84 +units=m +no_defs" Z_Num$level\_DeZoom1_STD-MALT.tif tmp_geo.tif
gdal_translate -a_srs "+proj=utm +zone=$UTM +ellps=WGS84 +datum=WGS84 +units=m +no_defs" -a_nodata 0 AutoMask_STD-MALT_Num_$clevel.tif tmp_msk.tif

gdal_calc.py -A tmp_msk.tif -B tmp_geo.tif --outfile=../post_processed/$name\_Z.tif --calc="B*(A>0)" --NoDataValue=-9999
gdaldem hillshade ../post_processed/$name\_Z.tif ../post_processed/$name\_HS.tif
gdal_calc.py -A tmp_corr.tif --outfile=../post_processed/$name\_CORR.tif --calc="((A.astype(float)-127)/128)*100" --NoDataValue=-9999

rm -v tmp_msk.tif tmp_geo.tif tmp_corr.tif

cd $orig_dir
