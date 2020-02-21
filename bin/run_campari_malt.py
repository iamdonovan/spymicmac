import subprocess
from glob import glob


flist = sorted(glob('OIS*.tif'))

for i, im1 in enumerate(flist[:-1]):
    im2 = flist[i+1]
    print("{}, {}".format(im1, im2))
    im_pattern = '{}|{}'.format(im1, im2)
    
    # run campari
    subprocess.Popen(['mm3d', 'Campari', im_pattern, 'TerrainRelAuto',
                      'TerrainFinal_block{}'.format(i), 'GCP=[AutoGCPs.xml,5,AutoMeasures-S2D.xml,2]',
                      'SH=Homol', 'AllFree=1']).wait()

    # check residuals, re-run campari with large errors removed
    

    # run malt
    subprocess.Popen(['mm3d', 'Malt', 'Ortho', im_pattern, 'TerrainFinal_block{}'.format(i),
                      'NbVI=2', 'ZoomF=1', 'DirMEC=MEC-Malt_block{}'.format(i), 'DefCor=0',
                      'CostTrans=4', 'EZA=1'])
