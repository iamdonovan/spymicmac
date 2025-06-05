the one script
==============

Below is a sample "one big script" to process a set of KH-9 mapping camera images from start (as ``.tgz`` files) to
finish (georeferenced, masked DEM and Orthomosaic).

All you have to do is copy the script below, save it to a file, and fill in these lines:

.. code-block:: python

    local_crs = None # need to fill in the crs to use
    glacmask = None # fill in path to optional glacier mask
    landmask = None # fill in path to optional land mask
    out_name = '' # fill in prefix for final files

Then, run the script:

.. code-block:: sh

    python the_one_script.py


the script
----------

.. code-block:: python

    #!/usr/bin/env python
    from glob import glob
    from spymicmac import data, preprocessing, micmac, register


    local_crs = None # need to fill in the crs to use
    glacmask = None # fill in path to optional glacier mask
    landmask = None # fill in path to optional land mask
    out_name = '' # fill in prefix for final files

    imlist = [fn.split('.tgz')[0] for fn in glob('*.tgz')]

    # download + reproject cop30 DEM
    data.download_cop30_vrt(imlist, crs=local_crs)

    # preprocess images
    preprocessing.preprocess_kh9_mc(
        skip='balance',
        nproc='max',
        blend=True,
        add_sfs=True,
        res_low=400,
        res_high=8000,
        add_params=True
    )

    # create the relative dem/orthophoto
    micmac.malt(
        'OIS.*tif',
        'Relative',
        dirmec='MEC-Relative',
        zoomf=2,
        cost_trans=4,
        szw=3,
        regul=0.1
    )

    # register the images to the reference DEM
    register.register_relative(
        'MEC-Relative',
        'Copernicus_DSM_ell.tif',
        glacmask=glacmask,
        landmask=landmask,
        density=250,
        strategy='chebyshev'
    )

    # create the absolute dem/orthophotos
    micmac.malt(
        'OIS.*tif',
        'TerrainFinal',
        dirmec='MEC-Malt',
        zoomf=1,
        cost_trans=4,
        szw=3,
        regul=0.1,
        resol_terr=10,
        resol_ort=2
    )

    # create the orthomosaic
    micmac.tawny('MEC-Malt')

    # clean up the outputs
    micmac.post_process(
        projstr=local_crs,
        out_name=out_name,
        dirmec='MEC-Malt',
        do_ortho=True,
        ind_ortho=True
    )
