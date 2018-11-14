import matplotlib.pyplot as plt
import numpy as np
import fitsio
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

from collapse_unwise_bitmask import collapse_unwise_bitmask

def render_remapped_bitmask(im):

    cmap = mpl.cm.Set1
    cmap.set_bad('w', 1.0)

    im_masked = np.ma.masked_array(im, mask=(im == -1))

    fig = plt.figure(frameon=False)

    # would be good to generalize this to handle non-square images properly
    fig.set_size_inches(6,6)

    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    vmin = 0
    vmax = 9

    ims = plt.imshow(im_masked[::-1, :], vmin=vmin, vmax=vmax, cmap=cmap,
                     interpolation='nearest', aspect='normal')

def remap_unwise_bitmask(mask_w1, mask_w2, band=None):

    # mask_w1 and mask_w2 should be the *collapsed* bitmasks !!!
    # could add an assert to help check for this

    # band should either be None (take the OR of both bands), 
    # 1 (use only W1 info), or 2 (use only W2 info)
    assert((band == None) or (band == 1) or (band == 2))

    if band == None:
        assert((mask_w1 is not None) and (mask_w2 is not None))
        assert(mask_w1.shape == mask_w2.shape)
        mask = np.bitwise_or(mask_w1, mask_w2)
    elif band == 1:
        assert(mask_w1 is not None)
        mask = mask_w1
    else:
        assert(mask_w2 is not None)
        mask = mask_w2

    sh = mask.shape

    # red will be diff spikes (combining psf-based and geometric)         -> 0
    # blue will be core+wings                                             -> 1
    # green will be halos                                                 -> 2
    # purple/grey will be latents (combining first and second latents)    -> 3
    # orange will be ghosts                                               -> 4

    # precedence, "bottom" to "top"
    # halos
    # spikes
    # core+wings
    # ghosts
    # latents

    # in the future use a dictionary to make the powers of 2 below
    # less confusing/obscure

    result = np.zeros(sh) - 1
    result[np.bitwise_and(mask, 2**5) != 0] = 2 # halos
    result[np.bitwise_and(mask, 2**1 + 2**7) != 0] = 0 # diff spikes
    result[np.bitwise_and(mask, 2**0) != 0] = 1 # core+wings
    result[np.bitwise_and(mask, 2**2) != 0] = 4 # ghosts
    result[np.bitwise_and(mask, 2**3 + 2**4) != 0] = 3 # latents

    return result

def remap_unwise_bitmask_native(mask, band=None):

    assert((band == None) or (band == 1) or (band == 2))

    mask_w1 = (collapse_unwise_bitmask(mask, 1) if (band != 2) else None)
    mask_w2 = (collapse_unwise_bitmask(mask, 2) if (band != 1) else None)

    return remap_unwise_bitmask(mask_w1, mask_w2, band=band)

def test_native():
    fname = '/n/fink2/ameisner/unwise-coadds/fulldepth_neo4/028/0287p272/unwise-0287p272-msk.fits.gz'
    msk = fitsio.read(fname)
    result = remap_unwise_bitmask_native(msk, band=None)
    render_remapped_bitmask(result)
    plt.savefig('unwise-0287p272-msk.png', dpi=500)

def render_unwise_bitmask(coadd_id, outdir=None, dpi=500):
    fname = '/n/fink2/ameisner/unwise-coadds/fulldepth_neo4/' + coadd_id[0:3] + '/' + coadd_id + '/unwise-' + coadd_id + '-msk.fits.gz'

    msk = fitsio.read(fname)

    outname = os.path.split(fname)[-1]
    if outdir is not None:
        outname = os.path.join(outdir, outname)

    outname = outname.replace('.fits.gz', '.png')
    print outname

    result = remap_unwise_bitmask_native(msk, band=None)
    render_remapped_bitmask(result)
    plt.savefig(outname, dpi=dpi)
    plt.close()

def decals_region():
    # (dec < 32) and (dec > -20) and (|b_gal| > 18)
    atlas = fitsio.read('/n/home09/ameisner/unwise/pro/decals-atlas.fits')

    outdir = '/n/fink2/ameisner/bitmask_renderings'

    for coadd_id in atlas['COADD_ID']:
        render_unwise_bitmask(coadd_id, outdir=outdir)
