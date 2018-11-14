import numpy as np
from collections import OrderedDict

# output bits :

# 2^0 = bright star core and wings
# 2^1 = PSF-based diffraction spike
# 2^2 = optical ghost
# 2^3 = first latent
# 2^4 = second latent
# 2^5 = AllWISE-like circular halo
# 2^6 = bright star saturation
# 2^7 = geometric diffraction spike

def collapse_unwise_bitmask(bitmask, band):

    assert((band == 1) or (band == 2))

    bits_w1 = OrderedDict([('core_wings', 2**0 + 2**1),
                           ('psf_spike', 2**27),
                           ('ghost', 2**25 + 2**26),
                           ('first_latent', 2**13 + 2**14), 
                           ('second_latent', 2**17 + 2**18),
                           ('circular_halo', 2**23),
                           ('saturation', 2**4),
                           ('geom_spike', 2**29)])

    bits_w2 = OrderedDict([('core_wings', 2**2 + 2**3),
                           ('psf_spike', 2**28),
                           ('ghost', 2**11 + 2**12),
                           ('first_latent', 2**15 + 2**16), 
                           ('second_latent', 2**19 + 2**20),
                           ('circular_halo', 2**24),
                           ('saturation', 2**5),
                           ('geom_spike', 2**30)])

    bits = (bits_w1 if (band == 1) else bits_w2)

    # hack to handle both scalar and array inputs
    result = 0*bitmask

    for i, feat in enumerate(bits.keys()):
        result += (2**i)*(np.bitwise_and(bitmask, bits[feat]) != 0)

    return result.astype('uint8')
