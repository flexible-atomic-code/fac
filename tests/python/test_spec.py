from pfac.atom import atomic_data
from pfac.spm import *
from pfac import const
import os


THIS_DIR = os.path.abspath(os.path.join(__file__, os.pardir))


def test_spectrum():
    # test with hydrogen like Lithium
    asym = 'Li'

    # generate atomic data for H-like to Na-like ions.
    dir0 = THIS_DIR + '/data/'
    if not os.path.exists(dir0):
        os.mkdir(dir0)

    neles = list(range(1, 2))
    atomic_data(neles, asym, iprint=1, dir=dir0)

    z = 3
    den = [1.0]

    for nele in neles:
        dir1 = THIS_DIR + '/Li%02d/'%nele
        if not os.path.exists(dir1):
            os.mkdir(dir1)

        (temp, logt, population) = get_tgrid(z, nele, dt=1)
        print('NION = 1')
        spectrum([nele], temp, den, population,
                 'Li', dir0=dir0, dir1=dir1, nion=1)
