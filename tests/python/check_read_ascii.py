import os
import numpy as np
from . import read_ascii


THIS_DIR = os.path.abspath(os.path.join(__file__, os.pardir))


def test_lev():
    # path = THIS_DIR + '/../../demo/aidr/ref/se.lev'
    path = THIS_DIR + '/../../demo/aidr/se.lev'
    header, blocks = read_ascii.lev(path)

    assert header['FAC'] == '1.1.4'
    
    assert header['NBlocks'] == len(blocks)
    assert header['E0'] == -5.33765518E+04
        
    # block 1
    block = blocks[1]
    
    assert block['NELE'] == 10
    assert np.allclose(block['ILEV'][:2], [2, 3])
    assert np.allclose(block['IBASE'][:2], [-1, -1])
    assert np.allclose(block['ENERGY'][:2], [3.86302556E+02, 4.05994119E+02])
    assert np.allclose(block['P'][:2], [0, 0])
    assert (block['VNL'][:2] == np.array(['300', '300'])).all()
    assert np.allclose(block['2J'][:2], [4, 0])
    assert (block['ncomplex'][:2] == ['1*2 2*6 3*2', '1*2 2*6 3*2']).all()
    assert (block['sname'][:3] == ['2p4', '2p4', '2p4 3s1 3p1']).all()
    assert (block['name'][:3] ==
            ['2p+2(4)4', '2p+2(0)0', '2p+2(4)4 3s+1(1)5 3p-1(1)4']).all()
                       

def test_tr():
    # path = THIS_DIR + '/../../demo/aidr/ref/se.lev'
    path = THIS_DIR + '/../../demo/structure/ne.tr'
    header, blocks = read_ascii.tr(path)

    assert header['FAC'] == '1.1.4'
    assert header['NBlocks'] == len(blocks)
        
    # block 1
    block = blocks[0]
    
    assert block['NELE'] == 10
    assert np.allclose(block['upper_index'][:2], [1, 2])
    assert np.allclose(block['upper_2J'][:2], [4, 2])
    assert np.allclose(block['lower_index'][:2], [0, 0])
    assert np.allclose(block['lower_2J'][:2], [0, 0])
    assert np.allclose(block['Delta E'][:2], [7.243437E+02, 7.263864E+02])
    assert np.allclose(block['gf'][:2], [4.481998E-08, 1.118383E-01])
    assert np.allclose(block['rate'][:2], [2.040800E+05, 8.535211E+11])
    assert np.allclose(block['multipole'][:2], [4.481998E-08, 1.118383E-01])

    
def test_ai():
    # path = THIS_DIR + '/../../demo/aidr/ref/se.lev'
    path = THIS_DIR + '/../../demo/aidr/se.ai'
    header, blocks = read_ascii.ai(path)

    assert header['FAC'] == '1.1.4'
    assert header['NBlocks'] == len(blocks)
        
    # block 1
    block = blocks[0]
    
    assert block['NELE'] == 10
    assert block['EMIN'] == 0.0
    assert np.allclose(block['bound_index'][:2], [2, 2])
    assert np.allclose(block['bound_2J'][:2], [4, 4])
    assert np.allclose(block['free_index'][:2], [0, 1])
    assert np.allclose(block['free_2J'][:2], [3, 1])
    assert np.allclose(block['rate'][0], [3.8630E+02, 1.3904E+13, 1.1136E+01])
    assert np.allclose(block['rate'][1], [3.4343E+02, 1.8221E+11, 3.2829E-01])


def test_ce():
    # path = THIS_DIR + '/../../demo/aidr/ref/se.lev'
    path = THIS_DIR + '/../../demo/excitation/ne.ce'
    header, blocks = read_ascii.ce(path)

    assert header['FAC'] == '1.1.4'
    assert header['NBlocks'] == len(blocks)
        
    # block 1
    block = blocks[0]
    
    assert block['NELE'] == 10
    assert np.allclose(block['TEGRID'], [7.24471154E+02, 9.43883429E+02])
    assert np.allclose(block['lower_index'][:2], [0, 0])
    assert np.allclose(block['lower_2J'][:2], [0, 0])
    assert np.allclose(block['upper_index'][:2], [1, 2])
    assert np.allclose(block['upper_2J'][:2], [4, 2])
    assert np.allclose(block['bethe'][:2], [-1E0, 8.3122E-03])
    assert np.allclose(block['born'][:2, 0], [0, -6.1412e-3])
    assert np.allclose(block['born'][:2, 1], [0, 9.4388e4])
    assert np.allclose(block['collision strength'][:2, 0], 
                       [1.5599e-3, 1.7267e-3])
    assert np.allclose(block['crosssection'][:2, 0],
                       [2.4352E-01, 2.6885E-01])
    
    
def test_ci():
    # path = THIS_DIR + '/../../demo/aidr/ref/se.lev'
    path = THIS_DIR + '/../../demo/ionization/ne.ci'
    header, blocks = read_ascii.ci(path)
    
    assert header['FAC'] == '1.1.4'
    assert header['NBlocks'] == len(blocks)
        
    # block 1
    block = blocks[0]
    
    assert block['NELE'] == 10
    assert np.allclose(block['EGRID'], 
                       [6.63657689E+01,
                        1.02302245E+03,
                        2.47985383E+03,
                        4.54487006E+03,
                        7.26452555E+03,
                        1.06185230E+04])
    assert np.allclose(block['bound_index'][:2], [0, 0])
    assert np.allclose(block['bound_2J'][:2], [0, 0])
    assert np.allclose(block['free_index'][:2], [1, 2])
    assert np.allclose(block['free_2J'][:2], [3, 1])
    assert np.allclose(block['Delta E'][:2], [1.2603e3, 1.2730e3])
    assert np.allclose(block['Delta L'][:2], [1, 1])
    assert np.allclose(block['parameters'][:2, 0], [1.3579e-1, 6.8298e-2])
    assert np.allclose(block['collision strength'][:2, 0], 
                       [1.0927e-1, 5.3465E-02])
    assert np.allclose(block['crosssection'][:2, 0],
                       [2.3611E+00, 1.1469E+00])
    
    