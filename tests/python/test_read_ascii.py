import numpy as np
import os

from . import read_ascii


THIS_DIR = os.path.abspath(os.path.join(__file__, os.pardir))


def test_read_en():
    filename = THIS_DIR + '/../../demo/aidr/se.lev'
    header, blocks = read_ascii.en(filename)
