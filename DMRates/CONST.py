import numpy as np
from math import erf

""" Mean Velocity """
V_MEAN = 220E3
""" Earth Velocity """
V_EARTH = 232E3
""" Escape Velocity """
V_ESCAPE = 544E3
""" Dimensionless Velocity """
A = V_ESCAPE / V_MEAN
B = V_EARTH / V_MEAN
""" Mass of Neutron in kg """
U = 1.660539040E-27
""" Energy of Neutron in GeV """
E_U = 0.9314940954
""" Reduced Planck Constant in J*s """
HBAR = 1.054571800E-34
""" Speed of Light in m/s """
C = 2.99792458E8
""" Avogaro's Number in 1/mol """
A_V = 6.022140857E23
""" Number of Neutrons per kg """
N_0 = A_V * 1E3
""" Local Density of DM in GeV/m^3 """
P_X = 0.3E6
""" Charge of Electron in C """
CH_E = 1.6021766208E-19
""" keV in J """
KEV = 1E3 * CH_E
""" GeV in J """
GEV = 1E9 * CH_E
""" MeV in J """
MEV = 1E6 * CH_E
""" Hbar times C in GeV*m """
HBARC = HBAR * C / GEV
""" Mass of GeV """
M_GEV = 1E9 * CH_E / C / C
""" One Year in s """
YEAR = 31536000
DAY = 24 * 60 * 60
""" Fermi Constant in 1/J^2 """
GF = 4.5437957E14
""" Fermi Constant in 1/GeV^2 """
GF_GEV = 1.1663787E-5
""" kg """
KG = 1
TON = 1E3 * KG
""" M """
M = 1
CM = 1E-2
KM = 1E3
""" M2 """
M2 = M ** 2
CM2 = CM ** 2
KM2 = KM ** 2
""" Pion Mass """
M_PI = 134.9766 * M_GEV * 1E-3

""" M-B Distribution Norm Factor """
N_MB = (np.sqrt(np.pi) * V_MEAN ** 2) ** 1.5
""" Sharp Cutoff Norm Factor """
N_SH = N_MB * (erf(A) - 2 / np.sqrt(np.pi) * A * np.exp(- A ** 2))



""" 
Data from https://en.wikipedia.org/wiki/Isotopes_of_xenon
"""
Xe124 = {
'Type'		: 'Xe124',
'MassNum'	: 124,
'AtomicNum'	: 54,
'Mass'		: 123.905893,
'Spin'		: 0,
'Fraction'	: 9.52E-4
}

Xe126 = {
'Type'		: 'Xe126',
'MassNum'	: 126,
'AtomicNum'	: 54,
'Mass'		: 125.904274,
'Spin'		: 0,
'Fraction'	: 8.90E-4
}

Xe128 = {
'Type'		: 'Xe128',
'MassNum'	: 128,
'AtomicNum'	: 54,
'Mass'		: 127.9035313,
'Spin'		: 0,
'Fraction'	: 0.019102
}

Xe129 = {
'Type'		: 'Xe129',
'MassNum'	: 129,
'AtomicNum'	: 54,
'Mass'		: 128.9047794,
'Spin'		: 0.5,
'Fraction'	: 0.264006
}

Xe130 = {
'Type'		: 'Xe130',
'MassNum'	: 130,
'AtomicNum'	: 54,
'Mass'		: 129.9035080,
'Spin'		: 0,
'Fraction'	: 0.040710
}

Xe131 = {
'Type'		: 'Xe131',
'MassNum'	: 131,
'AtomicNum'	: 54,
'Mass'		: 130.9050824,
'Spin'		: 1.5,
'Fraction'	: 0.212324
}

Xe132 = {
'Type'		: 'Xe132',
'MassNum'	: 132,
'AtomicNum'	: 54,
'Mass'		: 131.9041535,
'Spin'		: 0,
'Fraction'	: 0.269086
}

Xe134 = {
'Type'		: 'Xe134',
'MassNum'	: 134,
'AtomicNum'	: 54,
'Mass'		: 133.9053945,
'Spin'		: 0,
'Fraction'	: 0.104357
}

Xe136 = {
'Type'		: 'Xe136',
'MassNum'	: 136,
'AtomicNum'	: 54,
'Mass'		: 135.907219,
'Spin'		: 0,
'Fraction'	: 0.088573
}

ATOM_TABLE = {
'Xe124'	: Xe124,
'Xe126'	: Xe126,
'Xe128'	: Xe128,
'Xe129'	: Xe129,
'Xe130'	: Xe130,
'Xe131'	: Xe131,
'Xe132'	: Xe132,
'Xe134'	: Xe134,
'Xe136'	: Xe136
}

