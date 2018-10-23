import numpy as np
from math import erf
from scipy.special import spherical_jn

from DMRates.CONST import *

class TARGET(object):
	def __init__(self, ATOM, M_X, pure = False):
		super(TARGET, self).__init__()
		""" Target Type """
		self.TAR = ATOM
		""" Spin of Nucleus """
		self.SPIN = ATOM['Spin']
		""" Isotopic Abundance """
		self.WEIGHT = 1 if pure else ATOM['Fraction']
		""" Atomic Number of Target """
		self.A_TAR = ATOM['MassNum']
		""" Mass Numbert of Target """
		self.Z_TAR = ATOM['AtomicNum']
		""" Mass of Target Atom """
		self.M_TAR = ATOM['Mass'] * U
		""" Mass of Dark Matter Atom """
		self.M_X = M_X
		""" Reduced Mass """
		self.M_X_TAR = self.M_TAR * self.M_X / (self.M_TAR + self.M_X)
		self.M_TAR_X = self.M_TAR * self.M_X / (self.M_TAR + self.M_X)
		self.M_X_P = self.M_X * U / (self.M_X + U)
		self.M_P_X = self.M_X * U / (self.M_X + U)
		""" Mass Factor """
		self.r = 4 * self.M_X * self.M_TAR / (self.M_X + self.M_TAR) ** 2

		""" Mass in GeV """
		self.E_TAR = self.M_TAR / U * E_U
		self.E_X = self.M_X / M_GEV
		self.E_X_TAR = self.E_TAR * self.E_X / (self.E_TAR + self.E_X)
		self.E_TAR_X = self.E_TAR * self.E_X / (self.E_TAR + self.E_X)
		self.E_X_P = self.M_X_P / U * E_U
		self.E_P_X = self.M_X_P / U * E_U
		""" Mean Energy in J """
		self.E_MEAN = 0.5 * self.M_X * V_MEAN ** 2
		""" Max Recoil Energy in J """
		self.E_R_MAX = 2 * self.M_X_TAR ** 2 * (V_ESCAPE + V_EARTH) ** 2 / self.M_TAR

		""" Radius of Target in m """
		self.R_TAR = (self.A_TAR * 1.0) ** (1.0 / 3.0) * 1.0E-15

		""" Local Num Density of DM in 1/m^3 """
		self.N_X = P_X / self.E_X

		""" Minimum Velocity """
		self.V_MIN = lambda E_R: np.sqrt(E_R / self.E_MEAN / self.r) * V_MEAN

		""" Dimensionless Minimum Velocity """
		self.C = lambda E_R: self.V_MIN(E_R) / V_MEAN


	""" Momentum Transferred """
	def TransMoment(self, E_R):
		return np.sqrt(2 * self.M_TAR * E_R)

	""" 2pi * de Broglie Wavelength """
	def MomentToLength(self, q):
		if q == 0:
			return 1E99
		return HBAR / q

	""" SD Cross Section to a0 and a1 """
	def XSecToCoef(self, sigmap, sigman):
		ap = np.sqrt(sigmap / (3 * GF_GEV ** 2 * self.E_X_P ** 2) * (2 * np.pi) / HBARC ** 2)
		an = np.sqrt(sigman / (3 * GF_GEV ** 2 * self.E_X_P ** 2) * (2 * np.pi) / HBARC ** 2)
		a0 = ap + an
		a1 = ap - an
		return a0, a1

	""" SI Form Factor using Helm Factor """
	def SI_FormFactor(self, q):
		if q == 0:
			return 1.0
		a = 0.52E-15
		c = (1.23 * self.A_TAR ** (1.0 / 3.0) - 0.60) * 1E-15
		s = 0.9E-15
		rn = np.sqrt(c ** 2 + 7 / 3 * np.pi ** 2 * a ** 2 - 5 * s ** 2)
		Q = q / HBAR
		qrn = Q * rn
		F_TAR2 = (3 * spherical_jn(1, qrn) / qrn) ** 2 * np.exp(- (Q * s) ** 2 / 2)
		return F_TAR2

	""" SD Structure Function For XENON From arxiv:1304.7684 """
	def SD_StructureFactor(self, q):
		lam = self.MomentToLength(q)

		if self.TAR['Type'] == 'Xe131':
			u = (2.2905E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			Sp_min = coef * np.sum(np.array([1.59352E-3, -2.07344E-3, 5.67412E-3, -6.05643E-3, 3.37794E-3, -6.88135E-4, -3.42717E-5, 3.13222E-5, -4.02617E-6, 1.72711E-7]) * \
						   		   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			Sp_max = coef * np.sum(np.array([5.29643E-3, -5.28808E-3, -6.27452E-3, 2.27436E-2, -1.92229E-2, 8.44826E-3, -2.12755E-3, 3.03972E-4, -2.27893E-5, 7.05661E-7]) * \
								   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			Sn_min = coef * np.sum(np.array([1.11627E-1, -3.08602E-1, 4.74842E-1, -3.75201E-1, 1.82382E-1, -5.39711E-2, 9.44180E-3, -9.34456E-4, 4.73386E-5, -9.01514E-7]) * \
								   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			Sn_max = coef * np.sum(np.array([1.36735E-1, -3.93930E-1, 6.17924E-1, -4.88443E-1, 2.34645E-1, -6.81357E-2, 1.16393E-2, -1.11487E-3, 5.34878E-5, -9.03594E-7]) * \
								   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S00 = coef * np.sum(np.array([4.17857E-2, -1.11132E-1, 1.71306E-1, -1.32481E-1, 6.30161E-2, -1.77684E-2, 2.82192E-3, -2.32247E-4, 7.81471E-6, 1.25984E-9]) * \
								np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S11_min = coef * np.sum(np.array([1.67361E-2, -4.72853E-2, 6.84924E-2, -5.14413E-2, 2.37858E-2, -6.92778E-3, 1.24370E-3, -1.31617E-4, 7.46669E-6, -1.73484E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S11_max = coef * np.sum(np.array([2.71052E-2, -8.12985E-2, 1.22960E-1, -9.40491E-2, 4.39746E-2, -1.28013E-2, 2.27407E-3, -2.35642E-4, 1.28691E-5, -2.77011E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S01_min = coef * np.sum(np.array([-6.75438E-2, 1.95710E-1, -3.06688E-1, 2.43678E-1, -1.18395E-1, 3.51428E-2, -6.22577E-3, 6.31685E-4, -3.33272E-5, 6.82500E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S01_max = coef * np.sum(np.array([-5.29487E-2, 1.46987E-1, -2.25003E-1, 1.79499E-1, -8.88278E-2, 2.71514E-2, -4.99280E-3, 5.31148E-4, -2.99162E-5, 6.81902E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
		
		elif self.TAR['Type'] == 'Xe129':
			u = (2.2853E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			Sp_min = coef * np.sum(np.array([1.96369E-3, -1.19154E-3, -3.24210E-3, 6.22602E-3, -4.96653E-3, 2.24469E-3, -5.74412E-4, 8.31313E-5, -6.41114E-6, 2.07744E-7]) * \
						   		   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			Sp_max = coef * np.sum(np.array([7.15281E-3, -1.34790E-2, 7.88823E-3, 3.11153E-3, -6.53771E-3, 3.75478E-3, -1.05558E-3, 1.59440E-4, -1.25055E-5, 4.04987E-7]) * \
								   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			Sn_min = coef * np.sum(np.array([1.46535E-1, -4.09290E-1, 5.21423E-1, -3.74011E-1, 1.62155E-1, -4.24842E-2, 6.74911E-3, -6.33434E-4, 3.20266E-5, -6.54245E-7]) * \
								   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			Sn_max = coef * np.sum(np.array([1.79056E-1, -5.08334E-1, 6.57560E-1, -4.77988E-1, 2.09437E-1, -5.54186E-2, 8.89251E-3, -8.42977E-4, 4.30517E-5, -8.88774E-7]) * \
								   np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S00 = coef * np.sum(np.array([5.47144E-2, -1.46407E-1, 1.80603E-1, -1.25526E-1, 5.21484E-2, -1.26363E-2, 1.76284E-3, -1.32501E-4, 4.23423E-6, -1.68052E-9]) * \
								np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S11_min = coef * np.sum(np.array([2.21559E-2, -6.56100E-2, 8.63920E-2, -6.31729E-2, 2.78792E-2, -7.56661E-3, 1.26767E-3, -1.27755E-4, 7.10322E-6, -1.67272E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S11_max = coef * np.sum(np.array([3.57742E-2, -1.07895E-1, 1.45055E-1, -1.08549E-1, 4.90401E-2, -1.36169E-2, 2.33283E-3, -2.39926E-4, 1.35553E-5, -3.21404E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S01_min = coef * np.sum(np.array([-8.85644E-2, 2.54049E-1, -3.32322E-1, 2.44981E-1, -1.09298E-1, 2.96705E-2, -4.92657E-3, 4.88467E-4, -2.65022E-5, 5.98909E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
			S01_max = coef * np.sum(np.array([-6.96691E-2, 1.97380E-1, -2.54839E-1, 1.85896E-1, -8.25294E-2, 2.24322E-2, -3.75109E-3, 3.77179E-4, -2.09510E-5, 4.92362E-7]) * \
									np.array([1, u, u**2, u**3, u**4, u**5, u**6, u**7, u**8, u**9]))
		else:
			return {'p': 0, 'n': 0, '00': 0, '11': 0, '01': 0}
		return {'p': np.mean([Sp_min, Sp_max]), 'n': np.mean([Sn_min, Sn_max]), '00': S00, '11': np.mean([S11_min, S11_max]), '01': np.mean([S01_min, S01_max]), 'p_up': Sp_max, 'p_dn': Sp_min, 'n_up': Sn_max, 'n_dn': Sn_min, '11_up': S11_max, '11_dn': S11_min, '01_up': S01_max, '01_dn': S01_min}

	""" SD Cross Section """
	def SD_CrossSection(self, E_R, type_, sigma0, external_StructureFactor = None):
		q = self.TransMoment(E_R)
		if external_StructureFactor == None:
			StructureFactor = self.SD_StructureFactor
		else:
			StructureFactor = external_StructureFactor

		isospin = type_.split('_')[0]
		type___ = type_.split('_')[1]

		a = {'p': self.XSecToCoef(sigma0, 0)[0], 'n': self.XSecToCoef(0, sigma0)[0]}[isospin]
		if type___ == 'mean':
			S = a ** 2 * StructureFactor(q)[isospin]
		elif type___ == 'upper':
			S = a ** 2 * StructureFactor(q)[isospin + '_up']
		elif type___ == 'lower':
			S = a ** 2 * StructureFactor(q)[isospin + '_dn']

		return 2 * GF_GEV ** 2 / (2 * self.SPIN + 1) * self.E_X_TAR ** 2 * S * HBARC ** 2 * self.WEIGHT

	""" SI Cross Section """
	def SI_CrossSection(self, E_R, XSEC0 = 1E-49, external_FormFactor = None):
		q = self.TransMoment(E_R)
		if external_FormFactor == None:
			FormFactor = self.SI_FormFactor
		else:
			FormFactor = external_FormFactor
		return self.WEIGHT * XSEC0 * (self.M_X_TAR / self.M_X_P) ** 2 * self.A_TAR ** 2 * FormFactor(q)

	""" Mean Reverse Relative Velocity """
	def MRVel(self, E_R):
		N = erf(A) - 2 * A * np.exp(- A ** 2) / np.sqrt(np.pi)
		C = self.C(E_R)
		if A < B and C < np.abs(A - B):
			return 1 / V_MEAN / B
		elif A > B and C < np.abs(A - B):
			return 1 / (2 * N * V_MEAN * B) * (erf(C + B) - erf(C - B) - 4 / np.sqrt(np.pi) * B * np.exp(- A ** 2))
		elif np.abs(B - C) < C and C < B + A:
			return 1 / (2 * N * V_MEAN * B) * (erf(A) - erf(C - B) - 2 / np.sqrt(np.pi) * (A + B - C) * np.exp(- A ** 2))
		else:
			return 0

	""" SD Rate := dR / dE_R """
	def SD_DifRate(self, E_R, type_, sigma0, external_StructureFactor = None):
		return self.SD_CrossSection(E_R, type_, sigma0, external_StructureFactor) * self.N_X * self.MRVel(E_R) / 2 / self.M_X_TAR ** 2
	
	""" SI Rate := dR / dE_R """
	def SI_DifRate(self, E_R, XSEC0, external_FormFactor = None):
		return self.SI_CrossSection(E_R, XSEC0, external_FormFactor) * self.N_X * self.MRVel(E_R) / 2 / self.M_X_TAR ** 2


def XeSpinDependentRate(MChi_GeV, SigmaNucleon, Type):
	Xe129 = TARGET(ATOM_TABLE['Xe129'], MChi_GeV * M_GEV, False)
	Xe131 = TARGET(ATOM_TABLE['Xe131'], MChi_GeV * M_GEV, False)
	return lambda Er: Xe129.SD_DifRate(Er, Type, SigmaNucleon) + Xe131.SD_DifRate(Er, Type, SigmaNucleon)

def XeSpinIndependentRate(MChi_GeV, SigmaNucleon):
	Xe124 = TARGET(ATOM_TABLE['Xe124'], MChi_GeV * M_GEV, False)
	Xe126 = TARGET(ATOM_TABLE['Xe126'], MChi_GeV * M_GEV, False)
	Xe128 = TARGET(ATOM_TABLE['Xe128'], MChi_GeV * M_GEV, False)
	Xe129 = TARGET(ATOM_TABLE['Xe129'], MChi_GeV * M_GEV, False)
	Xe130 = TARGET(ATOM_TABLE['Xe130'], MChi_GeV * M_GEV, False)
	Xe131 = TARGET(ATOM_TABLE['Xe131'], MChi_GeV * M_GEV, False)
	Xe132 = TARGET(ATOM_TABLE['Xe132'], MChi_GeV * M_GEV, False)
	Xe134 = TARGET(ATOM_TABLE['Xe134'], MChi_GeV * M_GEV, False)
	Xe136 = TARGET(ATOM_TABLE['Xe136'], MChi_GeV * M_GEV, False)
	return lambda Er: Xe124.SI_DifRate(Er, SigmaNucleon) + \
					  Xe126.SI_DifRate(Er, SigmaNucleon) + \
					  Xe128.SI_DifRate(Er, SigmaNucleon) + \
					  Xe129.SI_DifRate(Er, SigmaNucleon) + \
					  Xe130.SI_DifRate(Er, SigmaNucleon) + \
					  Xe131.SI_DifRate(Er, SigmaNucleon) + \
					  Xe132.SI_DifRate(Er, SigmaNucleon) + \
					  Xe134.SI_DifRate(Er, SigmaNucleon) + \
					  Xe136.SI_DifRate(Er, SigmaNucleon)

class Pion(TARGET):
	def __init__(self, ATOM, M_X, pure = False):
		TARGET.__init__(self, ATOM, M_X, pure = False)
		self.M_X_PI = M_PI * self.M_X / (M_PI + self.M_X)
		self.M_PI_X = M_PI * self.M_X / (M_PI + self.M_X)

	def PionFormFactor(self, q):
		lam = self.MomentToLength(q)
		if self.TAR['Type'] == 'Xe128':
			u = (2.2847E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.42605, 2.01883, -0.576294, 0.077613, -0.00519097, 1.39081E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		elif self.TAR['Type'] == 'Xe129':
			u = (2.2873E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.44233, 2.03693, -0.579809, 0.0775201, -0.00512894, 1.35327E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		elif self.TAR['Type'] == 'Xe130':
			u = (2.2899E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.45715, 2.063, -0.594377, 0.0810307, -0.0055788, 1.59249E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		elif self.TAR['Type'] == 'Xe131':
			u = (2.2925E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.47546, 2.08643, -0.602812, 0.0824072, -0.00570646, 1.65335E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		elif self.TAR['Type'] == 'Xe132':
			u = (2.2950E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.49308, 2.11087, -0.612728, 0.0844652, -0.00597987, 1.82198E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		elif self.TAR['Type'] == 'Xe134':
			u = (2.3001E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.52965, 2.15556, -0.62789, 0.0863288, -0.00602651, 1.78002E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		elif self.TAR['Type'] == 'Xe136':
			u = (2.3051E-15 / lam) ** 2 / 2.0
			coef = np.exp(-u)
			FF = np.sum(np.array([-2.56752, 2.19645, -0.642445, 0.0883411, -0.00611004, 1.75076E-4]) * \
						np.array([1, u, u ** 2, u ** 3, u ** 4, u ** 5])) * coef
		return FF

	def Pion_DifRate(self, E_R, XSEC0):
		q = self.TransMoment(E_R)
		
		# return 4 * self.N_X * XSEC0 * self.PionFormFactor(q) ** 2 * np.exp(-self.V_MIN(E_R) ** 2 / V_MEAN ** 2) / np.sqrt(np.pi) / self.M_X_PI ** 2 / V_MEAN
		return 2 * self.N_X * XSEC0 / self.M_X_PI ** 2 * self.PionFormFactor(q) ** 2 * self.MRVel(E_R) * self.WEIGHT
		# return self.N_X * self.PionFormFactor(q) ** 2 * self.MRVel(E_R) * XSEC0 / (2 * self.M_X_P ** 2) * self.WEIGHT

def XePionWIMPRate(MChi_GeV, SigmaNucleon):
	Xe128 = Pion(ATOM_TABLE['Xe128'], MChi_GeV * M_GEV, False)
	Xe129 = Pion(ATOM_TABLE['Xe129'], MChi_GeV * M_GEV, False)
	Xe130 = Pion(ATOM_TABLE['Xe130'], MChi_GeV * M_GEV, False)
	Xe131 = Pion(ATOM_TABLE['Xe131'], MChi_GeV * M_GEV, False)
	Xe132 = Pion(ATOM_TABLE['Xe132'], MChi_GeV * M_GEV, False)
	Xe134 = Pion(ATOM_TABLE['Xe134'], MChi_GeV * M_GEV, False)
	Xe136 = Pion(ATOM_TABLE['Xe136'], MChi_GeV * M_GEV, False)
	return lambda Er: Xe128.Pion_DifRate(Er, SigmaNucleon) + \
					  Xe129.Pion_DifRate(Er, SigmaNucleon) + \
					  Xe130.Pion_DifRate(Er, SigmaNucleon) + \
					  Xe131.Pion_DifRate(Er, SigmaNucleon) + \
					  Xe132.Pion_DifRate(Er, SigmaNucleon) + \
					  Xe134.Pion_DifRate(Er, SigmaNucleon) + \
					  Xe136.Pion_DifRate(Er, SigmaNucleon)


