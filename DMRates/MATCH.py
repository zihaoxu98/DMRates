""" Match WIMP-quark EFT and WIMP-nucleon EFT up to NLO """
""" Still need to be verified """
import numpy as np

from DMRates.QCD import *

class NucleonFormFactor(object):
	def __init__(self):
		pass

	'''arxiv:1707.06998'''
	def EMFF(self, isospin):
		delta = {'p': 1, 'n': 0}[isospin]
		mN = m_nucleon[isospin]
		if isospin == 'n':
			F1 = lambda t: delta + rEn2 / 6 * t + (delta - mun) * t / 4 / mN ** 2
			F2 = lambda t: mun - delta
		elif isospin =='p':
			F1 = lambda t: delta + rEp2 / 6 * t + (delta - mup) * t / 4 / mN ** 2
			F2 = lambda t: mup - delta
		return {'F1': F1, 'F2': F2}

	'''arxiv:1605.08043'''
	def ScalarFF(self, quark, isospin):
		xi = (m_quark['d'] - m_quark['u']) / (m_quark['d'] + m_quark['u'])
		dsigma = 0.27 / GeV
		dsigmas = 0.3 / GeV
		if quark == 'u':
			if isospin == 'n':
				return lambda t: m_nucleon['n'] * fnu + (1 - xi) / 2 * dsigma * t
			elif isospin == 'p':
				return lambda t: m_nucleon['p'] * fpu + (1 - xi) / 2 * dsigma * t
		elif quark == 'd':
			if isospin == 'n':
				return lambda t: m_nucleon['n'] * fnd + (1 + xi) / 2 * dsigma * t
			elif isospin == 'p':
				return lambda t: m_nucleon['p'] * fpd + (1 + xi) / 2 * dsigma * t
		elif quark == 's':
			if isospin == 'n':
				return lambda t: m_nucleon['n'] * fns + dsigmas * t
			elif isospin == 'p':
				return lambda t: m_nucleon['p'] * fps + dsigmas * t

	'''arxiv:1707.06998'''
	def PseudoscalarFF(self, quark, isospin):
		mN = m_nucleon[isospin]
		pipole = lambda t: mN ** 2 / mpi ** 2 + t * mN ** 2 / mpi ** 4# - 2 * t ** 2 * mN ** 2 / mpi ** 6
		etpole = lambda t: mN ** 2 / met ** 2 + t * mN ** 2 / met ** 4# - 2 * t ** 2 * mN ** 2 / met ** 6
		if quark == 'u':
			if isospin == 'n':
				return self.PseudoscalarFF('d', 'p')
			elif isospin == 'p':
				api = m_quark['u'] / mN * B0 * gA
				aet = m_quark['u'] * B0 * (gS - 2 * deltas) / (3 * mN)
				return lambda t: pipole(t) * api + etpole(t) * aet
		elif quark == 'd':
			if isospin == 'n':
				return self.PseudoscalarFF('u', 'p')
			elif isospin == 'p':
				api = -m_quark['d'] / mN * B0 * gA
				aet = m_quark['d'] * B0 * (gS - 2 * deltas) / (3 * mN)
				return lambda t: pipole(t) * api + etpole(t) * aet
		elif quark == 's':
			if isospin == 'n':
				return self.PseudoscalarFF('s', 'p')
			elif isospin == 'p':
				api = 0
				aet = -2 * m_quark['s'] * B0 * (gS - 2 * deltas) / (3 * mN)
				return lambda t: pipole(t) * api + etpole(t) * aet

	'''arxiv:1707.06998'''
	def VectorFF(self, quark, isospin):
		if quark == 'u':
			if isospin == 'n':
				return self.VectorFF('d', 'p')
			elif isospin == 'p':
				F1 = lambda t: 2 + 5.56 / GeV ** 2 * t
				F2 = lambda t: 1.609
		elif quark == 'd':
			if isospin == 'n':
				return self.VectorFF('u', 'p')
			elif isospin == 'p':
				F1 = lambda t: 1 + 2.84 / GeV ** 2 * t
				F2 = lambda t: -2.097
		elif quark == 's':
			if isospin == 'n':
				return self.VectorFF('s', 'p')
			elif isospin == 'p':
				F1 = lambda t: 0 - 0.018 / GeV ** 2 * t
				F2 = lambda t: -0.064
		return {'F1': F1, 'F2': F2}

	'''arxiv:1707.06998'''
	def AxialvectorFF(self, quark, isospin):
		mN = m_nucleon[isospin]
		pipole = lambda t: mN ** 2 / mpi ** 2 + t * mN ** 2 / mpi ** 4# - 2 * t ** 2 * mN ** 2 / mpi ** 6
		etpole = lambda t: mN ** 2 / met ** 2 + t * mN ** 2 / met ** 4# - 2 * t ** 2 * mN ** 2 / met ** 6
		if quark == 'u':
			if isospin == 'n':
				return self.AxialvectorFF('d', 'p')
			elif isospin == 'p':
				FA = lambda t: deltau * (1 + 1.47 / GeV ** 2 * t)
				FP = lambda t: pipole(t) * 2 * gA + etpole(t) * 2 * (gS - 2 * deltas) / 3 - 4.65
		elif quark == 'd':
			if isospin == 'n':
				return self.AxialvectorFF('u', 'p')
			elif isospin == 'p':
				FA = lambda t: deltad * (1 + 2.47 / GeV ** 2 * t)
				FP = lambda t: -pipole(t) * 2 * gA + etpole(t) * 2 * (gS - 2 * deltas) / 3 + 3.28
		elif quark == 's':
			if isospin == 'n':
				return self.AxialvectorFF('s', 'p')
			elif isospin == 'p':
				FA = lambda t: deltas * (1 + 3.0 / GeV ** 2 * t)
				FP = lambda t: -4 / 3 * (gS - 2 * deltas) * etpole(t) - 11 * deltas
		return {'FA': FA, 'FP': FP}

	'''Phys. Lett. B. 78 (1978) 443.'''
	def CPEvenGluonFF(self, isospin):
		return lambda t: -2 / 27 * (m_nucleon[isospin] - sum([self.ScalarFF(quark, isospin)(t) for quark in ['u', 'd', 's']]))

	'''arxiv:1707.06998, Ward Identity'''
	def CPOddGluonFF(self, isospin):
		mu = 1 / (1 / m_quark['u'] + 1 / m_quark['d'] + 1 / m_quark['s'])
		mN = m_nucleon[isospin]
		return lambda t: mu * sum([self.PseudoscalarFF(quark, isospin)(t) / m_quark[quark] - \
								   self.AxialvectorFF(quark, isospin)['FA'](t) * mN / m_quark[quark] - \
								   self.AxialvectorFF(quark, isospin)['FP'](t) * t / (4 * mN * m_quark[quark]) \
								   for quark in ['u', 'd', 's']])

	'''arxiv:1707.06998'''
	def TensorFF(self, quark, isospin):
		mq = m_quark[quark]
		gq = {'u': 0.794, 'd': -0.204, 's': 3.2E-4}[quark]
		if quark == 'u':
			if isospin == 'n':
				return self.TensorFF('d', 'p')
			elif isospin == 'p':
				F0 = lambda t: mq * gq * (1 + 0.8 / GeV ** 2 * t)
				F1 = lambda t: -mq * 3.0 * (1 + 1.0 / GeV ** 2 * t)
				F2 = lambda t: mq / 2 * (-0.5) * (1 + 1.2 / GeV ** 2 * t)
		elif quark == 'd':
			if isospin == 'n':
				return self.TensorFF('u', 'p')
			elif isospin == 'p':
				F0 = lambda t: mq * gq * (1 + 0.7 / GeV ** 2 * t)
				F1 = lambda t: -mq * 0.24 * (1 - 0.1 / GeV ** 2 * t)
				F2 = lambda t: mq / 2 * 0.46 * (1 + 1.0 / GeV ** 2 * t)
		elif quark == 's': # Not Precise
			if isospin == 'n':
				return self.TensorFF('s', 'p')
			elif isospin == 'p':
				F0 = lambda t: mq * gq * (1 + 1.0 / GeV ** 2 * t)
				F1 = lambda t: -mq * (-0.2)
				F2 = lambda t: mq / 2 * 0.2
		return {'F0': F0, 'F1': F1, 'F2': F2}


class Match(NucleonFormFactor):
	"""docstring for Match"""
	def __init__(self, mchi_GeV):
		NucleonFormFactor.__init__(self)
		self.mchi = mchi_GeV * GeV
		self.mu = {'n': mun, 'p': mup}
		self.QuarkList = ['u', 'd', 's']
		self.QuarkBasis = {'Q5,1': 0,\
						   'Q5,2': 0,\
						   'Q6,1_u': 0, 'Q6,1_d': 0, 'Q6,1_s': 0,\
						   'Q6,2_u': 0, 'Q6,2_d': 0, 'Q6,2_s': 0,\
						   'Q6,3_u': 0, 'Q6,3_d': 0, 'Q6,3_s': 0,\
						   'Q6,4_u': 0, 'Q6,4_d': 0, 'Q6,4_s': 0,\
						   'Q7,1': 0,\
						   'Q7,2': 0,\
						   'Q7,3': 0,\
						   'Q7,4': 0,\
						   'Q7,5_u': 0, 'Q7,5_d': 0, 'Q7,5_s': 0,\
						   'Q7,6_u': 0, 'Q7,6_d': 0, 'Q7,6_s': 0,\
						   'Q7,7_u': 0, 'Q7,7_d': 0, 'Q7,7_s': 0,\
						   'Q7,8_u': 0, 'Q7,8_d': 0, 'Q7,8_s': 0,\
						   'Q7,9_u': 0, 'Q7,9_d': 0, 'Q7,9_s': 0,\
						   'Q7,10_u': 0, 'Q7,10_d': 0, 'Q7,10_s': 0}
		self.NeutronBasis = [lambda t: 0 for i in range(13)]
		self.ProtonBasis = [lambda t: 0 for i in range(13)]
		self.c1p = lambda t: 0
		self.c4p = lambda t: 0
		self.c5p = lambda t: 0
		self.c6p = lambda t: 0
		self.c7p = lambda t: 0
		self.c8p = lambda t: 0
		self.c9p = lambda t: 0
		self.c10p = lambda t: 0
		self.c11p = lambda t: 0
		self.c12p = lambda t: 0
		self.c1n = lambda t: 0
		self.c4n = lambda t: 0
		self.c5n = lambda t: 0
		self.c6n = lambda t: 0
		self.c7n = lambda t: 0
		self.c8n = lambda t: 0
		self.c9n = lambda t: 0
		self.c10n = lambda t: 0
		self.c11n = lambda t: 0
		self.c12n = lambda t: 0

	def SetCoef(self, key, value):
		self.QuarkBasis[key] = value

	def ZeroCoef(self):
		for key in self.QuarkBasis:
			self.QuarkBasis[key] = 0
		self.NeutronBasis = [lambda t: 0 for i in range(13)]
		self.ProtonBasis = [lambda t: 0 for i in range(13)]

	def Calculate(self):
		self.c1p = lambda t: -alpha / np.pi / 2 / self.mchi * self.QuarkBasis['Q5,1'] * self.EMFF('p')['F1'](t) +\
						sum([self.VectorFF(quark, 'p')['F1'](t) * self.QuarkBasis['Q6,1_' + quark] + self.ScalarFF(quark, 'p')(t) * self.QuarkBasis['Q7,5_' + quark] for quark in self.QuarkList]) +\
						self.QuarkBasis['Q7,1'] * self.CPEvenGluonFF('p')(t) -\
						t / 2 / self.mchi / m_nucleon['p'] * sum([(self.TensorFF(quark, 'p')['F0'](t) - self.TensorFF(quark, 'p')['F1'](t)) * self.QuarkBasis['Q7,9_' + quark] for quark in self.QuarkList])
		self.c4p = lambda t: -2 * alpha * self.mu['p'] * self.QuarkBasis['Q5,1'] / np.pi / m_nucleon['p'] +\
						sum([8 * self.TensorFF(quark, 'p')['F0'](t) * self.QuarkBasis['Q7,9_' + quark] - 4 * self.AxialvectorFF(quark, 'p')['FA'](t) * self.QuarkBasis['Q6,4_' + quark] for quark in self.QuarkList])
		self.c5p = lambda t: 2 * alpha * m_nucleon['p'] * self.QuarkBasis['Q5,1'] / np.pi / t * self.EMFF('p')['F1'](t) if self.QuarkBasis['Q5,1'] else 0.0
		self.c6p = lambda t: (2 * alpha * self.mu['p'] * m_nucleon['p'] * self.QuarkBasis['Q5,1'] / np.pi / t if self.QuarkBasis['Q5,1'] else 0) +\
						sum([self.AxialvectorFF(quark, 'p')['FP'](t) * self.QuarkBasis['Q6,4_' + quark] + m_nucleon['p'] / self.mchi * self.PseudoscalarFF(quark, 'p')(t) * self.QuarkBasis['Q7,8_' + quark] for quark in self.QuarkList]) +\
						m_nucleon['p'] / self.mchi * self.CPOddGluonFF('p')(t) * self.QuarkBasis['Q7,4']
		self.c7p = lambda t: -2 * sum([self.AxialvectorFF(quark, 'p')['FA'](t) * self.QuarkBasis['Q6,3_' + quark] for quark in self.QuarkList])
		self.c8p = lambda t: 2 * sum([self.VectorFF(quark, 'p')['F1'](t) * self.QuarkBasis['Q6,2_' + quark] for quark in self.QuarkList])
		self.c9p = lambda t: 2 * sum([(self.VectorFF(quark, 'p')['F1'](t) + self.VectorFF(quark, 'p')['F2'](t)) * self.QuarkBasis['Q6,2_' + quark] +\
								 m_nucleon['p'] / self.mchi * self.AxialvectorFF(quark, 'p')['FA'](t) * self.QuarkBasis['Q6,3_' + quark] for quark in self.QuarkList])
		self.c10p = lambda t: self.CPOddGluonFF('p')(t) * self.QuarkBasis['Q7,3'] + sum([self.PseudoscalarFF(quark, 'p')(t) * self.QuarkBasis['Q7,7_' + quark] - \
																					2 * m_nucleon['p'] / self.mchi * self.TensorFF(quark, 'p')['F0'](t) * self.QuarkBasis['Q7,10_' + quark] for quark in self.QuarkList])
		self.c11p = lambda t: (2 * alpha * m_nucleon['p'] / np.pi / t * self.QuarkBasis['Q5,2'] * self.EMFF('p')['F1'](t) if self.QuarkBasis['Q5,2'] else 0) +\
						 sum([2 * (self.TensorFF(quark, 'p')['F0'](t) - self.TensorFF(quark, 'p')['F1'](t)) * self.QuarkBasis['Q7,10_' + quark] -\
							  m_nucleon['p'] / self.mchi * self.ScalarFF(quark, 'p')(t) * self.QuarkBasis['Q7,6_' + quark] for quark in self.QuarkList]) -\
						 m_nucleon['p'] / self.mchi * self.CPEvenGluonFF('p')(t) * self.QuarkBasis['Q7,2']
		self.c12p = lambda t: -8 * sum([self.TensorFF(quark, 'p')['F0'](t) * self.QuarkBasis['Q7,10_' + quark] for quark in self.QuarkList])


		self.c1n = lambda t: -alpha / np.pi / 2 / self.mchi * self.QuarkBasis['Q5,1'] * self.EMFF('n')['F1'](t) +\
						sum([self.VectorFF(quark, 'n')['F1'](t) * self.QuarkBasis['Q6,1_' + quark] + self.ScalarFF(quark, 'n')(t) * self.QuarkBasis['Q7,5_' + quark] for quark in self.QuarkList]) +\
						self.QuarkBasis['Q7,1'] * self.CPEvenGluonFF('n')(t) -\
						t / 2 / self.mchi / m_nucleon['n'] * sum([(self.TensorFF(quark, 'n')['F0'](t) - self.TensorFF(quark, 'n')['F1'](t)) * self.QuarkBasis['Q7,9_' + quark] for quark in self.QuarkList])
		self.c4n = lambda t: -2 * alpha * self.mu['n'] * self.QuarkBasis['Q5,1'] / np.pi / m_nucleon['n'] +\
						sum([8 * self.TensorFF(quark, 'n')['F0'](t) * self.QuarkBasis['Q7,9_' + quark] - 4 * self.AxialvectorFF(quark, 'n')['FA'](t) * self.QuarkBasis['Q6,4_' + quark] for quark in self.QuarkList])
		self.c5n = lambda t: 2 * alpha * m_nucleon['n'] * self.QuarkBasis['Q5,1'] / np.pi / t * self.EMFF('n')['F1'](t) if self.QuarkBasis['Q5,1'] else 0.0
		self.c6n = lambda t: (2 * alpha * self.mu['n'] * m_nucleon['n'] * self.QuarkBasis['Q5,1'] / np.pi / t if self.QuarkBasis['Q5,1'] else 0) +\
						sum([self.AxialvectorFF(quark, 'n')['FP'](t) * self.QuarkBasis['Q6,4_' + quark] + m_nucleon['n'] / self.mchi * self.PseudoscalarFF(quark, 'n')(t) * self.QuarkBasis['Q7,8_' + quark] for quark in self.QuarkList]) +\
						m_nucleon['n'] / self.mchi * self.CPOddGluonFF('n')(t) * self.QuarkBasis['Q7,4']
		self.c7n = lambda t: -2 * sum([self.AxialvectorFF(quark, 'n')['FA'](t) * self.QuarkBasis['Q6,3_' + quark] for quark in self.QuarkList])
		self.c8n = lambda t: 2 * sum([self.VectorFF(quark, 'n')['F1'](t) * self.QuarkBasis['Q6,2_' + quark] for quark in self.QuarkList])
		self.c9n = lambda t: 2 * sum([(self.VectorFF(quark, 'n')['F1'](t) + self.VectorFF(quark, 'n')['F2'](t)) * self.QuarkBasis['Q6,2_' + quark] +\
								 m_nucleon['n'] / self.mchi * self.AxialvectorFF(quark, 'n')['FA'](t) * self.QuarkBasis['Q6,3_' + quark] for quark in self.QuarkList])
		self.c10n = lambda t: self.CPOddGluonFF('n')(t) * self.QuarkBasis['Q7,3'] + sum([self.PseudoscalarFF(quark, 'n')(t) * self.QuarkBasis['Q7,7_' + quark] - \
																					2 * m_nucleon['n'] / self.mchi * self.TensorFF(quark, 'n')['F0'](t) * self.QuarkBasis['Q7,10_' + quark] for quark in self.QuarkList])
		self.c11n = lambda t: (2 * alpha * m_nucleon['n'] / np.pi / t * self.QuarkBasis['Q5,2'] * self.EMFF('n')['F1'](t) if self.QuarkBasis['Q5,2'] else 0) +\
						 sum([2 * (self.TensorFF(quark, 'n')['F0'](t) - self.TensorFF(quark, 'n')['F1'](t)) * self.QuarkBasis['Q7,10_' + quark] -\
							  m_nucleon['n'] / self.mchi * self.ScalarFF(quark, 'n')(t) * self.QuarkBasis['Q7,6_' + quark] for quark in self.QuarkList]) -\
						 m_nucleon['n'] / self.mchi * self.CPEvenGluonFF('n')(t) * self.QuarkBasis['Q7,2']
		self.c12n = lambda t: -8 * sum([self.TensorFF(quark, 'n')['F0'](t) * self.QuarkBasis['Q7,10_' + quark] for quark in self.QuarkList])

		self.ProtonBasis[1] = lambda qGeV: self.c1p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[4] = lambda qGeV: self.c4p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[5] = lambda qGeV: self.c5p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[6] = lambda qGeV: self.c6p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[7] = lambda qGeV: self.c7p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[8] = lambda qGeV: self.c8p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[9] = lambda qGeV: self.c9p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[10] = lambda qGeV: self.c10p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[11] = lambda qGeV: self.c11p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.ProtonBasis[12] = lambda qGeV: self.c12p((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[1] = lambda qGeV: self.c1n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[4] = lambda qGeV: self.c4n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[5] = lambda qGeV: self.c5n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[6] = lambda qGeV: self.c6n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[7] = lambda qGeV: self.c7n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[8] = lambda qGeV: self.c8n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[9] = lambda qGeV: self.c9n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[10] = lambda qGeV: self.c10n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[11] = lambda qGeV: self.c11n((qGeV * GeV) ** 2) * 246.2 ** 2
		self.NeutronBasis[12] = lambda qGeV: self.c12n((qGeV * GeV) ** 2) * 246.2 ** 2
