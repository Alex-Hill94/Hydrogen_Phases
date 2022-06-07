import numpy as np
from decimal import *

getcontext().prec = 40


class Compute_H_Fractions():

	def __init__(self,
				Redshift = 0,
				H_abun = [],
				Temp = [],
				Density = [],
				SFR = [],
				Mass = None,
				PIDs = None):
		
		self.Redshift = Redshift 	# Redshift of the computation
		self.H_abun = H_abun   		# Fraction of gas particle which is Hydrogen
		self.Temp = Temp      		# Kelvin
		self.Density = Density   	# 10^10 Msun  / Mpc^3 
		self.SFR = SFR        		# Solar masses per year
		self.Mass = Mass      		# 1e10 solar masses
		self.PIDs = PIDs      		# IDs of each gas particle

	def compute(self):
		
		###Best-fit paramters from Rahmati et al. 2013 (their Fig. A1) https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2427R/exportcitation

		z_array = np.array((0.0,1.0,2.0,3.0,4.0,5.0))
		uvb     = np.array((8.34e-14,   7.39e-13,    1.50e-12,
					1.16e-12,   7.92e-13,    5.43e-13))
		alph_1  = np.array((-3.98,   -2.94,   -2.22,    -1.99,   -2.05,    -2.63))
		alph_2  = np.array((-1.09,   -0.90,   -1.09,    -0.88,   -0.75,    -0.57))
		beta    = np.array(( 1.29,    1.21,    1.75,     1.72,    1.93,     1.77))
		f_array  = np.array(( 0.01,    0.03,    0.03,     0.04,    0.02,     0.01))
		log_n0  = np.array((-2.94,   -2.29,   -2.06,    -2.13,   -2.23,    -2.35))

		#Conversion factors
		mu       = 1.22
		m_H      = 1.67e-24 # Mass of hydrogen atom in grammes
		msun     = 1e10 * 1.989e33 # Densities are stored in units of 10^10 Msun (g) / Mpc^3 (cm^3)
		mpc_cube = (3.086e24)**3.0
		conv     = msun/mpc_cube

		Z = self.Redshift
		H_abun = self.H_abun
		temp = self.Temp
		density = self.Density
		SFR = self.SFR
		M = self.Mass
		PIDs = self.PIDs		

		#density = np.array(density)
		#M = np.array(M)
		#temp = np.array(temp)
		#H_abun = np.array(H_abun)

		density = density * conv # converts density into g/cm3
		n_H = (H_abun*density)/(m_H) # Number density of hydrogen n = Xp/m

		def interpolate_constants(z):
			try:
				if z<=1:
					offset = 0
				elif z<=2:
					offset = 1
				elif z<=3:
					offset = 2
				elif z<=4:
					offset = 3
				elif z<=5:
					offset = 4

				idx_low = offset
				idx_high = offset+1

				z_l = z_array[idx_low]
				z_h = z_array[idx_high]
				uv_l = uvb[idx_low]
				uv_h = uvb[idx_high]
				a1_l = alph_1[idx_low]
				a1_h = alph_1[idx_high]
				a2_l = alph_2[idx_low]
				a2_h = alph_2[idx_high]
				bt_l = beta[idx_low]
				bt_h = beta[idx_high]
				f_l = f_array[idx_low]
				f_h = f_array[idx_high]
				ln0_l = log_n0[idx_low]
				ln0_h = log_n0[idx_high]

				dz  = (z-z_l)/(z_h-z_l)
				uv  = uv_l+(uv_h-uv_l)*dz
				a1  = a1_l+(a1_h-a1_l)*dz
				a2  = a2_l+(a2_h-a2_l)*dz
				bt  = bt_l+(bt_h-bt_l)*dz
				F   = f_l+(f_h-f_l)*dz
				ln0 = ln0_l+(ln0_h-ln0_l)*dz
				return dz, uv, a1, a2, bt, F, ln0
			except:
				print('Invalid redshift')


		dz, uv, a1, a2, bt, f, ln0 = interpolate_constants(Z)

		#### Computation of the neutral fraction as in Rahmati et al. 2013####

		def eta(N):

			def alp_A(T):
				lambd = 315614.0/T
				ans    = (1.269e-13)*lambd**1.503
				ans    = ans/(1.0+(lambd/0.522)**0.47)**1.923
				return ans

			def LAMBDA(T): 
				ans = (1.17e-10)*np.sqrt(T)*np.exp(-157809.0/T)
				ans = ans/(1.0+np.sqrt(T/100000.0))
				return ans

			def gam_phot(N):
				n0  = 10.0**ln0
				ans = (1.0-f)
				ans = ans*(1.0 + (N/n0)**bt)**a1
				ans = ans + f*(1.0 + N/n0)**a2
				ans = ans
				return ans*uv

			aA     = alp_A(temp)
			lam    = LAMBDA(temp)
			phot   = gam_phot(n_H)

			A = np.array([aA   + lam], dtype=np.dtype(Decimal))
			B = np.array([2.0*aA + (phot/N) + lam], dtype=np.dtype(Decimal))
			C = np.array([aA], dtype=np.dtype(Decimal))
			D = np.array([B*B - 4.0*A*C], dtype=np.dtype(Decimal))
			ans = B - (D)**0.5
			ans = ans/(2.0*A)
			return ans[0][0].astype(float)

		f_neutral =  eta(n_H)

		### Computation of H2 fraction via Blitz and Rosowolski 2006 ###

		a = 1/(1+Z)
		P = (n_H*temp)/(mu*H_abun) # Analytic recomputation of pressure is needed due to artificial temperature floor
		P0    = 4.3e4 
		alpha = 0.92
		R_mol = (P/P0)**alpha
		h2_fraction = 1.0/(1.0+(1.0/R_mol))
		h2_fraction[SFR < 10**-10] = 0.0
		HI_fraction  = f_neutral - h2_fraction
		HII_fraction = 1 - HI_fraction

		self.n_H = n_H
		self.Neutral_Fraction = f_neutral
		self.H2_Fraction = h2_fraction
		self.HI_Fraction = HI_fraction
		self.HII_fraction = HII_fraction
		
		if self.Mass is not None:
			self.Neutral_Mass = self.Neutral_Fraction * self.Mass
			self.H2_Mass = self.H2_Fraction * self.Mass
			self.HI_Mass = self.HI_Fraction * self.Mass
			self.HII_Mass = self.HII_fraction * self.Mass

	def plot_example(self):

		n_H  = np.log10(self.n_H).astype('float64')
		temp = np.log10(self.Temp).astype('float64')
		n_edges = np.linspace(np.nanmin(n_H), np.nanmax(n_H), 512) # Set edges for the bins
		t_edges = np.linspace(np.nanmin(temp), np.nanmax(temp), 512)

		def dense(n, t, W): ### Function that creates 2D prob dist of n and T data
			hist, dum1, dum2 = np.histogram2d(n, t, bins = [n_edges, t_edges], weights = W, normed = True) # 2D histogram
			hist = matrix.transpose(hist) # Sets number density to be the x-axis and temperature to be the y-axis
			hist = np.ma.masked_where(hist == 0, hist) # Set 0s to be 'bad data'
			hist = np.log10(hist) # Moves histogram into log space
			return hist



		hist_all = dense(n_H, temp, self.HI_Mass)

		plt.figure()
		cmap = plt.cm.rainbow
		cmap.set_bad(color='white')
		#plt.contour(hist1)
		plt.imshow(hist_all, cmap = 'jet', origin = "lower", interpolation = 'nearest', extent = [np.min(n_H), np.max(n_H),np.min(temp), np.max(temp)], aspect = 'auto')
		plt.colorbar(label = '$\\mathrm{log_{10}\ (HI\ Mass\ Fraction)}$')
		plt.ylabel('$\\mathrm{log_{10}\ T\ [K]}$')
		plt.xlabel('$\\mathrm{log_{10}\ n_{H}\ [cm^{-3}]}$')
		plt.savefig('Test_Image.png')
