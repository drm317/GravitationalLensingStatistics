#!/usr/bin/env python3
"""
Cosmological Parameter Constraints from Gravitational Lensing Statistics

This program calculates constraints on the matter density parameter (Omega_m) 
and the cosmological constant (Omega_Lambda) using gravitational lensing 
statistics with a singular isothermal sphere (SIS) lens model.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

class CosmologyCalculator:
    """Handles cosmological calculations and distance measures."""
    
    def __init__(self, omega_m=0.3, omega_lambda=0.7, h=0.7):
        self.omega_m = omega_m
        self.omega_lambda = omega_lambda
        self.omega_k = 1.0 - omega_m - omega_lambda
        self.h = h
        self.H0 = 100 * h  # km/s/Mpc
        self.c = 299792.458  # km/s
    
    def E(self, z):
        """Dimensionless Hubble parameter E(z) = H(z)/H0"""
        return np.sqrt(self.omega_m * (1 + z)**3 + 
                      self.omega_k * (1 + z)**2 + 
                      self.omega_lambda)
    
    def comoving_distance(self, z, n_points=100):
        """Comoving distance in Mpc"""
        z_array = np.linspace(0, z, n_points)
        dz = z / (n_points - 1)
        integrand = 1.0 / self.E(z_array)
        integral = np.trapz(integrand, dx=dz)
        return (self.c / self.H0) * integral
    
    def angular_diameter_distance(self, z):
        """Angular diameter distance in Mpc"""
        dc = self.comoving_distance(z)
        if self.omega_k == 0:
            return dc / (1 + z)
        elif self.omega_k > 0:
            sqrt_ok = np.sqrt(self.omega_k)
            dh = self.c / self.H0
            return dh * np.sinh(sqrt_ok * dc / dh) / (sqrt_ok * (1 + z))
        else:
            sqrt_ok = np.sqrt(-self.omega_k)
            dh = self.c / self.H0
            return dh * np.sin(sqrt_ok * dc / dh) / (sqrt_ok * (1 + z))

class SISLensModel:
    """Singular Isothermal Sphere lens model implementation."""
    
    def __init__(self, sigma_v, z_lens, z_source, cosmology):
        self.sigma_v = sigma_v  # velocity dispersion in km/s
        self.z_lens = z_lens
        self.z_source = z_source
        self.cosmology = cosmology
        self.einstein_radius = self._calculate_einstein_radius()
    
    def _calculate_einstein_radius(self):
        """Calculate Einstein radius in arcseconds"""
        c_km_s = 299792.458
        
        D_s = self.cosmology.angular_diameter_distance(self.z_source)
        
        # Simplified calculation for D_ls/D_s ratio
        if self.z_lens >= self.z_source:
            return 0
        
        # Approximate D_ls/D_s for the given redshifts
        D_ls_over_D_s = (self.z_source - self.z_lens) / self.z_source
        
        theta_E = 4 * np.pi * (self.sigma_v / c_km_s)**2 * D_ls_over_D_s
        return theta_E * 206265  # convert to arcseconds
    
    def optical_depth(self, theta_threshold):
        """Calculate optical depth for lensing above threshold"""
        if self.einstein_radius == 0:
            return 0
        return (self.einstein_radius / theta_threshold)**2 if theta_threshold > 0 else np.inf

class LensSurvey:
    """Simulates a gravitational lens survey."""
    
    def __init__(self, n_sources=2300, n_lenses=4):
        self.n_sources = n_sources
        self.n_lenses = n_lenses
        self.sources = []
        self.lenses = []
        self._generate_survey()
    
    def _generate_survey(self):
        """Generate simulated survey data"""
        np.random.seed(42)
        
        # Generate source redshifts (median ~1.0)
        self.source_redshifts = np.random.gamma(2, 0.5, self.n_sources)
        self.source_redshifts = np.clip(self.source_redshifts, 0.1, 5.0)
        
        # Generate lens properties
        self.lens_redshifts = np.random.uniform(0.1, 1.5, self.n_lenses)
        self.lens_sigma_v = np.random.normal(250, 50, self.n_lenses)  # km/s
        self.lens_sigma_v = np.clip(self.lens_sigma_v, 150, 400)
        
        # Detection threshold (minimum Einstein radius in arcseconds)
        self.detection_threshold = 1.0
    
    def calculate_lensing_statistics(self, cosmology):
        """Calculate expected and observed lensing statistics"""
        stats = {
            'n_expected': 0,
            'n_observed': self.n_lenses,
            'einstein_radii': [],
            'optical_depths': []
        }
        
        # Calculate expected number of lenses
        for i in range(self.n_sources):
            z_s = self.source_redshifts[i]
            
            # Simplified expected number calculation
            # Typical lensing probability per source
            if z_s > 0.5:  # Only consider sources at reasonable redshift
                lens_prob = 0.001 * (z_s / 1.0)**2  # Simple scaling
                stats['n_expected'] += lens_prob
        
        # Calculate Einstein radii for observed lenses
        for i in range(self.n_lenses):
            z_s = np.random.choice(self.source_redshifts)
            lens = SISLensModel(self.lens_sigma_v[i], self.lens_redshifts[i], z_s, cosmology)
            stats['einstein_radii'].append(lens.einstein_radius)
            stats['optical_depths'].append(lens.optical_depth(self.detection_threshold))
        
        return stats
    
    def _comoving_volume_element(self, z, cosmology):
        """Comoving volume element per unit redshift and solid angle"""
        c = 299792.458
        D_A = cosmology.angular_diameter_distance(z)
        D_H = c / cosmology.H0
        E_z = cosmology.E(z)
        return D_H * (1 + z)**2 * D_A**2 / E_z

class ParameterConstraints:
    """Calculate constraints on cosmological parameters."""
    
    def __init__(self, survey):
        self.survey = survey
    
    def chi_squared(self, params):
        """Calculate chi-squared for given parameters"""
        omega_m, omega_lambda = params
        
        if omega_m < 0 or omega_lambda < 0 or omega_m + omega_lambda > 1.2:
            return 1e10
        
        cosmology = CosmologyCalculator(omega_m, omega_lambda)
        stats = self.survey.calculate_lensing_statistics(cosmology)
        
        # Simple chi-squared based on expected vs observed lens counts
        expected = max(stats['n_expected'], 0.1)  # Avoid division by zero
        observed = stats['n_observed']
        
        # Add Poisson uncertainty
        sigma_expected = np.sqrt(expected)
        chi2 = ((observed - expected) / sigma_expected)**2
        
        # Add constraint from Einstein radius distribution
        if len(stats['einstein_radii']) > 0:
            mean_einstein = np.mean(stats['einstein_radii'])
            expected_mean = 2.0  # Expected mean Einstein radius
            sigma_mean = 0.5
            chi2 += ((mean_einstein - expected_mean) / sigma_mean)**2
        
        return chi2
    
    def find_best_fit(self):
        """Find best-fit parameters using minimization"""
        initial_guess = [0.3, 0.7]
        
        result = minimize(self.chi_squared, initial_guess, 
                         bounds=[(0.05, 0.95), (0.05, 0.95)],
                         method='L-BFGS-B')
        
        return result
    

def main():
    """Main execution function"""
    print("Cosmological Parameter Constraints from Gravitational Lensing")
    print("=" * 60)
    
    # Initialize survey
    print("\n1. Initializing simulated lens survey...")
    survey = LensSurvey(n_sources=2300, n_lenses=4)
    print(f"   - {survey.n_sources} sources with median z = {np.median(survey.source_redshifts):.2f}")
    print(f"   - {survey.n_lenses} lenses with <σ_v> = {np.mean(survey.lens_sigma_v):.0f} km/s")
    
    # Test with fiducial cosmology
    print("\n2. Testing with fiducial ΛCDM cosmology (Ω_m=0.3, Ω_Λ=0.7)...")
    fiducial_cosmo = CosmologyCalculator(0.3, 0.7)
    fiducial_stats = survey.calculate_lensing_statistics(fiducial_cosmo)
    print(f"   - Expected number of lenses: {fiducial_stats['n_expected']:.2f}")
    print(f"   - Observed number of lenses: {fiducial_stats['n_observed']}")
    print(f"   - Mean Einstein radius: {np.mean(fiducial_stats['einstein_radii']):.2f}\"")
    
    # Find best-fit parameters
    print("\n3. Finding best-fit cosmological parameters...")
    constraints = ParameterConstraints(survey)
    best_fit = constraints.find_best_fit()
    
    print(f"   - Best-fit Ω_m = {best_fit.x[0]:.3f}")
    print(f"   - Best-fit Ω_Λ = {best_fit.x[1]:.3f}")
    print(f"   - Minimum χ² = {best_fit.fun:.2f}")
    print(f"   - Optimization success: {best_fit.success}")
    
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY:")
    print(f"Best-fit parameters: Ω_m = {best_fit.x[0]:.3f}, Ω_Λ = {best_fit.x[1]:.3f}")
    print(f"Flat universe parameter: Ω_k = {1 - best_fit.x[0] - best_fit.x[1]:.3f}")

if __name__ == "__main__":
    main()