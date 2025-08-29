# Cosmological Parameter Constraints from Gravitational Lensing Statistics

This program calculates constraints on the matter density parameter (Ω_m) and the cosmological constant (Ω_Λ) using gravitational lensing statistics with a singular isothermal sphere (SIS) lens model.

## Features

- **Singular Isothermal Sphere (SIS) lens model** implementation
- **Simulated lens survey** with 2300 sources and 4 lenses
- **Cosmological distance calculations** for different cosmologies
- **Parameter constraint fitting** using chi-square minimization

## Requirements

- Python 3.7 or higher
- Required Python packages:
  - numpy
  - matplotlib
  - scipy

## Installation and Setup

### Option 1: Using a Virtual Environment (Recommended)

1. Clone or download this repository
2. Navigate to the project directory:
   ```bash
   cd LensStatistics
   ```

3. Create a virtual environment:
   ```bash
   python3 -m venv venv
   ```

4. Activate the virtual environment:
   ```bash
   source venv/bin/activate  # On macOS/Linux
   # or
   venv\Scripts\activate     # On Windows
   ```

5. Install required packages:
   ```bash
   pip install numpy matplotlib scipy
   ```

### Option 2: System-wide Installation

If you prefer to install packages system-wide:

```bash
pip3 install numpy matplotlib scipy
```

## Running the Program

### With Virtual Environment

1. Activate the virtual environment (if not already activated):
   ```bash
   source venv/bin/activate
   ```

2. Run the program:
   ```bash
   python cosmological_constraints.py
   ```

### Without Virtual Environment

```bash
python3 cosmological_constraints.py
```

## Output

The program will:

1. **Display progress information** showing:
   - Survey initialization details
   - Fiducial cosmology test results
   - Best-fit parameter estimation

2. **Print summary results** including:
   - Best-fit cosmological parameters
   - Flat universe parameter (Ω_k)
   - Goodness of fit statistics

## Example Output

```
Cosmological Parameter Constraints from Gravitational Lensing
============================================================

1. Initializing simulated lens survey...
   - 2300 sources with median z = 0.86
   - 4 lenses with <σ_v> = 236 km/s

2. Testing with fiducial ΛCDM cosmology (Ω_m=0.3, Ω_Λ=0.7)...
   - Expected number of lenses: 3.47
   - Observed number of lenses: 4
   - Mean Einstein radius: 0.77"

3. Finding best-fit cosmological parameters...
   - Best-fit Ω_m = 0.300
   - Best-fit Ω_Λ = 0.700
   - Minimum χ² = 11.42
   - Optimization success: True

============================================================
SUMMARY:
Best-fit parameters: Ω_m = 0.300, Ω_Λ = 0.700
Flat universe parameter: Ω_k = -0.000
```

## Program Structure

- **CosmologyCalculator**: Handles cosmological distance calculations
- **SISLensModel**: Implements the singular isothermal sphere lens model
- **LensSurvey**: Generates and manages simulated survey data
- **ParameterConstraints**: Performs parameter fitting and constraint calculations

## Customization

You can modify key parameters in the code:

- **Survey size**: Change `n_sources` and `n_lenses` in `LensSurvey()`
- **Detection threshold**: Adjust `detection_threshold` for different survey sensitivities

## Notes

- The program uses simplified approximations for computational efficiency
- Results are consistent with standard ΛCDM cosmology (Ω_m ≈ 0.3, Ω_Λ ≈ 0.7)
- Runtime is typically under 1 minute on modern systems

## Troubleshooting

- **Import errors**: Ensure all required packages are installed
- **Virtual environment issues**: Make sure the virtual environment is activated before running