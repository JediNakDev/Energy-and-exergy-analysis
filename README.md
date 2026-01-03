# Energy and Exergy Analysis of PCM-Integrated Solar Panels

[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1742--6596%2F2653%2F1%2F012038-blue)](https://doi.org/10.1088/1742-6596/2653/1/012038)

Simulation code for the paper: **"Energy and exergy analysis with cost-efficiency comparison of modified phase change material (PCM) integrated with solar panel thermal system"**

Published in: *Journal of Physics: Conference Series*

## Abstract

This study investigates the efficiency and cost-effectiveness of integrating Phase Change Materials (PCMs) into solar panels to optimize power output and power conversion efficiency under varying temperature conditions. By applying PCMs, temperature fluctuations can be reduced, maintaining solar panels at an optimal temperature for power generation.

**Key Findings:**
- Solar panels integrated with PCM generate up to **11.1% more power** than conventional panels
- Best cost-efficiency achieved with paraffin wax and polyethylene glycol combination: **637.7% improvement**

## Methodology

The simulation models thermal dynamics of a multi-layer PV panel structure:

```
┌─────────────────────────────┐  ← Solar Irradiance
│  g: Tempered Glass          │
├─────────────────────────────┤
│  b: EVA                     │
├─────────────────────────────┤
│  c: PV Cells                │  ← Power Generation
├─────────────────────────────┤
│  d: EVA                     │
├─────────────────────────────┤
│  e: Tedlar Foil             │
├─────────────────────────────┤
│  p: Phase Change Material   │  ← Thermal Regulation (optional)
├─────────────────────────────┤
│  f: Transparency Acrylic    │
└─────────────────────────────┘
```

### Numerical Approach

- **Coupled ODEs**: Models heat transfer between layers using energy balance equations
- **Time-stepping scheme**: Custom numerical integration with `dT/dt = (T[i+1] - T[i]) / Δt`
- **Matrix inversion**: Per-step solution of `A·X = B` via `X = A⁻¹·B` to propagate temperature states
- **Phase change handling**: Tracks PCM state (solid/liquid) and latent heat absorption

## Project Structure

```
Energy-and-exergy-analysis/
│
├── Core Modules
│   ├── constants.py        # Physical constants & simulation parameters
│   ├── materials.py        # Material/PCM classes with thermal properties
│   ├── utils.py            # Utility functions (smoothing, thermal conductance)
│   └── power.py            # PV power calculation (single-diode model)
│
├── Data
│   ├── data.py             # Raw measurement data (temperature, irradiance)
│   ├── data_smooth.py      # Data smoothing and preprocessing
│   ├── temperaturedata.py  # Simulation results (temperature)
│   └── powerdata.py        # Simulation results (power)
│
├── Simulation Scripts
│   ├── findtemp_NoPCM.py   # Conventional PV simulation (no PCM)
│   ├── findtemp.py         # Single PCM layer simulation
│   └── findtemp_modifyPCM.py # Dual PCM mixture optimization
│
├── Visualization
│   ├── temperaturegraph.py # Temperature profile plots
│   ├── powergraph.py       # Power output plots
│   ├── effgraph.py         # Efficiency comparison plots
│   └── I-V curve.py        # Current-voltage characteristics
│
├── Reference
│   ├── properties.py       # PCM properties lookup table
│   └── test.py             # Testing utilities
│
└── requirements.txt        # Python dependencies
```

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd Energy-and-exergy-analysis

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Run Temperature Simulation

```bash
# Conventional PV (without PCM)
python findtemp_NoPCM.py

# PV with single PCM
python findtemp.py

# PV with dual PCM optimization
python findtemp_modifyPCM.py
```

### Generate Plots

```bash
python temperaturegraph.py   # Temperature profiles
python powergraph.py         # Power output comparison
python effgraph.py           # Efficiency analysis
python "I-V curve.py"        # I-V characteristics
```

### Configure PCM Type

Edit the PCM configuration in `findtemp.py` or `findtemp_modifyPCM.py`:

```python
from materials import create_pcm

# Available PCMs:
# - "paraffin_wax"
# - "sodium_sulfate_decahydrate"
# - "polyethylene_glycol"
# - "calcium_chloride_hexahydrate"
# - "octadecane"
# - "lauric_acid"
# - "capric_acid"

p = create_pcm("paraffin_wax")
```

## Data Source

Temperature and light intensity data collected in **Wangchan District, Rayong Province, Thailand** using:
- Light intensity sensor
- Thermocouple

## PCM Properties

| Material | Density (kg/m³) | Melting Temp (K) | Latent Heat (kJ/kg) |
|----------|-----------------|------------------|---------------------|
| Paraffin Wax | 800 | 296 | 140 |
| Sodium Sulfate Decahydrate | 1460 | 305.55 | 61.6 |
| Polyethylene Glycol | 1127 | 310.65 | 159 |
| Calcium Chloride Hexahydrate | 1470 | 297.15 | 140 |
| Octadecane | 778.6 | 301.15 | 240 |
| Lauric Acid | 883 | 298.6 | 187.2 |
| Capric Acid | 890 | 304.15 | 163.4 |

## Citation

If you use this code in your research, please cite:

```bibtex
@article{pcm_solar_2023,
  title={Energy and exergy analysis with cost-efficiency comparison of modified phase change material (PCM) integrated with solar panel thermal system},
  journal={Journal of Physics: Conference Series},
  volume={2653},
  number={1},
  pages={012038},
  year={2023},
  doi={10.1088/1742-6596/2653/1/012038}
}
```

## License

This project is for academic and research purposes.

## Contact

For questions regarding the simulation methodology or paper, please refer to the [published paper](https://doi.org/10.1088/1742-6596/2653/1/012038).

