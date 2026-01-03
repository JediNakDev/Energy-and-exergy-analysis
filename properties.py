"""
PCM (Phase Change Material) properties lookup.

Reference data for various PCM materials used in PV-PCM systems.
"""

# =============================================================================
# Imports
# =============================================================================
import chemicals
from chemicals import CAS_from_any, Tm


# =============================================================================
# PCM Properties Reference
# =============================================================================
# Format: (density, c_p_solid, c_p_liquid, latent_heat, melting_temp, emissivity, price)
# Units: (kg/m³, J/(kg·K), J/(kg·K), J/kg, K, -, $/kg)

PCM_PROPERTIES = {
    1: {
        "name": "Paraffin Wax",
        "p": 800, "c_ps": 2500, "c_pl": 2500, "L": 140000, "Tm": 296, "emi": 0.91, "price": 95.5696
    },
    2: {
        "name": "Sodium Sulfate Decahydrate",
        "p": 1460, "c_ps": 1631.65, "c_pl": 919.16, "L": 61630, "Tm": 305.55, "emi": 0.96, "price": 9.87252
    },
    3: {
        "name": "Polyethylene Glycol",
        "p": 1127, "c_ps": 2142, "c_pl": 2490, "L": 159000, "Tm": 310.65, "emi": 0.1, "price": 86.368772
    },
    4: {
        "name": "Calcium Chloride Hexahydrate",
        "p": 1470, "c_ps": 1400, "c_pl": 2100, "L": 140000, "Tm": 297.15, "emi": 0.95, "price": 854.023695
    },
    5: {
        "name": "Octadecane",
        "p": 778.6, "c_ps": 1910, "c_pl": 2230, "L": 240000, "Tm": 301.15, "emi": 0.954, "price": 7019.8576
    },
    6: {
        "name": "Lauric Acid",
        "p": 883, "c_ps": 2200, "c_pl": 1950, "L": 187200, "Tm": 298.6, "emi": 0.903, "price": 118.421779
    },
    7: {
        "name": "Capric Acid",
        "p": 890, "c_ps": 2760, "c_pl": 2760, "L": 163370, "Tm": 304.15, "emi": 0.92, "price": 1034.625445
    },
}


# =============================================================================
# Example Usage
# =============================================================================
if __name__ == "__main__":
    # Lookup example using chemicals library
    CAS_PCM = CAS_from_any('Capric Acid')
    cp_l = chemicals.miscdata.lookup_VDI_tabular_data('7782-41-4', 'Cp (l)')
    print(f"CAS: {CAS_PCM}")
    print(f"Cp (liquid): {cp_l}")
    
    # Print all PCM properties
    print("\nPCM Properties Reference:")
    print("-" * 80)
    for idx, props in PCM_PROPERTIES.items():
        print(f"{idx}. {props['name']}")
        print(f"   ρ={props['p']} kg/m³, Tm={props['Tm']} K, L={props['L']} J/kg")
