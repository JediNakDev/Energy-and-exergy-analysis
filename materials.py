"""
Material and PCM classes with predefined instances.

Layer Structure (top to bottom):
    g: Tempered glass
    b: EVA
    c: PV cells
    d: EVA
    e: Tedlar foil
    p: PCM (optional)
    f: Transparency acrylic glass
"""

from constants import PANEL_AREA, PCM_THICKNESS


# =============================================================================
# Material Class
# =============================================================================
class Material:
    """Material properties for thermal simulation layers."""
    
    def __init__(self, c_p=None, d=None, emi=None, p=None, k=None):
        """
        Initialize material with optional properties.
        
        Args:
            c_p: Specific heat capacity [J/(kg·K)]
            d: Thickness [m]
            emi: Emissivity
            p: Density [kg/m³]
            k: Thermal conductivity [W/(m·K)]
        """
        if all(v is not None for v in [c_p, d, emi, p, k]):
            self.properties(c_p, d, emi, p, k)
    
    def properties(self, c_p, d, emi, p, k):
        """Set material properties."""
        self.c_p = c_p
        self.d = d
        self.emi = emi
        self.m = p * PANEL_AREA * d
        self.k = k
        self.M = self.m * c_p
        self.kperd = k / d


class PCM:
    """Phase Change Material properties."""
    
    def __init__(self, p=None, c_ps=None, c_pl=None, L=None, emi=None, Tm=None, price=None):
        """
        Initialize PCM with optional properties.
        
        Args:
            p: Density [kg/m³]
            c_ps: Specific heat capacity (solid) [J/(kg·K)]
            c_pl: Specific heat capacity (liquid) [J/(kg·K)]
            L: Latent heat of fusion [J/kg]
            emi: Emissivity
            Tm: Melting temperature [K]
            price: Cost [$/kg] (optional)
        """
        if all(v is not None for v in [p, c_ps, c_pl, L, emi, Tm]):
            self.properties(p, c_ps, c_pl, L, emi, Tm, price)
    
    def properties(self, p, c_ps, c_pl, L, emi, Tm, price=None):
        """Set PCM properties."""
        self.m = p * PANEL_AREA * PCM_THICKNESS
        self.c_ps = c_ps
        self.c_pl = c_pl
        self.L = L
        self.Tm = Tm
        self.M_s = self.m * c_ps
        self.M_sl = self.m * L
        self.M_l = self.m * c_pl
        self.emi = emi
        self.price = price


# =============================================================================
# Predefined Layer Materials
# =============================================================================
# Layer g: Tempered glass
tempered_glass = Material(c_p=840, d=0.0032, emi=0.85, p=2520, k=0.96)

# Layer b: EVA (front)
eva_front = Material(c_p=3135, d=0.0005, emi=0.53, p=934, k=0.24)

# Layer c: PV cells
pv_cells = Material(c_p=710.08, d=0.0003, emi=0.75, p=2330, k=148)

# Layer d: EVA (back)
eva_back = Material(c_p=3135, d=0.0005, emi=0.53, p=934, k=0.24)

# Layer e: Tedlar foil
tedlar_foil = Material(c_p=1050, d=0.0005, emi=0.08, p=2700, k=0.1583)

# Layer f: Transparency acrylic glass
acrylic_glass = Material(c_p=840, d=0.005, emi=0.85, p=2500, k=0.96)


# =============================================================================
# Predefined PCM Materials
# =============================================================================
PCM_CATALOG = {
    "paraffin_wax": {
        "p": 800, "c_ps": 2500, "c_pl": 2500, "L": 140000, 
        "Tm": 296, "emi": 0.91, "price": 95.5696
    },
    "sodium_sulfate_decahydrate": {
        "p": 1460, "c_ps": 1631.65, "c_pl": 919.16, "L": 61630, 
        "Tm": 305.55, "emi": 0.96, "price": 9.87252
    },
    "polyethylene_glycol": {
        "p": 1127, "c_ps": 2142, "c_pl": 2490, "L": 159000, 
        "Tm": 310.65, "emi": 0.1, "price": 86.368772
    },
    "calcium_chloride_hexahydrate": {
        "p": 1470, "c_ps": 1400, "c_pl": 2100, "L": 140000, 
        "Tm": 297.15, "emi": 0.95, "price": 854.023695
    },
    "octadecane": {
        "p": 778.6, "c_ps": 1910, "c_pl": 2230, "L": 240000, 
        "Tm": 301.15, "emi": 0.954, "price": 7019.8576
    },
    "lauric_acid": {
        "p": 883, "c_ps": 2200, "c_pl": 1950, "L": 187200, 
        "Tm": 298.6, "emi": 0.903, "price": 118.421779
    },
    "capric_acid": {
        "p": 890, "c_ps": 2760, "c_pl": 2760, "L": 163370, 
        "Tm": 304.15, "emi": 0.92, "price": 1034.625445
    },
}


def create_pcm(name: str) -> PCM:
    """Create a PCM instance from the catalog."""
    if name not in PCM_CATALOG:
        raise ValueError(f"Unknown PCM: {name}. Available: {list(PCM_CATALOG.keys())}")
    props = PCM_CATALOG[name]
    pcm = PCM()
    pcm.properties(**props)
    return pcm


# Predefined PCM instances
paraffin_wax = create_pcm("paraffin_wax")
sodium_sulfate = create_pcm("sodium_sulfate_decahydrate")
