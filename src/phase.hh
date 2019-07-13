#ifndef phase_hh
#define phase_hh

#include <vector>

double const P0 = 1.0e-3; // in kbar
double const T0 = 298.15; // K
double const R = 0.0083144598;  // kJ/K
double const V0 = R*T0/P0;  // kJ/kbar = 10 cm^3

enum ELEMENT
{
  E_H,
  E_C,
  E_O,
  E_NA,
  E_MG,
  E_AL,
  E_SI,
  E_S,
  E_CL,
  E_K,
  E_CA,
  E_TI,
  E_CR,
  E_MN,
  E_FE,
  E_ZR,

  E_END
};

unsigned const MAX_Z = 7;

enum PHASE
{
	  // components
  P_H2O_LIQUID,
  P_CO2,
  P_SiO2_QUARTZ,
  P_TiO2,
  P_Al2O3,
  P_Fe2O3,
  P_FeO,
  P_MnO,
  P_MgO,
  P_CaO,
  P_Na2O,
  P_K2O,
  P_Cr2O3,
  P_ZrO2,
  P_HALITE,
  P_S,

      // alternate components
  P_O2,

      // pure phases
  P_ACMITE,
  P_AKERMANITE,
  P_AKIMOTOITE,
  P_FE_AKIMOTOITE,
  P_ALBITE,
  P_ALBITE_HIGH,
  P_ALBITE_LIQUID,
  P_ALMANDINE,
  P_AMESITE,
  P_ANALCITE,
  P_ANDALUSITE,
  P_ANDRADITE,
  P_ANHYDRITE,
  P_ANKERITE,
  P_ANNITE,
  P_ANORTHITE,
  P_ANORTHITE_LIQUID,
  P_ANTHOPHYLLITE,
  P_FE_ANTHOPHYLLITE,
  P_ANTIGORITE,
  P_ARAGONITE,
  P_MN_BIOTITE,
  P_BIXBYITE,
  P_BRUCITE,
  P_CALCITE,
  P_CARNEGIEITE,
  P_CARNEGIEITE_HI,
  P_CELADONITE,
  P_MN_CHLORITE,
  P_NOAL_CHLORITE,
  P_MG_CHLORITOID,
  P_MN_CHLORITOID,
  P_CHRYSOTILE,
  P_CLINOCHLORE,
  P_CLINOENSTATITE,
  P_CLINOENSTATITE_HP,
  P_CLINOHUMITE,
  P_CLINOZOISITE,
  P_CO,
  P_SiO2_COESITE,
  P_CORDIERITE,
  P_FE_CORDIERITE,
  P_OH_CORDIERITE,
  P_MN_CORDIERITE,
  P_MG_CORUNDUM,
  P_CORUNDUM_LIQUID,
  P_SiO2_CRISTOBALITE,
  P_CUMMINGTONITE,
  P_K_CYMRITE,
  P_DAPHNITE,
  P_DEERITE,
  P_DIAMOND,
  P_DIASPORE,
  P_DIOPSIDE,
  P_DIOPSIDE_LIQUID,
  P_DISULFUR,
  P_DOLOMITE,
  P_ENSTATITE,
  P_ENSTATITE_LIQUID,
  P_EASTONITE,
  P_EPIDOTE,
  P_FE_EPIDOTE,
  P_FAYALITE,
  P_FAYALITE_LIQUID,
  P_FE_CHLORITOID,
  P_FE_PUMPELLYITE,
  P_FE_STAUROLITE,
  P_FERROACTINOLITE,
  P_FERROCARPHOLITE,
  P_FERROCELADONITE,
  P_FERROSILITE,
  P_FERROSUDOITE,
  P_FORSTERITE,
  P_FORSTERITE_LIQUID,
  P_GEHLENITE,
  P_GEIKIELITE,
  P_FE_GLAUCOPHANE,
  P_MG_GLAUCOPHANE,
  P_GOETHITE,
  P_GRAPHITE,
  P_GREENALITE,
  P_GROSSULAR,
  P_GRUNERITE,
  P_H2S,
  P_HALITE_LIQUID,
  P_HEDENBERGITE,
  P_HERCYNITE,
  P_HEULANDITE,
  P_HOLLANDITE,
  P_HYDROGEN,
  P_ILMENITE,
  P_IRON,
  P_JADEITE,
  P_JULGOLDITE,
  P_K_FELDSPAR_LIQUID,
  P_KALSILITE,
  P_KAOLINITE,
  P_KNORRINGITE,
  P_KOSMOCHLORE,
  P_KYANITE,
  P_LARNITE,
  P_LAUMONTITE,
  P_LAWSONITE,
  P_LEUCITE,
  P_LEUCITE_LIQUID,
  P_LIME_LIQUID,
  P_LIZARDITE,
  P_MAGNESIOCARPHOLITE,
  P_MAGNESIOFERRITE,
  P_MAGNESITE,
  P_MAGNETITE,
  P_MAJORITE,
  P_MARGARITE,
  P_MEIONITE,
  P_MERWINITE,
  P_METHANE,
  P_MG_PUMPELLYITE,
  P_MG_STAUROLITE,
  P_MICROCLINE,
  P_MINNESOTAITE,
  P_MG_MINNESOTAITE,
  P_MONTICELLITE,
  P_AL_MULLITE,
  P_MUSCOVITE,
  P_NEPHELINE,
  P_NEPHELINE_LIQUID,
  P_FE_OSUMILITE,  
  P_OSUMILITE_1,
  P_OSUMILITE_2,
  P_PARAGONITE,
  P_PARGASITE,
  P_PERICLASE_LIQUID,
  P_AL_PEROVSKITE,
  P_CA_PEROVSKITE,
  P_FE_PEROVSKITE,
  P_MG_PEROVSKITE,
  P_PHA,
  P_PHLOGOPITE,
  P_NA_PHLOGOPITE,
  P_PICROCHROMITE,
  P_PIEMONTITE,
  P_PREHNITE,
  P_FE_PREHNITE,
  P_PROTOENSTATITE,
  P_PSEUDOWOLLASTONITE,
  P_PYRITE,
  P_PYROPE,
  P_PYROPHANITE,
  P_PYROPHYLLITE,
  P_CA_ESKOLA_PYROXENE,
  P_CA_TSCHERMAK_PYROXENE,
  P_MG_TSCHERMAK_PYROXENE,
  P_PYROXMANGITE,
  P_PYRRHOTITE_STOICHOMETRIC,
  P_PYRRHOTITE,
  P_RANKINITE,
  P_RHODOCHROSITE,
  P_RHODONITE,
  P_RIEBECKITE,
  P_FE_RINGWOODITE,
  P_MG_RINGWOODITE,
  P_FE_WADSLEYITE,
  P_MG_WADSLEYITE,
  P_SANIDINE,
  P_FE_SAPPHIRINE,
  P_SAPPHIRINE_4,
  P_SAPPHIRINE_5,
  P_SI_MULLITE,
  P_SIDERITE,
  P_SILLIMANITE,
  P_SILLIMANITE_LIQUID,
  P_SODALITE,
  P_SPESSARTINE,
  P_SPHENE,
  P_SPINEL,
  P_SPURRITE,
  P_SiO2_SHTISHOVITE,
  P_MN_STAUROLITE,
  P_STILBITE,
  P_STILPNOMELANE,
  P_MG_STILPNOMELANE,
  P_SUDOITE,
  P_SYLVITE,
  P_SYLVITE_LIQUID,
  P_TALC,
  P_FE_TALC,
  P_PRL_TALC,
  P_TSCHERMAK_TALC,
  P_TEPHROITE,
  P_TILLEYITE,
  P_CASI_TITANITE,
  P_TOPAZ,
  P_TREMOLITE,
  P_SiO2_TRIDYMITE,
  P_TROILITE,
  P_LO_TROLITE,
  P_TSCHERMAKITE,
  P_ULVOSPINEL,
  P_VESUVIANITE,
  P_WADEITE,
  P_WAIRAKITE,
  P_WALSTROMITE,
  P_WATER_VAPOR,
  P_WOLLASTONITE,
  P_WOLLASTONITE_LIQUID,
  P_ZIRCON,
  P_ZOISITE,

  P_SiO2_MELT,

  P_END
};

struct Phase;

class Model
{
  // Temperature in K, pressure in kbar, Gibbs free energy in kJ/mol
	
  // Note that Gf is the apparent free energy of formation, which ignores the 
  // entropy of the elements in their standard state.

  public:
    virtual double Gf(Phase const &phase, double T, double P) const = 0;
    virtual double volume(Phase const &phase, double T, double P) const = 0;
    virtual double P(Phase const &phase, double T, double Gf) const;
};

extern const class Solid : public Model
{
  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;
} *const SOLID;

extern const class Vapor: public Model
{
  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;
    virtual double P(Phase const &phase, double T, double Gf) const;
} *const VAPOR;

extern const class Melt : public Model
{
  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;
} *const MELT;

extern const class Aqueous : public Model
{
  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;
}  *const AQUEOUS;

extern const class Water : public Vapor
{
  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;

  private:

    static double constexpr c13 = 0.24657688e6;
    static double constexpr c14 = 0.51359951e2;
    static double constexpr c23 = 0.58638965;
    static double constexpr c24 = -0.28646939e-2;
    static double constexpr c25 = 0.31375577e-4;
    static double constexpr c33 = -0.62783840e1;
    static double constexpr c34 = 0.14791599e-1;
    static double constexpr c35 = 0.35779579e-3;
    static double constexpr c36 = 0.15432925e-7;
    static double constexpr c44 = -0.42719875;
    static double constexpr c45 = -0.16325155e-4;
    static double constexpr c53 = 0.56654978e4;
    static double constexpr c54 = -0.16580167e2;
    static double constexpr c55 = 0.76560762e-1;
    static double constexpr c64 = 0.10917883;
    static double constexpr c71 = 0.38878656e13;
    static double constexpr c72 = -0.13494878e9;
    static double constexpr c73 = 0.30916564e6;
    static double constexpr c74 = 0.75591105e1;
    static double constexpr c83 = -0.65537898e5;
    static double constexpr c84 = 0.18810675e3;
    static double constexpr c91 = -0.14182435e14;
    static double constexpr c92 = 0.18165390e9;
    static double constexpr c93 = -0.19769068e6;
    static double constexpr c94 = -0.23530318e2;
    static double constexpr c03 = 0.92093375e5;
    static double constexpr c04 = 0.12246777e3;

    static double p(double const rho);
    static double rho(Phase const &phase, double T, double P);
    static double rho(double P);

// cached
    static double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, T;
}  *const WATER;

extern const class Magma : public Model
{
  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;
}  *const MAGMA;

extern struct Phase
{
	char const *name;
	unsigned nz; // elements in formula
	unsigned z[MAX_Z]; // elements of formula
	double n[MAX_Z];  // quantities of each element in formula. double because these can be fractional for mineraloids or solid solutions.
   
	double W; // atomic weight
	double V; // molar volume at STP in kJ/kbar = 12.342 cm^3
	double Hf0;  // standard Gibbs free energy of formation at STP in kJ
	double S0; // entropy at STP in J/K
	
	// Cv expansion. Cp in kJ/K
	double a; 
	double b;
	double c;
	double d;

    // Compressibility
    double a0;   // in 1/K
    double k0;   // in kbar
    double k0p;
    double k0pp;

	Model const *model;
}
const phase[P_END];

extern double S0e[E_END];

void initialize_Gf(double T, double P, std::vector<Phase> const &phase, std::vector<double> &Gf);

template<typename Function>
void solve(double const y, double &x, Function);

#endif // phase_hh