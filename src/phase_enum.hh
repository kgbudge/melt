// phase_enum.hh
//
// Copyright (C) 2019 - Kent G. Budge
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef phase_enum_hh
#define phase_enum_hh

enum PHASE
{
	// basic components
  P_H2,
  P_GRAPHITE,
  P_O2,
  P_Na,
  P_Mg,
  P_Al,
  P_Si,
  P_P4,
  P_S,
  P_Cl2,
  P_K,
  P_Ca,
  P_Ti,
  P_Cr,
  P_Mn,
  P_Fe,
  P_Zr,


	  // alternate components
  P_H2O_LIQUID,
  P_CO2,
  P_QUARTZ,
  P_RUTILE,
  P_CORUNDUM,
  P_HEMATITE,
  P_WUSTITE,
  P_WUSTITE_LIQUID,
  P_MANGANOSITE,
  P_PERICLASE,
  P_LIME,
  P_Na2O,
  P_K2O,
  P_ESKOLAITE,
  P_BADDELEYITE,
  P_HALITE,

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
  P_COESITE,
  P_CORDIERITE,
  P_FE_CORDIERITE,
  P_OH_CORDIERITE,
  P_MN_CORDIERITE,
  P_MG_CORUNDUM,
  P_CORUNDUM_LIQUID,
  P_CRISTOBALITE,
  P_CUMMINGTONITE,
  P_K_CYMRITE,
  P_DAPHNITE,
  P_DIASPORE,
  P_DEERITE,
  P_DIAMOND,
  P_DIOPSIDE,
  P_DIOPSIDE_LIQUID,
  P_DISULFUR,
  P_DOLOMITE,
  P_EASTONITE,
  P_ENSTATITE,
  P_ENSTATITE_LIQUID,
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
  P_GREENALITE,
  P_GROSSULAR,
  P_GRUNERITE,
  P_H2S,
  P_HALITE_LIQUID,
  P_HEDENBERGITE,
  P_HERCYNITE,
  P_HEULANDITE,
  P_HOLLANDITE,
  P_ILMENITE,
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
  P_FE_PREHNITE,
  P_PREHNITE,
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
  P_STISHOVITE,
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
  P_TRIDYMITE,
  P_TROILITE,
  P_LOW_TROILITE,
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

  P_SiO2_LIQUID,
  P_Na2SiO3_LIQUID,
  P_KAlSiO4_LIQUID,

  P_Na2SiO3,

  P_HYDROXYAPATITE,

  P_END
};

#endif // phase_enum_hh
