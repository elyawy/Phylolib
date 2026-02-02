// PrecomputedChebyshevData.h
// Precomputed Chebyshev coefficients for all standard amino acid substitution models
// This eliminates the 30-60ms initialization bottleneck

#ifndef PRECOMPUTED_CHEBYSHEV_DATA_H
#define PRECOMPUTED_CHEBYSHEV_DATA_H

#include "definitions.h" // for VVVdouble

namespace PrecomputedCheby {
    // cpREV45
    extern const VVVdouble cpREV45_coff;
    extern const VVVdouble cpREV45_derv_coff;
    extern const VVVdouble cpREV45_sec_derv_coff;
    
    // dayhoff
    extern const VVVdouble dayhoff_coff;
    extern const VVVdouble dayhoff_derv_coff;
    extern const VVVdouble dayhoff_sec_derv_coff;
    
    // jones (JTT)
    extern const VVVdouble jones_coff;
    extern const VVVdouble jones_derv_coff;
    extern const VVVdouble jones_sec_derv_coff;
    
    // mtREV24
    extern const VVVdouble mtREV24_coff;
    extern const VVVdouble mtREV24_derv_coff;
    extern const VVVdouble mtREV24_sec_derv_coff;
    
    // WAG
    extern const VVVdouble wag_coff;
    extern const VVVdouble wag_derv_coff;
    extern const VVVdouble wag_sec_derv_coff;
    
    // HIVb
    extern const VVVdouble HIVb_coff;
    extern const VVVdouble HIVb_derv_coff;
    extern const VVVdouble HIVb_sec_derv_coff;
    
    // HIVw
    extern const VVVdouble HIVw_coff;
    extern const VVVdouble HIVw_derv_coff;
    extern const VVVdouble HIVw_sec_derv_coff;
    
    // LG
    extern const VVVdouble lg_coff;
    extern const VVVdouble lg_derv_coff;
    extern const VVVdouble lg_sec_derv_coff;
    
    // NOTE: empiriCodon is skipped (codon model with 61 states, not 20)
    
    // EX_BURIED
    extern const VVVdouble EX_BURIED_coff;
    extern const VVVdouble EX_BURIED_derv_coff;
    extern const VVVdouble EX_BURIED_sec_derv_coff;
    
    // EX_EXPOSED
    extern const VVVdouble EX_EXPOSED_coff;
    extern const VVVdouble EX_EXPOSED_derv_coff;
    extern const VVVdouble EX_EXPOSED_sec_derv_coff;
    
    // EHO_EXTENDED
    extern const VVVdouble EHO_EXTENDED_coff;
    extern const VVVdouble EHO_EXTENDED_derv_coff;
    extern const VVVdouble EHO_EXTENDED_sec_derv_coff;
    
    // EHO_HELIX
    extern const VVVdouble EHO_HELIX_coff;
    extern const VVVdouble EHO_HELIX_derv_coff;
    extern const VVVdouble EHO_HELIX_sec_derv_coff;
    
    // EHO_OTHER
    extern const VVVdouble EHO_OTHER_coff;
    extern const VVVdouble EHO_OTHER_derv_coff;
    extern const VVVdouble EHO_OTHER_sec_derv_coff;
    
    // EX_EHO_BUR_EXT
    extern const VVVdouble EX_EHO_BUR_EXT_coff;
    extern const VVVdouble EX_EHO_BUR_EXT_derv_coff;
    extern const VVVdouble EX_EHO_BUR_EXT_sec_derv_coff;
    
    // EX_EHO_BUR_HEL
    extern const VVVdouble EX_EHO_BUR_HEL_coff;
    extern const VVVdouble EX_EHO_BUR_HEL_derv_coff;
    extern const VVVdouble EX_EHO_BUR_HEL_sec_derv_coff;
    
    // EX_EHO_BUR_OTH
    extern const VVVdouble EX_EHO_BUR_OTH_coff;
    extern const VVVdouble EX_EHO_BUR_OTH_derv_coff;
    extern const VVVdouble EX_EHO_BUR_OTH_sec_derv_coff;
    
    // EX_EHO_EXP_EXT
    extern const VVVdouble EX_EHO_EXP_EXT_coff;
    extern const VVVdouble EX_EHO_EXP_EXT_derv_coff;
    extern const VVVdouble EX_EHO_EXP_EXT_sec_derv_coff;
    
    // EX_EHO_EXP_HEL
    extern const VVVdouble EX_EHO_EXP_HEL_coff;
    extern const VVVdouble EX_EHO_EXP_HEL_derv_coff;
    extern const VVVdouble EX_EHO_EXP_HEL_sec_derv_coff;
    
    // EX_EHO_EXP_OTH
    extern const VVVdouble EX_EHO_EXP_OTH_coff;
    extern const VVVdouble EX_EHO_EXP_OTH_derv_coff;
    extern const VVVdouble EX_EHO_EXP_OTH_sec_derv_coff;
}

#endif // PRECOMPUTED_CHEBYSHEV_DATA_H
