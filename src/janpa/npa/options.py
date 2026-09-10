# janpa/npa/options.py
"""Command-line options and parameters for NPA calculations. Replaces onpa/ono_options.java."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any


@dataclass
class OptionParameter:
    """A single configurable parameter with a default value."""
    name: str = ""
    default: Any = None
    description: str = ""
        
    def get_string(self) -> str:
        return str(self.default) if self.default is not None else ""
        
    def get_int(self) -> int:
        return int(self.default) if self.default is not None else 0
        
    def get_boolean(self) -> bool:
        return bool(self.default)
        
    def get_double(self) -> float:
        return float(self.default) if self.default is not None else 0.0


@dataclass
class NPAOptions:
    """Container for all NPA algorithmic and I/O options."""
    glbPrint: bool = False
    dont_print_matrices: bool = True
    
    Input_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-i", "", "input molden file"))
    Charges_File: OptionParameter = field(default_factory=lambda: OptionParameter("-npacharges", "", "file for exporting plain list of NPA charges"))
    NAO_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NAO_Molden_File", "", "Export NAOs in MOLDEN format (AO basis)"))
    LHO_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-LHO_Molden_File", "", "Export LHOs in MOLDEN format (AO basis)"))
    CLPO_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-CLPO_Molden_File", "", "Export CLPOs in MOLDEN format (AO basis)"))
    PrintGeometry: OptionParameter = field(default_factory=lambda: OptionParameter("-PrintGeom", False, "print geometry of the system."))
    WiebergBondOrders_File: OptionParameter = field(default_factory=lambda: OptionParameter("-WibergBondOrders_File", "", "file for exporting a matrix of Wiberg bond orders \n calculated in NAO basis"))
    Edges: OptionParameter = field(default_factory=lambda: OptionParameter("-edges", "", "file for exporting a matrix of Wiberg bond orders \n calculated in NAO basis"))
    D_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-D_Matrix_File", "", "a file name for D (density) matrix D export"))
    S_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-S_Matrix_File", "", "a file name for S (overlap) matrix export"))
    SDS_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-SDS_Matrix_File", "", "a file name for SDS matrix export"))
    PNAO_OverlapMatrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-PNAO_OverlapMatrix_File", "", "a file name for PNAO overlap matrix \n export"))
    PNAO_SDS_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-PNAO_SDS_Matrix_File", "", "a file name for exporting SDS matrix in \n PNAO basis"))
    NMB_old_Overlap_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NMB_old_Overlap_Matrix_File", "", ""))
    NMB_old_SDS_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NMB_old_SDS_Matrix_File", "", ""))
    NRB_old_Overlap_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NRB_old_Overlap_Matrix_File", "", ""))
    NRB_new_Overlap_Matrix_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NRB_new_Overlap_Matrix_File", "", ""))
    S_Matrix_after_ON2_File: OptionParameter = field(default_factory=lambda: OptionParameter("-S_Matrix_after_ON2_File", "", "orbital overlap matrix after 2-nd intracenter \n naturalizatoin transformation"))
    SDS_Matrix_after_ON2_File: OptionParameter = field(default_factory=lambda: OptionParameter("-SDS_Matrix_after_ON2_File", "", ""))
    NRB_Overlap_after_OW_heavy_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NRB_Overlap_after_OW_heavy_File", "", ""))
    NRB_Overlap_after_OS2_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NRB_Overlap_after_Schmidt2_File", "", ""))
    NRB_Overlap_after_OW2_final_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NRB_Overlap_after_OW2_final_File", "", ""))
    OW2_File: OptionParameter = field(default_factory=lambda: OptionParameter("-OW2_File", "", "2-nd WSW transformation matrix (within _full_ \n NRB subspace)"))
    S_Matrix_after_OW2_File: OptionParameter = field(default_factory=lambda: OptionParameter("-S_Matrix_after_OW2_File", "", ""))
    SDS_Matrix_after_OW2_File: OptionParameter = field(default_factory=lambda: OptionParameter("-SDS_Matrix_after_OW2_File", "", ""))
    SDS_NAO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-SDS_NAO_File", "", "a file for exporting SDS matrix in NAO basis"))
    NAO2AO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-NAO2AO_File", "", "a AO->NAO transformation matrix"))
    PNAO2AO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-PNAO2AO_File", "", "a AO->PNAO transformation matrix"))
    CLPO2LHO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-CLPO2LHO_File", "", "a LHO->CLPO transformation matrix"))
    LHO2NAO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-LHO2NAO_File", "", "a NAO->LHO transformation matrix"))
    LPO2AHO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-LPO2AHO_File", "", "a AHO->LPO transformation matrix"))
    AHO2NAO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-AHO2NAO_File", "", "a NAO->AHO transformation matrix"))
    do_Fock: OptionParameter = field(default_factory=lambda: OptionParameter("-doFock", False, "whether to perform Fock matrix processing"))
    Fock_AO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-Fock_AO_File", "", "Fock matrix in AO basis"))
    Fock_NAO_File: OptionParameter = field(default_factory=lambda: OptionParameter("-Fock_NAO_File", "", "Fock matrix in NAO basis"))
    PrintNBMSubmatrices: OptionParameter = field(default_factory=lambda: OptionParameter("-printnmbfock", False, "prints Fock and overlap matrices in NMB PNAO/NAO \n bases (works only when Fock AO matrix is available)"))
    PNAO_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-PNAO_Molden_File", "", "Export PNAOs in MOLDEN format (AO basis)"))
    AHO_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-AHO_Molden_File", "", "Export AHOs in MOLDEN format (AO basis)"))
    LPO_Molden_File: OptionParameter = field(default_factory=lambda: OptionParameter("-LPO_Molden_File", "", "Export LPOs in MOLDEN format (AO basis)"))
    MatrixLineWidth: OptionParameter = field(default_factory=lambda: OptionParameter("-MatrixMaxValuesPerLine", 0, "a maximum number of values per line for matrix \n export (0 = no limit)"))
    MatrixFloatNumberFormat: OptionParameter = field(default_factory=lambda: OptionParameter("-MatrixFloatNumberFormat", "%12.5f", "a printf-compatible format for matrix elements \n printing"))
    maxClpoBondIonicityThreshold: OptionParameter = field(default_factory=lambda: OptionParameter("-maximumBondIonicity", 0.90, "maximum allowable ionicity of CLPO bonding (BD) orbitals"))
    RyOccPrintThreshold: OptionParameter = field(default_factory=lambda: OptionParameter("-RyOccPrintThreshold", 1e-3, "occupancy threshold for RY orbital printing\n set to negative value to print all RY orbitals"))
    HybrOptOccConvThresh: OptionParameter = field(default_factory=lambda: OptionParameter("-HybrOptOccConvThresh", 1e-5, "numerical convergence threshold for hybrid optimization"))
    HybrOptMaxIter: OptionParameter = field(default_factory=lambda: OptionParameter("-HybrOptMaxIter", 1000, "maximum number of iterations during hybrids optimization"))
    VerbosePrint: OptionParameter = field(default_factory=lambda: OptionParameter("-verboseprint", False, "print some additional information"))
    
    bf_nrm2_dev_threshold: float = 1.0E-07
    mo_nrm2_dev_threshold: float = 1.0E-04
    mo_off_diag_threshold: float = 1.0E-04
    S_AO_lin_dep_thresh: float = 1.0E-9
    fail_if_lin_dep: bool = False
    mo_non1_renorm_thresh: float = 1.0E-04
    do_BS_renormalize: bool = True
    do_MO_renormalize: bool = False
    
    NRB_OW_Straightforward: OptionParameter = field(default_factory=lambda: OptionParameter("-directNRBwsw", False, "whether not to use Schmidt-Lowdin scheme for (WSW)^(1/2) \n calculation at spet 6"))
    Heavy_NRB_Threshold: OptionParameter = field(default_factory=lambda: OptionParameter("-heavyNRBthreshold", 1.0E-4, "a threshold value for 'heavily' occupied NRBs (valid only \n without -directNRBwsw)"))
    Overlap_Naive: OptionParameter = field(default_factory=lambda: OptionParameter("-overlap_straightforward", False, "whether to apply slow-but-reliable method to compute \n an overlap matrix (not recommended)"))
    Parameter_File: OptionParameter = field(default_factory=lambda: OptionParameter("-p", "", "a file with parameters (having lower priority than the \n command line)"))
