from tequila.quantumchemistry.pyscf_interface import QuantumChemistryPySCF
import os
import numpy
from pyscf import scf,mp
from pyscf.tools import molden
from copy import deepcopy
import subprocess
from tequila.quantumchemistry.qc_base import QuantumChemistryBase
from sunrise.molecules.hybrid_base import HybridBase
from sunrise.molecules.fermionic_base import FermionicBase
from sunrise.CLPO import call_molden2aim
from copy import deepcopy
from typing import Tuple
from numbers import Number
from sunrise import from_tequila
from tequila import TequilaException
import shlex
from numpy import array
from warnings import warn
from pathlib import Path

def __transform(modified:QuantumChemistryBase,original:QuantumChemistryBase=None,orbital_type='CLPO')->Tuple[QuantumChemistryBase,dict]:
    '''
    Procedure similar to what is done in use_native_orbitals but for arbitrary basis. Keeps frozen orbitals canHF
    orthogonalized with the active modified ones
    Returns modified molecule with the core orbitals of the original one
    And a dictionary with the form {active_orbital_index_before:active_orbital_index_after}
    The frozen orbitals will always be the N first on the orbital matrix
    '''
    def inner(a, b, s):
        return numpy.sum(numpy.multiply(numpy.outer(a, b), s))
    core = [i.idx_total for i in original.integral_manager.orbitals if i.idx is None]
    assert len(original.integral_manager.orbitals) == len(modified.integral_manager.orbitals)
    d = deepcopy(modified.integral_manager.orbital_coefficients).T
    c = deepcopy(original.integral_manager.orbital_coefficients).T
    s = original.integral_manager.overlap_integrals
    n_basis = len(d)
    ov = numpy.zeros(shape=(n_basis))
    for i in core:
        for j in range(n_basis):
            ov[j] += numpy.abs(inner(c[i], d[j], s))
    co = {}
    for i in core:
        idx = numpy.argmax(ov)
        co[i] = idx
        ov[idx] = 0
    active = [i for i in range(n_basis) if i not in co.values()]
    to_active =  [i for i in range(n_basis) if  i not in co.keys()]
    to_active = {active[i]:to_active[i] for i in range(len(active))}
    reference_orbitals = [*co.keys()]
    i =0
    while len(reference_orbitals)<original.parameters.total_n_electrons//2:
        if i not in reference_orbitals:
            reference_orbitals.append(i)
        i += 1
    sbar = numpy.zeros(shape=s.shape)
    for k in active:
        for i in core:
            sbar[i][to_active[k]] = inner(c[i], d[k], s)
    dbar = numpy.zeros(shape=s.shape)

    for j in active:
        dbar[to_active[j]] = d[j]
        for i in core:
            temp = sbar[i][to_active[j]] * c[i]
            dbar[to_active[j]] -= temp
    for i in to_active.values():
        norm = numpy.sqrt(inner(dbar[i], dbar[i], s.T))
        if not numpy.isclose(norm, 0):
            dbar[i] = dbar[i] / norm
    for j in to_active.values():
        c[j] = dbar[j]
    sprima = numpy.eye(len(c))
    for idx, i in enumerate(to_active.values()):
        for j in [*to_active.values()][idx:]:
            sprima[i][j] = inner(c[i], c[j], s)
            sprima[j][i] = sprima[i][j]
    lam_s, l_s = numpy.linalg.eigh(sprima)
    for ei,e in enumerate(lam_s):
        if numpy.isclose(e,0,atol=1.e-9):
            lam_s[ei] = max(e,1.e-9) # Note: Fix to avoid inestabilities 
    lam_s = lam_s * numpy.eye(len(lam_s))
    lam_sqrt_inv = numpy.sqrt(numpy.linalg.inv(lam_s))
    symm_orthog = numpy.dot(l_s, numpy.dot(lam_sqrt_inv, l_s.T))
    jcoef = symm_orthog.dot(c).T
    ref = [i.idx_total for i in original.integral_manager.reference_orbitals if i not in original.integral_manager.active_reference_orbitals]
    ref.extend([i for i in range(n_basis) if  i not in co.keys()][:len(original.integral_manager.active_reference_orbitals)])
    integral_manager = modified.initialize_integral_manager(one_body_integrals=original.integral_manager.one_body_integrals,
                    two_body_integrals=original.integral_manager.two_body_integrals,constant_term=original.integral_manager.constant_term,
                    active_orbitals= [i for i in range(n_basis) if  i not in co.keys()],frozen_orbitals=[*co.keys()],orbital_coefficients=jcoef,
                    overlap_integrals=original.integral_manager.overlap_integrals,reference_orbitals=ref,orbital_type=orbital_type)
    parameters = deepcopy(original.parameters)
    if isinstance(modified,FermionicBase):
        return FermionicBase(parameters=parameters,integral_manager=integral_manager,fermionic_backend=modified.fermionic_backend),to_active
    elif isinstance(modified,HybridBase):
        return HybridBase(parameters=parameters,integral_manager=integral_manager,transformation=modified.transformation,select=modified.select,two_qubit=modified.two_qubit, condense=modified.condense) ,to_active
    return QuantumChemistryBase(parameters=parameters,integral_manager=integral_manager,transformation=modified.transformation),to_active

def generate_molden(mol:QuantumChemistryBase,filename:str=None,output_dir:str=None,mo_occ:list=None,mo_energy:list=None,use_mp2:bool=False,option1:bool=True,use_active:bool=False):
    '''
    Interface with pyscf.tools molden file generation

    :param mol: Any kind of tequila/sunrise Molecule
    :param filename: The moldenfile will be saved as filename.molden. If None, the molecule.parameters.name is employed
    :param output_dir: default None = working file.
    :param mo_occ: Molecular Orbital electronic occupation. If None, hf or mp2 are used depending on use_mp2
    :param mo_energy: Molecular Orbital energy. If None, hf or mp2 are used depending on use_mp2
    :param use_mp2: Whether to us mp2 or hf if no mo_occ and mo_energy are provided
    :param option1: Whether to use the first or second molden generation alternatives propsed by pyscf. See 
                    https://github.com/pyscf/pyscf/blob/master/examples/tools/02-molden.py for more info
    :param use_active: Whether to include only the active orbitals.
    '''

    assert (mo_occ == None) and (mo_energy == None)
    size_basis = len(mol.integral_manager._orbital_coefficients)
    active = mol.integral_manager.active_space.active_orbitals
    pfmol = from_tequila(mol)
    if output_dir is None:
        output_dir = os.getcwd()
    if mo_occ is None:
        if use_mp2:
            mo_occ, mo_energy = __get_MP2_occ(mol)
        else:
            mf = scf.RHF(pfmol).run()
            mo_occ = mf.mo_occ
            mo_energy = mf.mo_energy
    else:
        assert len(mo_occ) == len(mo_energy)
        if (use_active and len(mo_occ) == len(active) or len(mo_occ) == len(size_basis)) or (not use_active and len(mo_occ) == len(size_basis)):
            pass
        else:
            raise TequilaException(f"{len(mo_occ)} Molecular Orbital Occupation but it doesn't fit not the total space size ({size_basis}) nor the active space size ({len(active)}). Use_active={use_active}")

    if filename is None: filename=mol.parameters.name


    mo_coeff = mol.integral_manager.orbital_coefficients
    if use_active:
        if len(mo_occ) == size_basis:
            mo_occ = [mo_occ[i] for i in active]
            mo_energy = [mo_energy[i] for i in active]
        mo_coeff = mo_coeff[:,active]

    if option1:
        ### OPTION 1
        with open(f'{output_dir}/{filename}.molden', 'w') as f1:
            molden.header(pfmol, f1)
            molden.orbital_coeff(pfmol, f1, mo_coeff, ene=mo_energy, occ=mo_occ)
    else:
        ### OPTION 2
        try:
            molden.from_mo(pfmol, f'{output_dir}/{filename}.molden',mo_coeff,ene=mo_energy,occ=mo_occ)
        except RuntimeError:
            print('    Found l=5 in basis.')
            molden.from_mo(pfmol, f'{output_dir}/{filename}.molden', mo_coeff,ene=mo_energy,occ=mo_occ,ignore_h=True)

def generate_CLPO_molecule_edges(mol:QuantumChemistryBase, output_dir:str=None, thres:Number=1.e-12, silent:bool=True, use_active:bool=True, rm_files:bool=True, **kwargs)->Tuple[QuantumChemistryBase,list]:
    '''
    Temporal function for generating a molecule with CLPO orbitals (10.1002/qua.25798) until integrated in Sunrise molecules
    
    :param mol: Any kind of Tequila/Sunrise Molecules.
    :param output_dir: default None = working file.
    :param thres: -HybrOptOccConvThresh from Janpa. Default 1.e-9.
    :param use_active: Whether to respect input molecule frozen/active space, letting frozen as HF
    :param kwargs: keywords accepted by 'generate_molden', see above.

    Return modified Molecule and SPA edges
    '''
    if 'filename' in kwargs:
        filename= kwargs['filename']
        kwargs.pop('filename')
    else: filename = mol.parameters.name
    if output_dir is None:
        output_dir = os.getcwd()
    generate_molden(mol=mol,filename=filename,output_dir=output_dir,use_active=False,**kwargs) #TODO: Janpa CLPO is bug for active space only, working on 
    call_molden2aim(moldenfile=filename+'.molden',output_dir=output_dir)
    call_molden2molden(command=f'-NormalizeBF -cart2pure  -i {filename}.molden -o {filename}.molden',silent=True)
    call_janpa(command=f'-i {filename}.molden -CLPO_Molden_File {filename}_CLPO.molden -HybrOptOccConvThresh {thres}',silent=silent)
    mo_matrix = read_molden_mo_matrix(f"{filename}_CLPO.molden")
    if rm_files:
        subprocess.call(f'rm {output_dir}/m2a.ini',shell=True)
        subprocess.call(f'rm {output_dir}/{filename}.molden',shell=True)
        subprocess.call(f'rm {output_dir}/{filename}_new.molden',shell=True)
        subprocess.call(f'rm {output_dir}/{filename}_CLPO.molden',shell=True)
    nmol = deepcopy(mol)
    nmol.integral_manager.orbital_coefficients = mo_matrix
    if use_active:
        mol,to_active = __transform(original=mol,modified=nmol)
    else: mol = nmol
    graph = extract_clpo_graph(f"{output_dir}/graph")
    if use_active:
        ncore = len(mol.integral_manager.orbital_coefficients)-mol.n_orbitals
        graph = [tuple([to_active[i] - ncore for i in edge if i in to_active.keys()]) for edge in graph]
        graph = [g for g in graph if len(g)]
    if rm_files:
        subprocess.call(f'rm {output_dir}/graph',shell=True)
    return mol,graph

def generate_HAO_molecule(mol:QuantumChemistryBase, output_dir:str=None, thres:Number=1.e-9, silent:bool=True, use_active:bool=True, rm_files:bool=True,**kwargs)->QuantumChemistryBase:
    '''
    Temporal function for generating a molecule with Hybrid Atomic Orbitals via janpa (10.1002/qua.25798) until integrated in Sunrise molecules
    
    :param mol: Any kind of Tequila/Sunrise Molecules.
    :param output_dir: default None = working file.
    :param thres: -HybrOptOccConvThresh from Janpa. Default 1.e-9.
    :param use_active: Whether to respect input molecule frozen/active space, letting frozen as HF
    :param kwargs: keywords accepted by 'generate_molden', see above.

    Return modified Molecule and SPA edges
    '''
    if 'filename' in kwargs:
        filename= kwargs['filename']
        kwargs.pop('filename')
    else: filename = mol.parameters.name
    if output_dir is None:
        output_dir = os.getcwd()

    generate_molden(mol=mol,filename=filename,output_dir=output_dir,**kwargs)
    call_molden2aim(moldenfile=filename+'.molden',output_dir=output_dir)
    call_molden2molden(command=f'-NormalizeBF -i {filename}.molden -o {filename}.molden',silent=True)
    call_janpa(command=f'-i {filename}.molden -AHO_Molden_File {filename}_HAO.molden -HybrOptOccConvThresh {thres}',silent=silent,output_dir=output_dir)
    mo_matrix = read_molden_mo_matrix(f"{output_dir}/{filename}_HAO.molden")
    if rm_files:
        subprocess.call(f'rm {output_dir}/m2a.ini',shell=True)
        subprocess.call(f'rm {output_dir}/{filename}.molden',shell=True) 
        subprocess.call(f'rm {output_dir}/{filename}_new.molden',shell=True)
        subprocess.call(f'rm {output_dir}/{filename}_HAO.molden',shell=True)
    nmol = deepcopy(mol)
    nmol.integral_manager.orbital_coefficients = mo_matrix
    if use_active:
        mol,to_active = __transform(original=mol,modified=nmol,orbital_type='HAO')
    else: mol = nmol
    return mol

def call_janpa(command: str = '', folder_x: str = None, output_dir=None, silent=False):
    """
    Wrapper for the janpa executable files.
    command: str Command string as it would be introduced on the terminal.
    folder_x: str Path to your specific folder to be processed as input.
    silent: bool Whether to silence janpa output.
    output_dir: default None = working directory.
    """
    
    # 1. Locate the jar file relative to this Python script
    # __file__ gets the path of the current script, .parent gets its directory
    jar_path = Path(__file__).parent / 'janpa.jar'
    
    if not jar_path.exists():
        raise FileNotFoundError(f"Could not find janpa.jar at {jar_path}. Make sure it is in the same folder as this script.")

    if output_dir is None:
        output_dir = os.getcwd()
            
    # 2. Construct the Java command. 
    cmd = ['java', '-jar', str(jar_path)]
    
    # 3. Append the specific command arguments
    if command:
        extra_args = shlex.split(command)
        cmd.extend(extra_args)
        
    # 4. Append your input folder X to the command arguments
    if folder_x:
        cmd.append(os.path.abspath(folder_x))

    # 5. Run the subprocess
    try:
        res = subprocess.run(
            cmd, 
            capture_output=silent, 
            text=True, 
            check=True,
            cwd=output_dir  # Ensures the command runs inside the specified output directory
        )
        return res
    except subprocess.CalledProcessError as e:
        print(f"Error executing JANPA jar: {e.stderr}")
        raise

def call_molden2molden(command: str = '', folder_x: str = None, output_dir=None, silent=False):
    """
    Wrapper for the janpa executable files.
    command: str Command string as it would be introduced on the terminal.
    folder_x: str Path to your specific folder to be processed as input.
    silent: bool Whether to silence janpa output.
    output_dir: default None = working directory.
    """
    
    # 1. Locate the jar file relative to this Python script
    # __file__ gets the path of the current script, .parent gets its directory
    jar_path = Path(__file__).parent / 'molden2molden.jar'
    
    if not jar_path.exists():
        raise FileNotFoundError(f"Could not find molden2molden.jar at {jar_path}. Make sure it is in the same folder as this script.")

    if output_dir is None:
        output_dir = os.getcwd()
            
    # 2. Construct the Java command. 
    cmd = ['java', '-jar', str(jar_path)]
    
    # 3. Append the specific command arguments
    if command:
        extra_args = shlex.split(command)
        cmd.extend(extra_args)
        
    # 4. Append your input folder X to the command arguments
    if folder_x:
        cmd.append(os.path.abspath(folder_x))

    # 5. Run the subprocess
    try:
        res = subprocess.run(
            cmd, 
            capture_output=silent, 
            text=True, 
            check=True,
            cwd=output_dir  # Ensures the command runs inside the specified output directory
        )
        return res
    except subprocess.CalledProcessError as e:
        print(f"Error executing JANPA jar: {e.stderr}")
        raise

def extract_clpo_graph(graph_file:str):
    """
    Reads a CLPO 'graph' file and extracts a graph representation.

    Returns:
        List of tuples:
        - (i,)      for lone pairs
        - (i, i+1)  for BD/NB pairs
    """
    nodes = []

    with open(graph_file, "r") as f:
        lines = f.readlines()

    # Skip header
    data_lines = [line.rstrip() for line in lines if line.strip()][1:]

    i = 0
    while i < len(data_lines):
        line = data_lines[i]
        # Lone pair
        if "(LP)" in line:
            occ = float(line.split()[3])
            if occ < 1.5 and occ > 0.5: # NOTE: arbitrary numbers
                warn(f'Lone Pair {i} found with occupation close to 1 => {occ}, take care.') 
            if occ < 1.: 
                i += 1
                continue
            nodes.append((i,))
            i += 1
            continue

        # Bonding orbital → must pair with next line
        if "(BD)" in line:
            if i + 1 >= len(data_lines):
                raise ValueError("BD entry without following NB line")

            bd_pair_occ = float(line.split()[6])
            at_pair_occ = float(data_lines[i+1].split()[4])
            if bd_pair_occ + at_pair_occ < 1.7: # NOTE: arbitrary number
                warn(f'Bond pair orbitals {[i,i+1]} population expected under expected 2e-, predicted: {bd_pair_occ + at_pair_occ}, take care with predicted edges.')

            nodes.append((i, i + 1))
            i += 2
            continue

        i += 1

    return nodes

def read_molden_mo_matrix(filename:str):
    """
    Reads a Molden file and returns a NumPy array where
    each column is an MO coefficient vector.
    """
    mo_vectors = []
    current_mo = []
    in_mo_section = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Detect MO section
            if line == "[MO]":
                in_mo_section = True
                continue

            if not in_mo_section:
                continue

            if line.startswith("Sym="):
                if current_mo:
                    mo_vectors.append(current_mo)
                    current_mo = []
                continue

            # Skip metadata lines
            if (
                line.startswith("Ene=") or
                line.startswith("Spin=") or
                line.startswith("Occup=") or
                not line
            ):
                continue

            parts = line.split()
            if len(parts) == 2:
                try:
                    coeff = float(parts[1])
                    current_mo.append(coeff)
                except ValueError:
                    pass

        if current_mo:
            mo_vectors.append(current_mo)

    matrix = array(mo_vectors).T
    return matrix


if __name__ == '__main__':
    import  sunrise as sun
    # geometry = '''
    # C	0.0000000	1.3894860	0.0000000
    # C	1.2033300	0.6947430	0.0000000
    # C	1.2033300	-0.6947430	0.0000000
    # C	0.0000000	-1.3894860	0.0000000
    # C	-1.2033300	-0.6947430	0.0000000
    # C	-1.2033300	0.6947430	0.0000000
    # H	0.0000000	2.4702840	0.0000000
    # H	2.1393290	1.2351420	0.0000000
    # H	2.1393290	-1.2351420	0.0000000
    # H	0.0000000	-2.4702840	0.0000000
    # H	-2.1393290	-1.2351420	0.0000000
    # H	-2.1393290	1.2351420	0.0000000
    # Cu  0.00000000  0.0000000   1.8000000 '''
    # geometry = 'Cu 0. 0. 0. '
    geometry = 'H 0. 0. 0. \n H 0. 0. 1. \n H 0. 0. 2. \n H 0. 0. 3.'
    mol = sun.Molecule(geometry=geometry, basis_set='sto-3g', backend='pyscf')
    # sun.plot_MO(mol,'HF')
    print('Total Ne',mol.parameters.total_n_electrons)
    print('Active Ne ',mol.n_electrons)
    print('Core Ne',mol.parameters.get_number_of_core_electrons())
    print('Total Norb ',len(mol.integral_manager.orbitals))
    print('Active Norb',mol.n_orbitals)
    print("Core Norb",len([i.idx_total for i in mol.integral_manager.orbitals  if i.idx is None]))
    mol1,edges = generate_CLPO_molecule_edges(mol,rm_files=False,silent=True)
    print('edges ',edges)
    print("Ne Edges",2*len(edges))
    print("CLPO DONE")
    # mol2 =  generate_HAO_molecule(mol,rm_files=True,silent=True)
    # print("HAO Done")