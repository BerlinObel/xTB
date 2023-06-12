#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python
import sys
import pandas as pd
import copy
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")

from qmmol import QMMol
from qmconf import QMConf
from calculator.xtb import xTB
from calculator.orca import ORCA
from calculator.gaussian import Gaussian

from conformers.create_conformers import RotatableBonds


def map_atoms(reactant, product):
    """Returns new order in product """

    reac = copy.deepcopy(reactant)
    prod = copy.deepcopy(product)

    Chem.Kekulize(reac,clearAromaticFlags=True)
    Chem.Kekulize(prod,clearAromaticFlags=True)
     
    reac = change_mol(reac)
    prod = change_mol(prod)
    
    # Break Bond in reactant, in order to compare to product.
    smarts_bond = Chem.MolFromSmarts('[CX4;H0;R]-[CX4;H1;R]')
    atom_idx = list(reac.GetSubstructMatch(smarts_bond))
    
    if len(atom_idx) != 0:
        bond = reac.GetBondBetweenAtoms(atom_idx[0], atom_idx[1])
        broken_bond_reac = Chem.FragmentOnBonds(reac, [bond.GetIdx()], addDummies=False)
        
        # find new atom order for product
        prod_order = prod.GetSubstructMatch(broken_bond_reac)
    
    else:
        prod_order = prod.GetSubstructMatch(reac)
    
    return prod_order


def change_mol(mol):
    """ Change molecule to only compare conectivity """

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            bond.SetBondType(Chem.BondType.SINGLE)

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(0)

    return mol


def reorder_product(reac, prod):
    """ change atom order of product, to match reactant  """
    new_product_order = map_atoms(reac, prod)
    reordered_product = Chem.RenumberAtoms(prod, new_product_order)

    return reordered_product


def reactant2product(reac_smi):
    """ create prodruct from reactant """
    
    print(reac_smi)
    input_string = reac_smi
    prod_smi = input_string.replace("/N=N/", "/N=N\\")
    print(prod_smi)

    return reac_smi, prod_smi




def gs_conformer_search(name, rdkit_conf, chrg, mult, cpus):
    """ ground state conformer search """
    
    charged = True # hard coded for mogens
    
    # create conformers
    qmmol = QMMol()
    qmmol.add_conformer(rdkit_conf, fmt='rdkit', label=name,
                       charged_fragments=charged, set_initial=True)

    # find #-- rot bond.
    mol = qmmol.get_rdkit_mol()
    triple_smart = '[c]:[c]-[CH0]#[CH0]'
    extra_rot_bond = int(len(mol.GetSubstructMatches(Chem.MolFromSmarts(triple_smart))) / 2)

    # print compute number of conformers to find
    rot_bonds = len(RotatableBonds(qmmol.initial_conformer.get_rdkit_mol())) + extra_rot_bond
    num_confs = 5 + 5*rot_bonds
    qmmol.create_random_conformers(threads=cpus, num_confs=num_confs)

    xtb_params = {'method': 'gfn1',
                  'opt': 'tight',
                  'cpus': 1}

    qmmol.calc = xTB(parameters=xtb_params)
    qmmol.optimize(num_procs=cpus, keep_files=False)

    # Get most stable conformer. If most stable conformer
    # not identical to initial conf try second lowest.
    initial_smi = Chem.MolToSmiles(Chem.RemoveHs(qmmol.initial_conformer.get_rdkit_mol()))
    low_energy_conf = qmmol.nlowest(1)[0]
    conf_smi = Chem.MolToSmiles(Chem.RemoveHs(low_energy_conf.get_rdkit_mol()))

    i = 1
    while initial_smi != conf_smi:
        low_energy_conf = qmmol.nlowest(i+1)[-1]
        conf_smi = Chem.MolToSmiles(Chem.RemoveHs(low_energy_conf.get_rdkit_mol()))
        i += 1
        
        if len(qmmol.conformers) < i:
            sys.exit('no conformers match the initial input')

    return low_energy_conf


def gs_mogens(name, smi, chrg, mult, cps):
    """GS conformers search given a smiles string  """
    print('gs_mogens reached')    
    reac_smi, prod_smi = reactant2product(smi)
    print(reac_smi,prod_smi)
    reac_mol = Chem.AddHs(Chem.MolFromSmiles(reac_smi))
    try:
        print("Trying reorder")
        prod_mol = reorder_product(reac_mol, Chem.AddHs(Chem.MolFromSmiles(prod_smi))) #THIS IS THE ACTUAL PLACE THAT FAILS!
    except:
        print("Reorder failed,using no reorder")
        prod_mol = Chem.AddHs(Chem.MolFromSmiles(prod_smi))
    for i, comp in enumerate([(reac_mol, "_r"), (prod_mol, "_p")]):
        mol, p_r = comp
        
        AllChem.EmbedMolecule(mol)
        rdkit_conf = mol.GetConformer()

        # create mol for QMMol using rdkit
        n = name + p_r
        
        if p_r == '_r':
            try:
                print("Starting reac conformer search")
                reac_qmconf = gs_conformer_search(n, rdkit_conf, chrg, mult, cps)
            except:
                print("Original qmconf for reac failed, thus algo turns to no reorder")
                count = 0
                reac_qmconf = None
                while reac_qmconf is None and count < 40:
                    try:
                        reac_qmconf = gs_retry(chrg, mult, cps, name, reac_smi,p_r)
                    except:
                        count += 1
                        print("Failed retry:",count)
                        pass
        if p_r == '_p':
            try:
                print("Starting prod conformer search")
                prod_qmconf = gs_conformer_search(n, rdkit_conf, chrg, mult, cps)
            except:
                print("Original qmconf for prod failed, thus algo turns to no reorder")
                count = 0
                prod_qmconf = None
                while prod_qmconf is None and count < 40:
                    try:
                        prod_qmconf = gs_retry(chrg, mult, cps, name, prod_smi,p_r)
                    except:
                        count += 1
                        print("Failed retry:",count)
                        pass
    storage = prod_qmconf.results['energy'] - reac_qmconf.results['energy']

    return reac_qmconf, prod_qmconf, storage



def gs_retry(chrg, mult, cps, name, prod_smi,p_r):
    
    mol = Chem.AddHs(Chem.MolFromSmiles(prod_smi))
    AllChem.EmbedMolecule(mol)
    rdkit_conf = mol.GetConformer()
    n = name + p_r
    return gs_conformer_search(n, rdkit_conf, chrg, mult, cps)


if __name__ == '__main__':
    
    import pandas as pd
    import sys

    cpus = 2

    data = pd.read_csv(sys.argv[1])
    print(data)
    # find storage energy
    compound_list = list()
    for idx, compound in data.iterrows():
        try: 
            reac_qmconf, prod_qmconf, storage = gs_mogens(str(compound.comp_name),
                                                          compound.smiles,
                                                          compound.charge,
                                                          compound.multiplicity,
                                                          cpus)

            print(str(compound.comp_name),compound.smiles,compound.charge,compound.multiplicity,cpus) 
            compound_list.append({'rep': str(compound.comp_name),
                                  'reac': reac_qmconf,
                                  'prod': prod_qmconf,
                                  'storage': storage})
        except:
            pass

    results_df = pd.DataFrame(compound_list)
    pickle_name = sys.argv[1].split('/')[0].split('.')[0] + '.pkl'
    print(f"pickle name: {pickle_name}")
    print(results_df)
    results_df.to_pickle(pickle_name)
    
