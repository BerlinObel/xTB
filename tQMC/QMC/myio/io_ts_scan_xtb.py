import numpy as np
import sys
sys.path.append("/groups/kemi/ree/opt/tQMC/QMC")
import xyz2mol.xyz2mol as x2m
#import xyz2mol as x2m

def read_xtb_out(content, quantity='ts_guess_energy'):
    """Reads gaussian output file

    - quantity = 'ts_guess_structure' - structure of the highest energy structure from the xtbscan.log file.
    - quantity = 'ts_guess_atomic_numbers' - atmoic numbers of the highest energy structure from the scan
    - quantity = 'ts_guess_energy' - energy of the maximum of the scan from xtbscan.log.
    """

    if quantity == 'ts_guess_structure':
        return read_ts_guess_structure(content)

    elif quantity == 'ts_guess_atomic_numbers':
        return read_ts_guess_atomic_numbers(content)

    elif quantity == 'ts_guess_energy':
        ts_guess_energy, ts_index = read_energy(content)
        return ts_guess_energy



def read_ts_guess_energy(content):
    """Read total electronic energy """
    
    energies = []

    for line in content.split('\n'):
        if 'SCF done' in line:
            energies.append(float(line.strip().split()[2]))

    energies = np.array(energies)
    ts_guess_energy = np.max(energies)
    ts_index = np.argmax(energies)
    
    return ts_guess_energy, ts_index


def read_ts_guess_structure(content):
    """Read structure from xtbscan.log file """
   
    atom_positions = []
    number_of_atoms = int(content.split('\n')[0].strip())
    ts_guess_energy, ts_index = read_ts_guess_energy(content)
    ts_index = int(ts_index)

    i = 0
    for line in content.split('\n'):
        i += 1
        #if i > (number_of_atoms+2)*ts_index and i < (number_of_atoms+2)*ts_index+number_of_atoms+3: if one only wants the xyz file
        if i > (number_of_atoms+2)*ts_index+2 and i < (number_of_atoms+2)*ts_index+number_of_atoms+3:
            line = line.strip()
            tmp_line = line.split()
            
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            try:
                atom_position = list(map(float, tmp_line[1:]))
                #atom_position = [i*0.529177249 for i in atom_position] #convert Bohr to Angstrom
                atom_positions.append(atom_position)
            except:
                raise ValueError('Expected a line with one string and three floats.')
            
    return atom_positions



def read_ts_guess_structure2(content):
    """Read xyz file from xtbscan.log """

    atom_positions = []
    number_of_atoms = int(content.split('\n')[0].strip())
    ts_guess_energy, ts_index = read_ts_guess_energy(content)
    ts_index = int(ts_index)
    
    i = 0
    for line in content.split('\n'):
        i += 1
        if i > (number_of_atoms+2)*ts_index and i < (number_of_atoms+2)*ts_index+number_of_atoms+3:
            line = line.strip()
            atom_positions.append(line+'\n')
    return ''.join(atom_positions)



def read_ts_guess_structure3(content):
    """Read xyz file from xtbscan.log """

    atom_positions = []
    number_of_atoms = int(content.split('\n')[0].strip())
    ts_guess_energy, ts_index = read_ts_guess_energy(content)
    ts_index = int(ts_index)
    
    text = open(sys.argv[1]).readlines() 
    del text[(number_of_atoms+2)*ts_index+number_of_atoms+2:], text[:(number_of_atoms+2)*ts_index]
    
    f = open(sys.argv[1]+"_new", "w")
    f.writelines(text)
    f.close()

    return open(sys.argv[1]+"_new", 'r').read()



def read_ts_guess_atomic_numbers(content):
    """Read structure from xtbscan.log file """

    atom_symbols = []
    number_of_atoms = int(content.split('\n')[0].strip())
    ts_guess_energy, ts_index = read_ts_guess_energy(content)
    ts_index = int(ts_index)

    i = 0
    for line in content.split('\n'):
        i += 1

        if i > (number_of_atoms+2)*ts_index+2 and i < (number_of_atoms+2)*ts_index+number_of_atoms+3:
            line = line.strip()
            tmp_line = line.split()

            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            try:
                atom_symbol = str(tmp_line[0])
                atom_symbols.append(atom_symbol)
            except:
                raise ValueError('Expected a line with one string and three floats.')

    # atom symbols to atom numbers
    atomic_numbers = list()
    for atom in atom_symbols:
        #atomic_numbers.append(get_atom(atom))
        atomic_numbers.append(x2m.get_atom(atom))

    return atomic_numbers


if __name__ == '__main__':
    import sys

    with open(sys.argv[1], 'r') as out:
        output = out.read()

    print(read_ts_guess_structure2(output))
    #print(read_ts_guess_energy(output))
    #print(read_ts_guess_structure(output))
    #print(read_ts_guess_atomic_numbers(output))
