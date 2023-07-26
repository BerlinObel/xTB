import numpy as np
import xyz2mol.xyz2mol as x2m
#import xyz2mol as x2m
#from xyz2mol import get_atom

def read_xtb_out(content, quantity='energy'):
    """Reads gaussian output file

    - quantity = 'structure' - final structure form output.
    - quantity = 'atomic_numbers' - atmoic numbers
    - quantity = 'energy' - final energy from output.
    """

    if quantity == 'structure':
        return read_structure(content)

    elif quantity == 'atomic_numbers':
        return read_atomic_numbers(content)

    elif quantity == 'energy':
        return read_energy(content)

    elif quantity == 'converged':
        return read_converged(content)


def read_converged(content):
    """Check if calculation terminated correctly"""
#    converged = False
#    for lines in content.split("\n"):

#        if 'normal termination of xtb' in lines:
#            converged = True
#    return converged
    return True


def read_energy(content):
    """Read total electronic energy """

    for line in content.split('\n'):

        if 'TOTAL ENERGY' in line:
            energy = float(line.strip().split()[3])

    return energy


def read_structure(content):
    """Read structure from output file """
    
    temp_items = content.split('$coord')[1:]
    
    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]

        atom_positions = []

        for line in lines:
            line = line.strip()
            # if line equals '$end' - mol block ends.
            if set(line) == set('$end'):
                break
            
            tmp_line = line.split()
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions
            try:
                atom_position = list(map(float, tmp_line[:3]))
                atom_position = [i*0.529177249 for i in atom_position] #convert Bohr to Angstrom
                atom_positions.append(atom_position)
            except:
                raise ValueError('Expected a line with one string and three floats.')
            
    return atom_positions


def read_atomic_numbers(content):
    """Read structure from output file """

    temp_items = content.split('$coord')[1:]
    
    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]

        atom_symbols = []

        for line in lines:
            line = line.strip()
            
            # is line equal '$end' - mol block ends.
            if set(line) == set('$end'):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions
            try:
                atom_symbol = str(tmp_line[-1])
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


    print(read_energy(output))
    print(read_structure(output))
    print(read_atomic_numbers(output))
    #for x in read_atomic_numbers(output):
    #    #print(x)
    #    pass
