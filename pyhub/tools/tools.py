import os
import numpy as np

def find_file(directory, filename):
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    raise FileNotFoundError(f"Le fichier {filename} n'a pas été trouvé dans le répertoire {directory} ni dans ses sous-répertoires.")


def count_bit(n):
    '''
    Count the number of ones in bitstring given by integer n
    '''
    count = 0
    while n:
        count += 1
        n &= (n - 1)
    return count


def int2str(integer,lengh):
    integer_str = str(integer)
    while len(integer_str)<lengh:
        integer_str = '0'+integer_str
    return integer_str

def fidelity(psi, phi):
    return np.absolute(np.dot(np.conjugate(psi), phi))**2