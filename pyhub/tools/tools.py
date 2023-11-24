import os

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