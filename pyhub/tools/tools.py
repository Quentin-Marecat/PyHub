import os

def find_file(directory, filename):
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    raise FileNotFoundError(f"Le fichier {filename} n'a pas été trouvé dans le répertoire {directory} ni dans ses sous-répertoires.")