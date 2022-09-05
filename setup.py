#!/usr/bin/env python3

from cx_Freeze import setup, Executable

base = None
# Remplacer "monprogramme.py" par le nom du script qui lance votre programme
executables = [Executable("evofond.py", base=base)]
# Renseignez ici la liste complète des packages utilisés par votre application
packages = ["numpy", "json", "os", "shutil", "scipy", "matplotlib", "datetime", "logging", "argparse"]
options = {
    'build_exe': {    
        'packages':packages,
    },
}
# Adaptez les valeurs des variables "name", "version", "description" à votre programme.
setup(
    name = "evofond",
    options = options,
    version = "1.0",
    description = "Evofond, programme de calcul d'evolution de fond de lit lors de crues torrentielles",
    executables = executables
)