import os
import itertools
from ase.io import read, write
from ase.build.tools import stack
from matminer.featurizers.conversions import ASEAtomstoStructure

def calculate_lattice_mismatch(lattice_parameters1, lattice_parameters2):
    a1, b1, c1 = lattice_parameters1[:3]
    a2, b2, c2 = lattice_parameters2[:3]
    # mismatch = abs(a1 - a2) / ((a1 + a2) / 2) * 100

    if a1 > a2:
        a1, a2 = a2, a1  # Swap a1 and a2

    mismatch = abs((a1 - a2) / a1) * 100
    return mismatch

def create_heterostructure(atoms1, atoms2):
    # Stacking along the z-axis
    heterostructure = stack(atoms1, atoms2, axis=2, maxstrain=None, distance=4)
    return heterostructure

processed_pairs = set()
aa2s = ASEAtomstoStructure()

# 'MoS2/AB2', 'CdI2/AB2', 'GaS/A2B2', 'GeSe/AB', 'BN/AB', 'GaSe/A2B2'

# Bilayer prototypes
prototypes = [
    ("MoS2", "MoS2"), ("CdI2", "CdI2"), ("GaS", "GaS"), ("GaSe", "GaSe"),
    ("GeSe", "GeSe"), ("CdI2", "BN"), ("GaS", "CdI2"), ("GaSe", "GaS"),
    ("GeSe", "BN"), ("GeSe", "CdI2"), ("GeSe", "GaS"), ("GeSe", "GaSe"),
    ("MoS2", "BN"), ("MoS2", "CdI2"), ("MoS2", "GaS"), ("MoS2", "GaSe"),
    ("MoS2", "GeSe"), ("GaSe", "BN"), ("GaSe", "CdI2"), ("BN", "BN"),
    ("GaS", "BN")
]

poscar_count = 0

for prototype1, prototype2 in prototypes:
    # Use prototype names to generate directories
    output_directory = f"Bilayers/{prototype1}-{prototype2}"
    os.makedirs(output_directory, exist_ok=True)

    # Paths to monolayers
    dir1 = os.path.abspath(prototype1)
    dir2 = os.path.abspath(prototype2)

    monolayer_structures1 = []
    monolayer_structures2 = []

    # Gather monolayer structures from prototype1 directory (all subdirectories)
    for root, dirs, files in os.walk(dir1):
        for file in files:
            if file == "structure.json":  # Look for the specific file
                json_path1 = os.path.join(root, file)
                atoms1 = read(json_path1)
                monolayer_structures1.append(atoms1)

    # Gather monolayer structures from prototype2 directory (all subdirectories)
    for root, dirs, files in os.walk(dir2):
        for file in files:
            if file == "structure.json":  # Look for the specific file
                json_path2 = os.path.join(root, file)
                atoms2 = read(json_path2)
                monolayer_structures2.append(atoms2)

    # Create heterostructures for all combinations between prototype1 and prototype2
    for atoms1, atoms2 in itertools.product(monolayer_structures1, monolayer_structures2):
        if atoms1 != atoms2:
            formula1, formula2 = sorted([atoms1.get_chemical_formula(), atoms2.get_chemical_formula()])
            if (formula1, formula2) not in processed_pairs and (formula2, formula1) not in processed_pairs:
                processed_pairs.add((formula1, formula2))
                processed_pairs.add((formula2, formula1))

                # Calculate lattice mismatch and create heterostructure
                mismatch = calculate_lattice_mismatch(atoms1.cell.cellpar(), atoms2.cell.cellpar())
                if mismatch <= 4:
                    heterostructures = []
                    heterostructure = create_heterostructure(atoms1, atoms2)
                    heterostructure.center(vacuum=15.0, axis=2)
                    heterostructures.append(heterostructure)

                    # Save the heterostructure as POSCAR in the corresponding directory
                    hetero_directory = os.path.join(output_directory, f'{formula1}-{formula2}')
                    os.makedirs(hetero_directory, exist_ok=True)
                    write_path = os.path.join(hetero_directory, 'POSCAR')
                    write(write_path, heterostructure, format='vasp')
                    poscar_count += 1
                    # print(f"Saved POSCAR for {formula1}-{formula2} in {write_path}")

print(f"Total number of POSCAR files created: {poscar_count}")

# python3 hetero_formation.py > bilayers.log 2>&1