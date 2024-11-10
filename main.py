# Define the particles
class Proton:
    def __init__(self):
        self.charge = +1  # positive charge
        self.mass = 1.6726219e-27  # mass in kg


class Neutron:
    def __init__(self):
        self.charge = 0  # neutral
        self.mass = 1.674929e-27  # mass in kg


class Electron:
    def __init__(self):
        self.charge = -1  # negative charge
        self.mass = 9.10938356e-31  # mass in kg


# Define the Atom class
class Atom:
    def __init__(self, atomic_number, neutrons):
        self.atomic_number = atomic_number
        self.neutrons = neutrons
        self.protons = [Proton() for _ in range(atomic_number)]
        self.electrons = [Electron() for _ in range(atomic_number)]
        self.neutrons_list = [Neutron() for _ in range(neutrons)]
        
        # Calculate atomic mass
        self.atomic_mass = atomic_number * self.protons[0].mass + neutrons * self.neutrons_list[0].mass
        self.bonds = []  # [(other_atom, bond_type), ...]

    def bond_with(self, other_atom, bond_type="single"):
        self.bonds.append((other_atom, bond_type))
        other_atom.bonds.append((self, bond_type))

    def electron_configuration(self):
        configuration = []
        remaining_electrons = self.atomic_number
        shell_capacity = [2, 8, 18, 32, 32, 18, 8, 2]
        
        for capacity in shell_capacity:
            if remaining_electrons > 0:
                electrons_in_shell = min(remaining_electrons, capacity)
                configuration.append(electrons_in_shell)
                remaining_electrons -= electrons_in_shell
            else:
                break
        return configuration

    def summary(self):
        print(f"Atomic Number: {self.atomic_number}")
        print(f"Number of Neutrons: {len(self.neutrons_list)}")
        print(f"Atomic Mass: {self.atomic_mass} kg")
        print(f"Electron Configuration: {self.electron_configuration()}")


# Define the Molecule class
class Molecule:
    def __init__(self, name="Unknown Molecule"):
        self.name = name
        self.atoms = []
    
    def add_atom(self, atom):
        self.atoms.append(atom)

    def calculate_molecular_mass(self):
        return sum(atom.atomic_mass for atom in self.atoms)

    def formula(self):
        atom_counts = {}
        for atom in self.atoms:
            atomic_number = atom.atomic_number
            atom_counts[atomic_number] = atom_counts.get(atomic_number, 0) + 1
        return " ".join(f"{atomic_number}{count if count > 1 else ''}" for atomic_number, count in atom_counts.items())

    def summary(self):
        print(f"Molecule: {self.name}")
        print(f"Chemical Formula: {self.formula()}")
        print(f"Molecular Mass: {self.calculate_molecular_mass()} kg")
        print("Bonding Structure:")
        for atom in self.atoms:
            for bond in atom.bonds:
                bonded_atom, bond_type = bond
                print(f"  Atom {atom.atomic_number} - {bond_type} bond - Atom {bonded_atom.atomic_number}")


# Example: Creating H2SO4 (Sulfuric Acid)
hydrogen1 = Atom(atomic_number=1, neutrons=0)
hydrogen2 = Atom(atomic_number=1, neutrons=0)
sulfur = Atom(atomic_number=16, neutrons=16)
oxygen1 = Atom(atomic_number=8, neutrons=8)
oxygen2 = Atom(atomic_number=8, neutrons=8)
oxygen3 = Atom(atomic_number=8, neutrons=8)
oxygen4 = Atom(atomic_number=8, neutrons=8)

sulfur.bond_with(oxygen1, bond_type="double")
sulfur.bond_with(oxygen2, bond_type="double")
sulfur.bond_with(oxygen3, bond_type="single")
sulfur.bond_with(oxygen4, bond_type="single")
oxygen3.bond_with(hydrogen1, bond_type="single")
oxygen4.bond_with(hydrogen2, bond_type="single")

sulfuric_acid = Molecule(name="Sulfuric Acid")
for atom in [hydrogen1, hydrogen2, sulfur, oxygen1, oxygen2, oxygen3, oxygen4]:
    sulfuric_acid.add_atom(atom)

sulfuric_acid.summary()
print("\n")


# Example: Creating C6H6 (Benzene)
carbon_atoms = [Atom(atomic_number=6, neutrons=6) for _ in range(6)]
hydrogen_atoms = [Atom(atomic_number=1, neutrons=0) for _ in range(6)]

# Form the ring structure with alternating single and double bonds
for i in range(6):
    bond_type = "double" if i % 2 == 0 else "single"
    carbon_atoms[i].bond_with(carbon_atoms[(i + 1) % 6], bond_type=bond_type)
    carbon_atoms[i].bond_with(hydrogen_atoms[i], bond_type="single")

benzene = Molecule(name="Benzene")
for atom in carbon_atoms + hydrogen_atoms:
    benzene.add_atom(atom)

benzene.summary()
print("\n")


# Example: Creating H2O2 (Hydrogen Peroxide)
hydrogen1 = Atom(atomic_number=1, neutrons=0)
hydrogen2 = Atom(atomic_number=1, neutrons=0)
oxygen1 = Atom(atomic_number=8, neutrons=8)
oxygen2 = Atom(atomic_number=8, neutrons=8)

oxygen1.bond_with(oxygen2, bond_type="single")
oxygen1.bond_with(hydrogen1, bond_type="single")
oxygen2.bond_with(hydrogen2, bond_type="single")

hydrogen_peroxide = Molecule(name="Hydrogen Peroxide")
for atom in [hydrogen1, hydrogen2, oxygen1, oxygen2]:
    hydrogen_peroxide.add_atom(atom)

hydrogen_peroxide.summary()
