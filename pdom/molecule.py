# -*- coding: utf-8 -*-
import urllib.request
from pathlib import Path
from textwrap import dedent

import numpy as np
import pubchempy as pcpug
from appdirs import user_cache_dir
from skimage import draw
from skimage.morphology import binary_erosion
from skimage.morphology import disk as morph_disk

from pdom import Parameter
from pdom import module_info


class Molecule(object):
    """This class should not be initialised directly.
    Use one of the following class methods instead:
    :meth:`from_chem_id`, :meth:`from_folder`,
    :meth:`from_inchi`, :meth:`from_inchi_key`,
    :meth:`from_iupac_name`, :meth:`from_name`

    :param identifier:
        * ``name`` common name
        * ``iupac_name`` IUPAC name (can be the same as common name)
        * ``chem_id`` compound ID from `PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_
        * ``inchi`` IUPAC International Chemical Identifier (human readable)
        * ``inchi_key`` IUPAC International Chemical Identifier (hash)
    :type identifier: dict
    :param properties:
        * ``chem_formula`` chemical formula e.g. 'C10H22'
        * ``chem_formdic`` chemical formula as dict
        * ``mol_weight`` molecular weight in g/mol
        * ``mol_volume`` molecular volume in cm^3/mol
        * ``mol_surface`` largest projected surface area m^2/molecule
        * ``excess_bonds`` number of bonds in excess of a simple carbon chain
        * ``structure_3d`` atom position list
    :type properties: dict
    :param save: save data to the user cache
    :type save: bool
    """

    db_folder = Path(user_cache_dir(module_info.app_name, module_info.app_group_id))

    def __init__(self, identifier, properties, save=True):
        """Constructor method
        """
        #: *dict* of identifiers as described in :class:`Molecule`
        self.identifier = identifier

        #: *dict* of properties as described in :class:`Molecule`
        self.properties = properties

        if self.properties['mol_surface'] is None:
            self._set_mol_surface()
        if self.properties['excess_bonds'] is None:
            self._set_excess_bonds()
        if save:
            self.save()

    @classmethod
    def from_folder(cls, folder):
        """create :class:`Molecule` instance from a folder created by :class:`save`

        :param folder: molecule folder
        :type folder: str, folder
        :return: :class:`Molecule` instance
        """

        identifier = {}
        save_identifier = (folder / 'identifier.txt').read_text()
        for line in save_identifier.splitlines():
            keydata = line.split(':')
            key = keydata[0].strip()
            data = keydata[1].strip()
            identifier[key] = data

        properties = {}
        save_properties = (folder / 'properties.txt').read_text()
        for line in save_properties.splitlines():
            keydata = line.split(':')
            key = keydata[0].strip()
            data = keydata[1].strip()
            if key in ['mol_surface', 'mol_weight', 'mol_volume']:
                properties[key] = float(data)
            if key in ['excess_bonds']:
                properties[key] = int(data)
            elif key in ['chem_formula']:
                properties[key] = data

        properties['structure_3d'] = []
        properties['chem_formdic'] = {}
        save_structure_3d = (folder / 'structure.xyz').read_text().splitlines()
        # Skip two header lines

        for atom_raw in save_structure_3d[2:]:
            atom = atom_raw.split()
            element = atom[0].title()
            properties['structure_3d'].append([element, float(atom[1]), float(atom[2]), float(atom[3])])
            if element in properties['chem_formdic']:
                properties['chem_formdic'][element] += 1
            else:
                properties['chem_formdic'][element] = 1

        return cls(identifier, properties, save=False)

    @classmethod
    def from_chem_id(cls, chem_id, name=None, inchi_key=None):
        """Create :class:`Molecule` instance identified by a chemID from `PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_ data.

        :note: if molecule is cached in ``Molecule.db_folder`` load from there

        :param chem_id: compound ID from `PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_
        :type chem_id: int
        :param name: overwrite compound name
        :type name: str, optional
        :param inchi_key: IUPAC International Chemical Identifier (to check cache quickly)
        :type inchi_key: str, optional
        :return: :class:`Molecule` instance
        """
        # If data is already stored on the system use it
        if inchi_key is None:
            received_data = pcpug.get_properties('InChIKey', chem_id, 'cid')
            if len(received_data) == 1:
                chem_id = int(received_data[0]['CID'])
                inchi_key = received_data[0]['InChIKey']
        search_folder = cls.db_folder / inchi_key
        if search_folder.is_dir():
            return cls.from_folder(search_folder)

        # Collect molecule information from the PubChem database or calculate them
        pcmol2d = pcpug.get_compounds(chem_id, record_type='2d')[0]

        # Get all identifier
        identifier = dict(chem_id=str(chem_id))
        try:
            identifier['iupac_name'] = str(pcmol2d.iupac_name)
        except:
            identifier['iupac_name'] = "-"

        identifier['inchi'] = str(pcmol2d.inchi)
        identifier['inchi_key'] = str(pcmol2d.inchikey)
        if name:
            identifier['name'] = name.lower()
        else:
            try:
                identifier['name'] = pcpug.get_synonyms(chem_id, 'cid')[0]['Synonym'][0].lower()
            except:
                identifier['name'] = pcmol2d.iupac_name.lower()

        # Get all properties
        properties = dict()
        properties['mol_weight'] = pcmol2d.molecular_weight

        properties['chem_formula'] = pcmol2d.molecular_formula
        properties['structure_3d'] = []
        properties['chem_formdic'] = {}

        try:
            pcmol3d = pcpug.get_compounds(chem_id, record_type='3d')[0]
            properties['mol_volume'] = pcmol3d.volume_3d
            for atom in pcmol3d.atoms:
                element = atom.element.title()
                properties['structure_3d'].append([element, float(atom.x), float(atom.y), float(atom.z)])
                if element in properties['chem_formdic']:
                    properties['chem_formdic'][element] += 1
                else:
                    properties['chem_formdic'][element] = 1

        except:
            properties['mol_volume'] = 0.0
            folder = cls.db_folder / 'helper_files' / '3d' / identifier['inchi_key']
            folder.mkdir(exist_ok=True, parents=True)
            try:
                save_structure_3d_file = (folder/'structure.xyz').open('r')
                next(save_structure_3d_file)
                next(save_structure_3d_file)
                for atom_raw in save_structure_3d_file:
                    atom = atom_raw.split()
                    element = atom[0].title()
                    properties['structure_3d'].append([element, float(atom[1]), float(atom[2]), float(atom[3])])
                    if element in properties['chem_formdic']:
                        properties['chem_formdic'][element] += 1
                    else:
                        properties['chem_formdic'][element] = 1
            except FileNotFoundError:
                folder2d = cls.db_folder / 'helper_files' / '2d' / identifier['inchi_key']
                folder2d.mkdir(exist_ok=True, parents=True)
                pcpug.download('SDF', str((folder2d/'structure.sdf').absolute()), chem_id, 'cid', overwrite=True)
                print(f'No 3D data found. Please convert {str((folder2d/"structure.sdf").absolute())} to a 3D .xyz file.')
                print(f'The result should be stored under {str((folder/"structure.xyz").absolute())}')
                return None

        properties['mol_surface'] = None
        properties['excess_bonds'] = None
        return cls(identifier, properties)

    @classmethod
    def from_name(cls, name):
        """Create :class:`Molecule` instance identified by a name if unique

        :note: queries ``chem_id`` and calls :meth:`from_chem_id`

        :param name: unique compound name
        :type name: str
        :return: :class:`Molecule` instance
        """
        received_data = pcpug.get_properties('IUPACName', name, 'name')
        if len(received_data) == 1:
            chem_id = received_data[0]['CID']
            return cls.from_chem_id(chem_id, name=name)

    def structure_3d(self, rotated=True):
        """Returned 3d structure of the molecule

        :param rotated: if ``True`` the structure is rotated to cover the maximum surface in the xy space
        :type rotated: bool
        :return: 3D structure (symbol, x, y, z)
        :rtype: list
        """
        if rotated:
            return self._rotate_structure(self.properties['structure_3d'])
        else:
            return self.properties['structure_3d']

    def save(self, folder=None, name=None):
        """saves the molecule information to disk, can be loaded with :meth:`from_folder`

        :note: If called without parameters the molecule is saved in an appropriate cache folder ``Molecule.db_folder``.

        :param folder: parent folder
        :type folder: str, Path
        :param name: molecule folder name (default INCHI key)
        :type name:  str
        """
        # Save all molecule data in a folder for later use

        # Select place to save
        if folder:
            folder = Path(folder)
            if name:
                folder = folder / name
            else:
                folder = folder / self.identifier['inchi_key']
        else:
            folder = self.db_folder / self.identifier['inchi_key']
        folder.mkdir(exist_ok=True, parents=True)

        # Write the data to separate files (human readable)

        identifier_text = '\n'.join([f"{key}: {str(value)}" for key, value in self.identifier.items()])
        (folder / 'identifier.txt').write_text(identifier_text)

        properties_text = '\n'.join([f"{key}: {str(value)}" for key, value
                                     in self.properties.items() if key not in ['chem_formdic', 'structure_3d']])
        (folder / 'properties.txt').write_text(properties_text)

        structure_3d_text = "%i\n" % (len(self.structure_3d()))
        structure_3d_text += "Structure of %s generate by pdom\n" % (len(self.identifier['name']))
        for atom in self.structure_3d():
            structure_3d_text += "%s %12.5f %12.5f %12.5f\n" % (atom[0], atom[1], atom[2], atom[3])
        (folder / 'structure.xyz').write_text(structure_3d_text)
        info = dedent("""\
            directory contains following data

            identifier
            name         - Common name 
            iupac_name   - IUPAC name (can be the same as common name)
            chem_id      - Compound ID from PubChem (http://pubchem.ncbi.nlm.nih.gov/)
            inchi        - IUPAC International Chemical Identifier (human readable)
            inchi_key    - IUPAC International Chemical Identifier (key for hash)

            properties
            chem_formula - Chemical formula e.g. 'C10H22'
            mol_weight    - Molecular weight in g/mol
            mol_volume   - Molecular volume in cm^3/mol
            mol_surface  - Largest projected surface area m^2/molecule
            excess_bonds - Number of bonds in excess of a simple carbon chain

            structure
            """)
        (folder / 'info.txt').write_text(info)

    @staticmethod
    def _rotate_structure(structure_3d):
        coords = []
        for atom in structure_3d:
            coords.append(atom[1:])
        coords = np.array(coords)
        coords -= coords.mean(axis=0)
        cov = coords.T.dot(coords)
        eigen_values, eigen_vectors = np.linalg.eigh(cov)
        eigen_vectors = eigen_vectors[:, eigen_values.argsort()[::-1]]
        new_structure_3d = []
        for atom in structure_3d:
            new_coordinates = np.dot(np.array(np.array(atom[1:])), eigen_vectors)
            new_structure_3d.append([atom[0], new_coordinates[0], new_coordinates[1], new_coordinates[2]])
        return new_structure_3d

    def _set_excess_bonds(self):
        start_row = 3
        sdf_link = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{self.identifier["chem_id"]}/SDF'
        sdf_raw = urllib.request.urlopen(sdf_link)
        sdf = sdf_raw.read().decode()
        sdf_raw.close()
        mol_data = sdf.partition('M  END')[0].splitlines()

        mol_meta = mol_data[start_row].split()
        atom_count = int(mol_meta[0])
        bond_count = int(mol_meta[1])

        carbon = []
        hetero_atoms = []
        hetero_count = {}
        backbone_bonds = 0
        carbon_count = 0

        for k, i in enumerate(range(start_row + 1, start_row + atom_count + 1)):
            atom = mol_data[i].split()
            if atom[3] == 'C':
                carbon.append(k + 1)
                carbon_count += 1

        for i in range(start_row + atom_count + 1, start_row + atom_count + bond_count + 1):
            bond = mol_data[i].split()
            atom_a = int(bond[0])
            atom_b = int(bond[1])

            if atom_a in carbon:
                if atom_b not in carbon:
                    try:
                        hetero_count[atom_b] += 1
                    except KeyError:
                        hetero_count[atom_b] = 1
            elif atom_b in carbon:
                try:
                    hetero_count[atom_a] += 1
                except KeyError:
                    hetero_count[atom_a] = 1

        for hetero_atom in hetero_count:
            if hetero_count[hetero_atom] > 1:
                hetero_atoms.append(hetero_atom)

        for i in range(start_row + atom_count + 1, start_row + atom_count + bond_count + 1):
            bond = mol_data[i].split()
            atom_a = int(bond[0])
            atom_b = int(bond[1])
            a_b_count = int(bond[2])

            if atom_a in carbon:
                if atom_b in (carbon + hetero_atoms):
                    backbone_bonds += a_b_count
            elif atom_a in hetero_atoms:
                if atom_b in carbon:
                    backbone_bonds += a_b_count

        self.properties['excess_bonds'] = backbone_bonds - self.properties['chem_formdic']['C'] + 1

    def _set_mol_surface(self):
        probe_size = 1.4  # radius water Ang
        scale = 30  # pixel per Ang
        cx = []
        cy = []
        cr = []
        for atom in self.structure_3d():
            cx.append(atom[1])
            cy.append(atom[2])
            cr.append(Parameter.van_der_waals_radii[atom[0]])
        x_scale = int(np.ceil(max(cx) - min(cx) + 2*max(cr) + 2*probe_size)*scale)
        x_offset = min(cx) - max(cr) - probe_size
        y_scale = int(np.ceil(max(cy) - min(cy) + 2*max(cr) + 2*probe_size)*scale)
        y_offset = min(cy) - max(cr) - probe_size
        img = np.zeros((x_scale, y_scale))
        for i in range(len(cx)):
            rr, cc = draw.disk(center=(int(round((cx[i] - x_offset) * scale)), int(round((cy[i] - y_offset) * scale))),
                               radius=int(round((cr[i]+probe_size)*scale)))
            img[rr, cc] = 1
        probe = morph_disk(int(round(probe_size * scale)))
        img = binary_erosion(img, probe)
        surface = np.sum(img)/(scale**2)  # Ang^2
        self.properties['mol_surface'] = surface * 1e-20  # m^2

    def __del__(self):
        pass
