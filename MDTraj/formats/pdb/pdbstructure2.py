from mdtraj import element as elem
from mdtraj.core.topology import Topology, Chain, Residue, Atom


class PDBReader(object):
    def __init__(self, input_stream):
        self._cmodel = None
        self._cchain = None
        self._cresidue = None
        self._cpositions = []

        self._unit_cell_lengths = None
        self._unit_cell_angles = None
        self._topology = None

        self.f = f

    def read_model(self):

        for line in f:
            # Look for atoms
            if (line.find("ATOM  ") == 0) or (line.find("HETATM") == 0):
                self._add_atom(line)

            # Notice MODEL punctuation, for the next level of detail
            # in the structure->model->chain->residue->atom->position hierarchy
            elif (line.find("MODEL") == 0):
                self._add_model()

            elif (line.find("ENDMDL") == 0):
                self._end_model()
                break

            elif (line.find("END") == 0):
                self._end_model()
                break

            elif (line.find("TER") == 0 and line.split()[0] == "TER"):
                self._end_chain()

            elif (line.find("CRYST1") == 0):
                self._unit_cell_lengths = (float(line[6:15]), float(line[15:24]), float(line[24:33]))
                self._unit_cell_angles = (float(line[33:40]), float(line[40:47]), float(line[47:54]))

            elif (line.find("CONECT") == 0):
                atoms = [int(line[7:12])]
                for pos in (12, 17, 22, 27):
                    try:
                        atoms.append(int(line[pos:pos+5]))
                    except:
                        pass

    def _add_model(self):
        pass

    def _add_chain(self):
        pass

    def _add_residue(self):
        pass

    def _add_atom(self, line):
        record_name = line[0:6].strip()
        serial_number = int(line[6:11])
        name_with_spaces = line[12:16]
        name = name_with_spaces.strip()
        residue_name_with_spaces = line[17:20]

        # In some MD codes, notably ffamber in gromacs, residue name has a fourth character in
        # column 21
        possible_fourth_character = line[20:21]

        if possible_fourth_character != " ":
            # Fourth character should only be there if official 3 are already full
            if len(residue_name_with_spaces.strip()) != 3:
                raise ValueError('Misaligned residue name: %s' % line)
            residue_name_with_spaces += possible_fourth_character
        residue_name = residue_name_with_spaces.strip()

        #chain_id = line[21]
        residue_number = int(line[22:26])
        insertion_code = line[26]
        # coordinates, occupancy, and temperature factor belong in Atom.Location object
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        segment_id = line[72:76].strip()
        element_symbol = line[76:78].strip()

        
        if self._cmodel is None:
            self._cmodel = Topology()
        if self._cchain is None:
            self._cchain = self._cmodel.add_chain()
        if (self._cresidue is None) or (self._cresidue.resSeq != residue_number):
            self._cresidue = self._cmodel.add_residue(residue_name, self._cchain, residue_number)
        
        atom = self._cmodel.add_atom(name, self._guess_element(element_symbol), self._cresidue)
            
        self._cpositions.append([x, y, z])
        


    def _end_model(self):
        pass

    def _end_chain(self):
        pass

    def _end_residue(self):
        pass

    def _guess_element(self, symbol):
        
        upper = symbol.upper()
        if upper.startswith('CL'):
            element = elem.chlorine
        elif upper.startswith('NA'):
            element = elem.sodium
        elif upper.startswith('MG'):
            element = elem.magnesium
        elif upper.startswith('BE'):
            element = elem.beryllium
        elif upper.startswith('LI'):
            element = elem.lithium
        elif upper.startswith('K'):
            element = elem.potassium
        elif upper.startswith('ZN'):
            element = elem.zinc
        #elif len(residue) == 1 and upper.startswith('CA'):
        #    element = elem.calcium

        # TJL has edited this. There are a few issues here. First,
        # parsing for the element is non-trivial, so I do my best
        # below. Second, there is additional parsing code in
        # pdbstructure.py, and I am unsure why it doesn't get used
        # here...
        #elif len(residue) > 1 and upper.startswith('CE'):
        #    element = elem.carbon  # (probably) not Celenium...
        #elif len(residue) > 1 and upper.startswith('CD'):
        #    element = elem.carbon  # (probably) not Cadmium...
        #if residue.name in ['TRP', 'ARG', 'GLN', 'HIS'] and upper.startswith('NE'):
        #    element = elem.nitrogen  # (probably) not Neon...
        #if residue.name in ['ASN'] and upper.startswith('ND'):
        #    element = elem.nitrogen  # (probably) not ND...
        #elif residue.name == 'CYS' and upper.startswith('SG'):
        #    element = elem.sulfur  # (probably) not SG...
        else:
            try:
                element = elem.get_by_symbol(symbol[0])
            except KeyError:
                try:
                    symbol = symbol[0:2].strip().rstrip("AB0123456789").lstrip("0123456789")
                    element = elem.get_by_symbol(symbol)
                except KeyError:
                    element = None
        return element


if __name__ == '__main__':
    fn = '/home/rmcgibbo/local/mdtraj/MDTraj/testing/reference/2EQQ.pdb'
    with open(fn) as f:
        g = PDBReader(f)
        g.read_model()

        A, _ = g._cmodel.to_dataframe()

    import mdtraj as md
    B, _ = md.load(fn).topology.to_dataframe()

    diff = A['name'] != B['name']
    print A['name'][diff]
    print B['name'][diff]


