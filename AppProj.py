import requests
# atom, sequence, sheet and helix lists to be used are filled here.
def read_pdb_file(file_path):
    try:
        with open(file_path, 'r') as file:
            atom_data = []
            sequences_dict = {}
            helix_data = []
            sheet_data = []
            current_id = None
            current_sequence = ""
            for line in file:
                if line.startswith("ATOM"):
                    atom_type = line[0:6].strip()
                    atom_serial = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    residue_serial = int(line[22:26].strip())
                    x_coord = float(line[30:38].strip())
                    y_coord = float(line[38:46].strip())
                    z_coord = float(line[46:54].strip())

                    atom_data.append({
                        "atom_type": atom_type,
                        "atom_serial": atom_serial,
                        "atom_name": atom_name,
                        "residue_name": residue_name,
                        "chain_id": chain_id,
                        "residue_serial": residue_serial,
                        "x_coord": x_coord,
                        "y_coord": y_coord,
                        "z_coord": z_coord,
                    })
                if line.startswith("SEQRES"):
                    parts = line.split()
                    chain_id = parts[2]
                    if current_id is None or chain_id != current_id:
                        if current_id is not None:
                            sequences_dict[current_id] = current_sequence
                        current_id = chain_id
                        current_sequence = ""
                    current_sequence += "".join(parts[4:])

                if line.startswith("HELIX"):
                    helix_type = line[0:6].strip()
                    helix_serial = int(line[7:11].strip())
                    init_residue_name = line[15:18].strip()
                    init_chain_id = line[19].strip()
                    init_residue_serial = int(line[21:25].strip())
                    end_residue_name = line[27:30].strip()
                    end_chain_id = line[31].strip()
                    end_residue_serial = int(line[33:37].strip())
                    helix_length = line[71:76].strip()

                    helix_data.append({
                        "helix_type": helix_type,
                        "helix_serial": helix_serial,
                        "init_residue_name": init_residue_name,
                        "init_chain_id": init_chain_id,
                        "init_residue_serial": init_residue_serial,
                        "end_residue_name": end_residue_name,
                        "end_chain_id": end_chain_id,
                        "end_residue_serial": end_residue_serial,
                        "helix_length": helix_length,
                    })

                if line.startswith("SHEET"):
                    sheet_type = line[0:6].strip()
                    strand = int(line[7:10].strip())
                    sheet_id = line[11:14].strip()
                    num_strands = line[14:16].strip()
                    init_residue_name = line[17:20].strip()
                    init_chain_id = line[21].strip()
                    init_residue_serial = int(line[22:26].strip())
                    end_residue_name = line[28:31].strip()
                    end_chain_id = line[32].strip()
                    end_residue_serial = int(line[33:37].strip())
                    sense_strand = line[38:40].strip()

                    sheet_data.append({
                        "sheet_type": sheet_type,
                        "strand": strand,
                        "sheet_id": sheet_id,
                        "num_strands": num_strands,
                        "init_residue_name": init_residue_name,
                        "init_chain_id": init_chain_id,
                        "init_residue_serial": init_residue_serial,
                        "end_residue_name": end_residue_name,
                        "end_chain_id": end_chain_id,
                        "end_residue_serial": end_residue_serial,
                        "sense_strand": sense_strand,
                    })
        
        return atom_data,  helix_data, sheet_data,sequences_dict
    except FileNotFoundError as e:
        print(f"Error: {e}")

    except PermissionError as e:
        print(f"Error: Permission denied. {e}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    
def count_atoms_in_amino_acid(atom_data, chain_id, aa_sequence_number):
    atom_count = 0

    for atom in atom_data:
        if atom["chain_id"] == chain_id and atom["residue_serial"] == aa_sequence_number:
            atom_count += 1

    return atom_count

def three_to_one_amino_acid(three_letter_code):
    amino_acid_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return amino_acid_map.get(three_letter_code, three_letter_code)

def convert_sequence_three_to_one(sequence):
    return ''.join(three_to_one_amino_acid(sequence[i:i+3]) for i in range(0, len(sequence), 3))

def count_side_chain_atoms(atom_data, chain_id, aa_sequence_number):
    side_chain_count = 0
    
    # Define the list of atom names to exclude
    exclusion_list = ['N', 'CA', 'C', 'O']
    
    for atom in atom_data:
        if atom["chain_id"] == chain_id and atom["residue_serial"] == aa_sequence_number and atom["atom_name"] not in exclusion_list:
            side_chain_count += 1

    return side_chain_count

def count_amino_acids_in_chain(atom_data, chain_id):
    amino_acids = set()

    for atom in atom_data:
        if atom["chain_id"] == chain_id:
            amino_acids.add((atom["residue_name"], atom["residue_serial"]))

    return len(amino_acids)

def count_helices_sheets():
    strand_per_chain = {}
    # Count helices per chain
    for strand in sheet_data:
        chain_id = strand["init_chain_id"]
        if chain_id in strand_per_chain:
            strand_per_chain[chain_id] += 1
        else:
            strand_per_chain[chain_id] = 1
    
    
    helices_per_chain = {}
    # Count helices per chain
    for helix in helix_data:
        chain_id = helix["init_chain_id"]
        if chain_id in helices_per_chain:
            helices_per_chain[chain_id] += 1
        else:
            helices_per_chain[chain_id] = 1
    return  helices_per_chain, strand_per_chain

def find_atom_coordinates(atom_data, chain_id, residue_number, atom_name):
    for atom in atom_data:
        if (atom["chain_id"] == chain_id and atom["residue_serial"] == residue_number and atom["atom_name"] == atom_name):
            return (
                atom["x_coord"],
                atom["y_coord"],
                atom["z_coord"]
            )

    # If the atom is not found, return None or an appropriate message
    return None

def find_amino_acids_for_helix(helix_data, helix_number):
    for helix in helix_data:
        if helix["helix_serial"] == helix_number:
            return helix["helix_length"], helix["init_residue_name"], helix["init_residue_serial"], helix["end_residue_name"], helix["end_residue_serial"]

    # If the helix number is not found, return None or an appropriate message
    return None

def calculate_distance(point1, point2):
    squared_sum = sum((a - b) ** 2 for a, b in zip(point1, point2))
    distance = squared_sum ** 0.5
    return distance

def find_distance_between_atoms(atom_data, chain_id1, residue_number1, atom_name1, chain_id2, residue_number2, atom_name2):
    atom1 = None
    atom2 = None

    for atom in atom_data:
        if (
            atom["chain_id"] == chain_id1
            and atom["residue_serial"] == residue_number1
            and atom["atom_name"] == atom_name1
        ):
            atom1 = (atom["x_coord"], atom["y_coord"], atom["z_coord"])
        elif (
            atom["chain_id"] == chain_id2
            and atom["residue_serial"] == residue_number2
            and atom["atom_name"] == atom_name2
        ):
            atom2 = (atom["x_coord"], atom["y_coord"], atom["z_coord"])

        if atom1 and atom2:
            break

    if atom1 is None or atom2 is None:
        return None  # If any of the atoms is not found

    # Calculate the Euclidean distance between the two atoms
    distance = calculate_distance(atom1, atom2)

    return distance
# To read sequences.txt that contains a list of IDs for proteins and a target chain
def read_sequences(filename):
    try:
        entries = []
        with open(filename, "r") as file:
            lines = file.readlines()
            for line in lines:
                protein_ID = line[0:4].strip()
                chain_ID = line[12:13].strip()
                sequence = line[14:].strip()

                entries.append({
                    "protein_ID": protein_ID,
                    "chain_ID": chain_ID,
                    "sequence": sequence,
                })
        return entries

    except FileNotFoundError as e:
        print(f"Error: {e}")

    except PermissionError as e:
        print(f"Error: Permission denied. {e}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    
def pairwise_alignment(seq1, seq2, gap_penalty=-2):
    # BLOSUM62 scoring matrix
    blosum62 = {
        ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
        ('R', 'A'): -1, ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
        ('N', 'A'): -2, ('N', 'R'): 0, ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
        ('D', 'A'): -2, ('D', 'R'): -2, ('D', 'N'): 1, ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
        ('C', 'A'): 0, ('C', 'R'): -3, ('C', 'N'): -3, ('C', 'D'): -3, ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
        ('Q', 'A'): -1, ('Q', 'R'): 1, ('Q', 'N'): 0, ('Q', 'D'): 0, ('Q', 'C'): -3, ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
        ('E', 'A'): -1, ('E', 'R'): 0, ('E', 'N'): 0, ('E', 'D'): 2, ('E', 'C'): -4, ('E', 'Q'): 2, ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
        ('G', 'A'): 0, ('G', 'R'): -2, ('G', 'N'): 0, ('G', 'D'): -1, ('G', 'C'): -3, ('G', 'Q'): -2, ('G', 'E'): -2, ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
        ('H', 'A'): -2, ('H', 'R'): 0, ('H', 'N'): 1, ('H', 'D'): -1, ('H', 'C'): -3, ('H', 'Q'): 0, ('H', 'E'): 0, ('H', 'G'): -2, ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
        ('I', 'A'): -1, ('I', 'R'): -3, ('I', 'N'): -3, ('I', 'D'): -3, ('I', 'C'): -1, ('I', 'Q'): -3, ('I', 'E'): -3, ('I', 'G'): -4, ('I', 'H'): -3, ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
        ('L', 'A'): -1, ('L', 'R'): -2, ('L', 'N'): -3, ('L', 'D'): -4, ('L', 'C'): -1, ('L', 'Q'): -2, ('L', 'E'): -3, ('L', 'G'): -4, ('L', 'H'): -3, ('L', 'I'): 2, ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
        ('K', 'A'): -1, ('K', 'R'): 2, ('K', 'N'): 0, ('K', 'D'): -1, ('K', 'C'): -3, ('K', 'Q'): 1, ('K', 'E'): 1, ('K', 'G'): -2, ('K', 'H'): -1, ('K', 'I'): -3, ('K', 'L'): -2, ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
        ('M', 'A'): -1, ('M', 'R'): -1, ('M', 'N'): -2, ('M', 'D'): -3, ('M', 'C'): -1, ('M', 'Q'): 0, ('M', 'E'): -2, ('M', 'G'): -3, ('M', 'H'): -2, ('M', 'I'): 1, ('M', 'L'): 2, ('M', 'K'): -1, ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
        ('F', 'A'): -2, ('F', 'R'): -3, ('F', 'N'): -3, ('F', 'D'):-3, ('F', 'C'): -2, ('F', 'Q'): -3, ('F', 'E'): -3, ('F', 'G'): -3, ('F', 'H'): -1, ('F', 'I'): 0, ('F', 'L'): 0, ('F', 'K'): -3, ('F', 'M'): 0, ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
        ('P', 'A'): -1, ('P', 'R'): -2, ('P', 'N'): -2, ('P', 'D'): -1, ('P', 'C'): -3, ('P', 'Q'): -1, ('P', 'E'): -1, ('P', 'G'): -2, ('P', 'H'): -2, ('P', 'I'): -3, ('P', 'L'): -3, ('P', 'K'): -1, ('P', 'M'): -2, ('P', 'F'): -4, ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
        ('S', 'A'): 1, ('S', 'R'): -1, ('S', 'N'): 1, ('S', 'D'): 0, ('S', 'C'): -1, ('S', 'Q'): 0, ('S', 'E'): 0, ('S', 'G'): 0, ('S', 'H'): -1, ('S', 'I'): -2, ('S', 'L'): -2, ('S', 'K'): 0, ('S', 'M'): -1, ('S', 'F'): -2, ('S', 'P'): -1, ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
        ('T', 'A'): 0, ('T', 'R'): -1, ('T', 'N'): 0, ('T', 'D'): -1, ('T', 'C'): -1, ('T', 'Q'): -1, ('T', 'E'): -1, ('T', 'G'): -2, ('T', 'H'): -2, ('T', 'I'): -1, ('T', 'L'): -1, ('T', 'K'): -1, ('T', 'M'): -1, ('T', 'F'): -2, ('T', 'P'): -1, ('T', 'S'): 1, ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
        ('W', 'A'): -3, ('W', 'R'): -3, ('W', 'N'): -4, ('W', 'D'): -4, ('W', 'C'): -2, ('W', 'Q'): -2, ('W', 'E'): -3, ('W', 'G'): -2, ('W', 'H'): -2, ('W', 'I'): -3, ('W', 'L'): -2, ('W', 'K'): -3, ('W', 'M'): -1, ('W', 'F'): 1, ('W', 'P'): -4, ('W', 'S'): -3, ('W', 'T'): -2, ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
        ('Y', 'A'): -2, ('Y', 'R'): -2, ('Y', 'N'): -2, ('Y', 'D'): -3, ('Y', 'C'): -2, ('Y', 'Q'): -1, ('Y', 'E'): -2, ('Y', 'G'): -3, ('Y', 'H'): 2, ('Y', 'I'): -1, ('Y', 'L'): -1, ('Y', 'K'): -2, ('Y', 'M'): -1, ('Y', 'F'): 3, ('Y', 'P'): -3, ('Y', 'S'): -2, ('Y', 'T'): -2, ('Y', 'W'): 2, ('Y', 'Y'): 7, ('Y', 'V'): -1,
        ('V', 'A'): 0, ('V', 'R'): -3, ('V', 'N'): -3, ('V', 'D'): -3, ('V', 'C'): -1, ('V', 'Q'): -2, ('V', 'E'): -2, ('V', 'G'): -3, ('V', 'H'): -3, ('V', 'I'): 3, ('V', 'L'): 1, ('V', 'K'): -2, ('V', 'M'): 1, ('V', 'F'): -1, ('V', 'P'): -2, ('V', 'S'): -2, ('V', 'T'): 0, ('V', 'W'): -3, ('V', 'Y'): -1, ('V', 'V'): 4,

    }

    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    # Initialize the score matrix
    score_matrix = [[0] * (len_seq2 + 1) for _ in range(len_seq1 + 1)]

    # Initialize the traceback matrix for backtracking
    traceback_matrix = [[0] * (len_seq2 + 1) for _ in range(len_seq1 + 1)]

    # Initialize the first row and column with gap penalties
    for i in range(1, len_seq1 + 1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
        traceback_matrix[i][0] = 1  # 1 represents an upward move (gap in seq2)

    for j in range(1, len_seq2 + 1):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty
        traceback_matrix[0][j] = 2  # 2 represents a leftward move (gap in seq1)

    # Fill in the score matrix using the BLOSUM62 scores
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            match_score = score_matrix[i-1][j-1] + blosum62.get((seq1[i-1], seq2[j-1]), 1)
            gap_penalty_seq1 = score_matrix[i-1][j] + gap_penalty
            gap_penalty_seq2 = score_matrix[i][j-1] + gap_penalty

            # Determine the maximum score for the current cell
            max_score = max(match_score, gap_penalty_seq1, gap_penalty_seq2)

            score_matrix[i][j] = max_score

            # Update the traceback matrix
            if max_score == match_score:
                traceback_matrix[i][j] = 0  # 0 represents a diagonal move (match/mismatch)
            elif max_score == gap_penalty_seq1:
                traceback_matrix[i][j] = 1  # 1 represents an upward move (gap in seq2)
            else:
                traceback_matrix[i][j] = 2  # 2 represents a leftward move (gap in seq1)

    # Backtrack to find the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = len_seq1, len_seq2

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 0:  # Diagonal move
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 1:  # Upward move (gap in seq2)
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:  # Leftward move (gap in seq1)
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    # Calculate alignment score
    alignment_score = score_matrix[len_seq1][len_seq2]

    # Calculate identity percentage
    identity_count = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2))
    identity_percentage = (identity_count / len(aligned_seq1)) * 100

    return aligned_seq1, aligned_seq2, alignment_score, identity_percentage

# Menu function to display options
def menu():
    print("Options:")
    print("1. Show the number of atoms in a particular Amino acid in a given chain")
    print("2. Return the sequence of the amino acids make up the protein including amino acids missing the structural information")
    print("3. Show the number of atoms in the side chain for particular aa in a given chain")
    print("4. Show the number of amino acids in a given chain")
    print("5. Show the number of helices or sheets in the chain")
    print("6. Return the coordinates of a particular atom")
    print("7. Show the name and number of the amino acid at the beginning or end of a given helix by its number")
    print("8. Find the distance between any two atoms in the chain")
    print("9. Report the alignment score for the chain of the protein against all chains in the list and identity percentage")
    print("0. Quit")


def download_pdb_file():
    pdb_id = input("Enter the ID of the protein to download the pdb file ")
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check for errors
    except requests.exceptions.HTTPError as errh:
        print(f"HTTP Error: {errh}")
    except requests.exceptions.ConnectionError as errc:
        print(f"Error Connecting: {errc}")
    except requests.exceptions.Timeout as errt:
        print(f"Timeout Error: {errt}")
    except requests.exceptions.RequestException as err:
        print(f"Something went wrong: {err}")
    
    # Save the content to a file
    with open('file.pdb', 'wb') as f:
        f.write(response.content)



def select_file_menu():  
    print("How do you want to upload the pdb file?")
    print("Options:")
    print("1. Reads a pdb file from the personal computer (1kam, 1bmc, 3uts)")
    print("2. Fetch a pdb file from online servers like www.rcsb.org")

# Main program
select_file_menu()
choice = input("Enter your choice (1/2): ")

if choice == '1':
    local_pdb_file = input("Please type the name of the pdb file tha saved in the personal computer: ") + ".pdb"
elif choice == '2':   
    download_pdb_file()  
    local_pdb_file = "file.pdb"
    
# Fill lists according to selected pdb file
atom_data,  helix_data, sheet_data,sequences_dict = read_pdb_file(local_pdb_file)
while True:
    menu()
    choice = input("Enter your choice (1/2/3/4/5/6/7/8/9/0): ")
    print('\n')
    if choice == '0':
        break  # Exit the loop if '0' is selected

    if choice in ('1', '2', '3', '4', '5', '6', '7', '8', '9'):
        try:
            
            if choice == '1':
                chain_id = input("Please type the chain you would like to work on: ").upper()
                aa_sequence_number = int(input("Please type the amino acid's sequence number you would like to work on: : "))
                atom_count = count_atoms_in_amino_acid(atom_data, chain_id, aa_sequence_number)
                print("\nThe number of atoms : ", atom_count)
            elif choice == '2':
                chain_id = input("Please type the chain you would like to work on: ").upper()
                print("The sequence of the amino acids (1-letter code) in Chain ", chain_id)
                print(convert_sequence_three_to_one(sequences_dict[chain_id]))
            elif choice == '3':
                chain_id = input("Please type the chain you would like to work on: ").upper()
                aa_sequence_number = int(input("Please type the amino acid's sequence number you would like to work on: : "))
                side_chain_count = count_side_chain_atoms(atom_data, chain_id, aa_sequence_number)
                print("\nThe number of atoms in the side chain : ", side_chain_count)
            elif choice == '4':
                chain_id = input("Please type the chain you would like to work on: ").upper()
                aa_count = count_amino_acids_in_chain(atom_data, chain_id)
                print("\nThe number of amino acids : ", aa_count)
            elif choice == '5':
                helices_per_chain, strand_per_chain = count_helices_sheets()
                print("\nNumber of Strand per Chain:")
                for chain, count in strand_per_chain.items():
                    print(f"Chain {chain}: {count} strand")
                print("\nNumber of Helices per Chain:")
                for chain, count in helices_per_chain.items():
                    print(f"Chain {chain}: {count} helices")
            elif choice == '6':
                chain_id = input("Please type the chain you would like to work on: ").upper()
                residue_number = int(input("Please type the amino acid's sequence number you would like to work on: : "))
                atom_name = input("Please type atom name you would like to work on: ").upper()
                coordinates = find_atom_coordinates(atom_data, chain_id, residue_number, atom_name)

                if coordinates is not None:
                    x, y, z = coordinates
                    print(f"\nCoordinates of Atom {atom_name} in Chain {chain_id}, Residue {residue_number}: X={x}, Y={y}, Z={z}")
                else:
                    print(f"\nAtom {atom_name} in Chain {chain_id}, Residue {residue_number} not found.")
            elif choice == '7':
                helix_number = int(input("Please type the helix number you would like to work on: "))
                helix_length, init_amino_acid_name, init_amino_acid_number, end_amino_acid_name, end_amino_acid_number = find_amino_acids_for_helix(helix_data, helix_number)

                if init_amino_acid_name is not None and init_amino_acid_number is not None:
                    print(f"\nInitial Amino Acid at Helix {helix_number}: Name={init_amino_acid_name}, Number={init_amino_acid_number}")
                    print(f"\nTerminal Amino Acid at Helix {helix_number}: Name={end_amino_acid_name}, Number={end_amino_acid_number}")
                    print(f"\nHelix Length {helix_length}")
                else:
                    print(f"Helix {helix_number} not found.")
            elif choice == '8':
                chain_id1 = input("Please type the chain of first atom you would like to work on: ").upper()
                residue_number1 = int(input("Please type the amino acid's sequence number of first atom you would like to work on: : "))
                atom_name1 =  input("Please type first atom you would like to work on: ").upper()
                chain_id2 = input("Please type the chain of second atom you would like to work on: ").upper()
                residue_number2 = int(input("Please type the amino acid's sequence number of second atom you would like to work on: : "))
                atom_name2 = input("Please type second atom you would like to work on: ").upper()

                distance = find_distance_between_atoms(atom_data, chain_id1, residue_number1, atom_name1, chain_id2, residue_number2, atom_name2)

                if distance is not None:
                    print(f"\nDistance between Atom {atom_name1} in Residue {residue_number1} in Chain {chain_id1} and Atom {atom_name2} in Residue {residue_number2} in Chain {chain_id2}: {distance:.2f} Ångstroms")
                else:
                    print(f"\nOne or both atoms not found.")
            elif choice == '9':
                chain_id = input("Please type the chain you would like to work on: ").upper()
                target_chain = convert_sequence_three_to_one(sequences_dict[chain_id])
                sequence_list = read_sequences('sequences.txt')
                for entry in sequence_list:
                    protein_id = entry['protein_ID']
                    chain_id = entry['chain_ID']
                    sequence = entry['sequence']
                    
                    print("In sequences.txt Protein ",protein_id, " Chain ", chain_id)
                    aligned_seq1, aligned_seq2, alignment_score, identity_percentage = pairwise_alignment(target_chain, sequence)
                    #print("Aligned Sequence 1:", aligned_seq1)
                    #print("Aligned Sequence 2:", aligned_seq2)
                    print("Alignment Score:", alignment_score)
                    print("Identity Percentage:", identity_percentage) 
                    print("\n")
                    
        except ValueError:
            print("Invalid input. Please enter valid numbers.")
    else:
        print("Invalid input. Please choose a valid option (1/2/3/4/5/6/7/8/9/0).")
