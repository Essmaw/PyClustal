"""This script is the main programm of PyClustal

Usage:
======
    python src/pyclustal.py --f [fasta_file_path] --seq-type [seq_type] --sub-matrix [sub_matrix] --gap-open [gap_open] --gap-ext [gap_ext] --job-name [job_name] --tag-log [tag_log]

Arguments:
==========
    --f : str
        The path to the fasta file containing the sequences to align.
    --seq-type : str (optional)
        The type of sequences to align. It can be either "dna" or "protein". The default value is stored in DEFAULT_SEQ_TYPE.
    --sub-matrix : str (optional)
        The substitution matrix to use for the alignment.The default value is stored in DEFAULT_SUB_MATRIX.
    --gap-open : int (optional)
        The gap opening penalty. The default value is stored in DEFAULT_GAP_OPEN.
    --gap-ext : int (optional)
        The gap extension penalty. The default value is stored in DEFAULT_GAP_EXT.
    --job-name : str (optional)
        The name of the output fasta file containing the aligned sequences. The default value is name of the input file with "_aligned" appended.
    --tag-log : bool (optional)
        The flag to enable or disable the logging of the alignment process. The default value is stored in DEFAULT_TAG_LOG.

Example:
========
    python src/pyclustal.py --f data/sequences.fasta --seq-type protein --sub-matrix BLOSOM62 --gap-open -5 --gap-ext -1 --job-name aligned_sequences.fasta --tag-log False

This command will align the sequences contained in the file "data/sequences.fasta" using the blosum62 substitution matrix, a gap opening penalty of 10 and a gap extension penalty of 1. 
The aligned sequences will be saved in the file "aligned_sequences.fasta" in the results folder. And the alignment process of each pair of sequences will not be logged.
"""

# METADATA
__authors__ = "Essmay Touami"
__contact__ = "essmay.touami@etu.u-paris.fr"
__date__ = "september 2024"
__version__ = "1.0.0"


# LIBRARY IMPORTS
import os
import sys
import argparse
import numpy as np
from tqdm import tqdm
from typing import Tuple, Dict, Union, List

import re
from ete3 import Tree
from Bio import SeqIO
from loguru import logger
from tabulate import tabulate
from Bio.Align import substitution_matrices
from concurrent.futures import ProcessPoolExecutor, as_completed


# CONSTANTS
DEFAULT_SEQ_TYPE = "protein"
DEFAULT_SUB_MATRIX = "BLOSUM62"
DEFAULT_GAP_OPEN = -5
DEFAULT_GAP_EXT = -1
DEFAULT_TAG_LOG = False


# FUNCTIONS
def get_available_sub_matrices() -> Dict[str, str]:
    """Get the available substitution matrices names from Biopython.

    Returns:
    --------
    sub_matrices_available : dict
        A dictionary containing the available substitution matrices.
        The keys are the sequence type and the values are the names of the matrices available for this sequence type.
    """
    # Get the name of the available substitution matrices
    all_sub_matrixes_available = substitution_matrices.load()
    # Initialize dictionaries for each sequence type
    sub_matrices_available = {
        "protein": [],
        "dna": []
    }
    # Regular expressions to detect matrix types
    dna_regex = re.compile(r"(TRANS|SCHNEIDER|NUC|MEGABLAST|HOXD70)", re.IGNORECASE)
    
    # Iterate through all available matrices and classify them
    for sub_matrix in all_sub_matrixes_available:
        if dna_regex.search(sub_matrix):
            sub_matrices_available["dna"].append(sub_matrix)
        else:
            sub_matrices_available["protein"].append(sub_matrix)
    
    return sub_matrices_available


def get_args() -> Tuple[str, str, str, int, int, str, bool]:
    """Parse the command line arguments.
    
    Returns:
    --------
    
    """
    logger.info("Parsing the command line arguments...")
    # Create the parser
    parser = argparse.ArgumentParser(description="Implementation of a multiple sequence alignment algorithm.")
    # Add the arguments
    parser.add_argument("--f", type=str, required=True, help="The path to the fasta file containing the sequences to align.")
    parser.add_argument("--seq-type", type=str, default=DEFAULT_SEQ_TYPE, help="The type of sequences to align. It can be either 'dna' or 'protein'. The default value is stored in DEFAULT_SEQ_TYPE.")
    parser.add_argument("--sub-matrix", type=str, default=DEFAULT_SUB_MATRIX, help="The substitution matrix to use for the alignment. The default value is stored in DEFAULT_SUB_MATRIX.")
    parser.add_argument("--gap-open", type=int, default=DEFAULT_GAP_OPEN, help="The gap opening penalty. The default value is strored in DEFAULT_GAP_OPEN.")
    parser.add_argument("--gap-ext", type=int, default=DEFAULT_GAP_EXT, help="The gap extension penalty. The default value is strored in DEFAULT_GAP_EXT.")
    parser.add_argument("--job-name", type=str, default=None, help="The name of the output fasta file containing the aligned sequences. The default value is name of the input file with '_aligned' appended.")
    parser.add_argument("--tag-log", type=bool, default=DEFAULT_TAG_LOG, help="Whether to log the pairwise alignment in precision mode or not. The default value is stored in DEFAULT_TAG_LOG.") 
    # Parse the arguments
    args = parser.parse_args()
    # Check the arguments
    # Fasta file does not exist
    if not os.path.isfile(args.f):
        logger.error(f"The file '{args.f}' does not exist.")
        sys.exit(1)
    # Sequence type is not valid
    if args.seq_type not in ["dna", "protein"]:
        logger.error(f"The sequence type '{args.seq_type}' is not valid. It must be either 'dna' or 'protein'.")
        sys.exit(1)
    # Substitution matrix is not available
    if args.sub_matrix not in SUB_MATRIXES_AVAILABLE[args.seq_type]:
        logger.error(f"The substitution matrix '{args.sub_matrix}' is not available. The available substitution matrices for the sequence type '{args.seq_type}' are: {', '.join(SUB_MATRIXES_AVAILABLE[args.seq_type])}.")
        sys.exit(1)
    # Gap penalties are negative
    if args.gap_open > 0:
        logger.error("The gap opening penalty must be a negative integer.")
        sys.exit(1)
    if args.gap_ext > 0:
        logger.error("The gap extension penalty must be a negative integer.")
        sys.exit(1)
    # Job name is not mentioned
    if args.job_name is None:
        # Add _aligned to the input file name
        args.job_name = f"{os.path.splitext(os.path.basename(args.f))[0]}"
    # tag_log is not boolean
    if not isinstance(args.tag_log, bool):
        logger.error("The tag-log flag must be a boolean value.")
    
    # Log the arguments
    logger.debug(f"fasta_file_path: {args.f}")
    logger.debug(f"seq_type: {args.seq_type}")
    logger.debug(f"sub_matrix: {args.sub_matrix}")
    logger.debug(f"gap_open: {args.gap_open}")
    logger.debug(f"gap_ext: {args.gap_ext}")
    logger.debug(f"job_name: {args.job_name}")
    logger.debug(f"tag_log: {args.tag_log}")
    logger.success("Command line arguments parsed successfully.\n")

    return args.f, args.seq_type, args.sub_matrix, args.gap_open, args.gap_ext, args.job_name, args.tag_log


def parse_fasta_to_dict(fasta_file_path: str) -> Dict[str, str]:
    """Parse a fasta file and return the sequences as a dictionary.
    
    Parameters:
    -----------
    fasta_file_path: str
        The path to the fasta file containing the sequences.
    
    Returns:
    --------
    seqs: dict
        A dictionary containing the sequences with the sequence id as key and the sequence as value.
    """
    logger.info("Parsing the fasta file...")
    # Load the sequences
    seqs_complete = SeqIO.to_dict(SeqIO.parse(fasta_file_path, "fasta"))
    # Extract the sequence id and sequence
    seqs = {seq_id: str(seq.seq) for seq_id, seq in seqs_complete.items()}
    logger.debug(f"Number of sequences: {len(seqs)}")
    logger.debug(f"Number of residues max : {max(len(seq) for seq in seqs.values())}")
    logger.success("Fasta file parsed successfully.\n")

    return seqs


def parse_sub_matrix_to_dict(sub_matrix_name: str) -> Dict[str, int]:
    """Parse a substitution matrix and return it as a dictionary.
    
    Parameters:
    -----------
    sub_matrix_name: str
        The name of the substitution matrix to parse.
    
    Returns:
    --------
    sub_matrix_dict: dict[str, int)]
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    """
    logger.info("Parsing the substitution matrix...")
    # Load the substitution matrix
    sub_matrix = substitution_matrices.load(sub_matrix_name)
    logger.debug(f"Substitution matrix:\n {sub_matrix}")
    # Convert the substitution matrix to a dictionary
    sub_matrix_dict = {}
    # Iterate over the rows
    for i, row in enumerate(sub_matrix):
        # Iterate over the columns
        for j, value in enumerate(row):
            # Add the value to the dictionary
            sub_matrix_dict[sub_matrix.alphabet[i] + sub_matrix.alphabet[j]] = value
    logger.success("Substitution matrix parsed successfully.\n")

    return sub_matrix_dict


def backtracking(seq1: str, seq2: str, scores: np.ndarray, sub_matrix: Dict[str, int], gap_open: int, gap_ext: int) -> Tuple[str, str]:
    """ Perform backtracking to find the aligned sequences.

    Parameters:
    -----------
    seq1: str
        The first sequence to align.
    seq2: str
        The second sequence to align.
    scores: np.ndarray
        The score matrix obtained from the Needleman-Wunsch algorithm.
    sub_matrix: dict[str, int]
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    gap_open: int
        The gap opening penalty.
    gap_ext: int
        The gap extension penalty.

    Returns:
    --------
    aligned_seq1: str
        The aligned first sequence.
    aligned_seq2: str
        The aligned second sequence.
    """
    aligned_seq1, aligned_seq2 = [], []
    # Begin at the bottom right corner of the score matrix
    i, j = len(seq1), len(seq2)
    # Backtrack to find the aligned sequences
    while i > 0 and j > 0:
        # Define the current score and the scores for the diagonal, left and top cell
        current_score = scores[i, j]
        diagonal_score = scores[i-1, j-1] 
        left_score = scores[i-1, j]
        top_score = scores[i, j-1]
        # Calculate scores possible from the current cell
        match_score = diagonal_score + sub_matrix[str(seq1[i - 1]) + str(seq2[j - 1])]
        gap_seq1_score = top_score + gap_open
        gap_seq2_score = left_score + gap_open
        # Check which cell the current score came from
        if current_score == match_score:
            # Adding residues in the 2 sequences
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == gap_seq1_score:
            # Adding the residue in seq2 and a gap in seq1
            aligned_seq2.append(seq2[j - 1])
            aligned_seq1.append('-')
            j -= 1
        elif current_score == gap_seq2_score:
            # Adding the residue in seq1 and a gap in seq2
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            raise ValueError(f"The score doesn't equal to one of the possible scores. {current_score} is different from {gap_seq1_score}, {gap_seq2_score}")
            
    # Reverse the aligned sequences (since we backtracked)
    # Join the lists into strings
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2


def pairwise_alignment(first_sequence: tuple, second_sequence: tuple, sub_matrix: Dict[str, int], gap_open: int, gap_ext: int, return_alignment: bool = False, tag_log: bool = False) -> Union[int, list]:
    """Perform a pairwise sequence alignment using the Needleman-Wunsch algorithm.
    
    Parameters:
    -----------
    first_sequence: tuple of str
        The id , sequence of the first sequence to align.
    second_sequence: tuple of str
        The id, sequence of the second sequence to align.
    sub_matrix: dict
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    gap_open: int
        The gap opening penalty. 
    gap_ext: int
        The gap extension penalty.
    return_alignment: bool
        If set to True, the aligned sequences are returned. Default is False.
    tag_log: bool
        If set to True, the log messages are tagged. Good if you want more precision about the alignment. Default is False.
    
    Returns:
    --------
    score_alignment or sequences_aligned : int or List[str, str]
        The alignment score if return_alignment is False, the aligned sequences if return_alignment is True.
    """
    # Get the sequence id and sequence
    id1, seq1 = first_sequence
    id2, seq2 = second_sequence
    if tag_log:
        logger.info(f"Performing a pairwise sequence alignment between {id1} and {id2}...")

    # Get the length of the sequences
    seq1_len, seq2_len = len(seq1), len(seq2)

    # Initialize the score matrix 
    scores = np.zeros((seq1_len + 1, seq2_len + 1), dtype=int)
    # Initialize the first row and column
    scores[:, 0] = np.arange(seq1_len + 1) * gap_open
    scores[0, :] = np.arange(seq2_len + 1) * gap_open
    # Fill the score matrix
    for i in range(1, seq1_len + 1):
        for j in range(1, seq2_len + 1):
            # Calculate the match score
            match = scores[i-1, j-1] + sub_matrix[str(seq1[i - 1]) + str(seq2[j - 1])]
            # Calculate the gap score for seq1 and seq2
            gap_seq1 = scores[i-1, j] + gap_open
            gap_seq2 = scores[i, j-1] + gap_open
            # Calculate the best score
            scores[i, j] = max(match, gap_seq1, gap_seq2)
    
    alignment_score = scores[seq1_len, seq2_len]
    if tag_log:
        logger.debug(f"Alignment score = {alignment_score}")
                 
    if not return_alignment:
        if tag_log:
            logger.success("Pairwise sequence alignment performed successfully.\n")
        return alignment_score
    else:
        # Return the aligned sequences by backtracking through the score matrix
        aligned_seq1, aligned_seq2 = backtracking(seq1, seq2, scores, sub_matrix, gap_open, gap_ext)
        
        if tag_log:
            logger.debug(f"Aligned sequences:")
            logger.debug(f"Sequence {id1}: {aligned_seq1}")
            logger.debug(f"Sequence {id2}: {aligned_seq2}")
            logger.success("Pairwise sequence alignment performed successfully.\n")
        
        return [aligned_seq1, aligned_seq2]


def align_and_update(seq1_id: str, seq2_id: str, seqs: Dict[str, str], sub_matrix: Dict[str, int], gap_open: int, gap_ext: int, tag_log: bool) -> Tuple[str, str, int]:
    """Function to perform pairwise alignment between two sequences.
    
    Parameters:
    -----------
    seq1_id: str
        The id of the first sequence.
    seq2_id: str
        The id of the second sequence.
    seqs: dict
        A dictionary containing the sequences with the sequence id as key and the sequence as value.
    sub_matrix: dict
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    gap_open: int
        The gap opening penalty.
    gap_ext: int
        The gap extension penalty.
    tag_log: bool
        If set to True, the log messages are tagged. Good if you want more precision about the alignment. Default is False.

    Returns:
    --------
    tuple of str, str, int
        The id of the first sequence, the id of the second sequence and the alignment score.
    """
    seq1 = (seq1_id, seqs[seq1_id])
    seq2 = (seq2_id, seqs[seq2_id])
    score = pairwise_alignment(seq1, seq2, sub_matrix, gap_open, gap_ext, False, tag_log)
    return (seq1_id, seq2_id, score)


def perform_parallel_alignment(seqs: Dict[str, str], sub_matrix: Dict[str, int], gap_open: int, gap_ext: int, tag_log: bool = False) -> np.ndarray:
    """Perform a pairwise sequence alignment for all non-redondant pairs of sequences in parallel.
    
    Parameters:
    -----------
    seqs: dict
        A dictionary containing the sequences with the sequence id as key and the sequence as value.
    sub_matrix: dict
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    gap_open: int
        The gap opening penalty. 
    gap_ext: int
        The gap extension penalty.
    tag_log: bool
        If set to True, the log messages are tagged. Good if you want more precision about the alignment. Default is False.
    
    Returns:
    --------
    score_matrix: numpy.ndarray 
        A numpy array contaigning the alignment scores done for all non-redondant pairs of sequences.
    """
    logger.info("Performing a pairwise sequence alignment for all pairs of sequences in parallel...")
    
    # Initialize the score matrix
    num_seqs = len(seqs)
    score_matrix = np.full((num_seqs, num_seqs), np.nan)
    # Create list of sequence pairs to compare (i, j) where i < j to avoid redundancy
    seq_ids = list(seqs.keys())
    seq_pairs = [(seq_ids[i], seq_ids[j]) for i in range(num_seqs) for j in range(i + 1, num_seqs)]
    # Create a loading bar
    with tqdm(total=len(seq_pairs), desc="Aligning sequences", unit="pair") as pbar:
        # Parallel execution
        with ProcessPoolExecutor() as executor:
            # Submit the tasks to the executor
            # Get result of the alignment as key and the pair of sequences as value
            future_to_pair = {executor.submit(align_and_update, seq1_id, seq2_id, seqs, sub_matrix, gap_open, gap_ext, tag_log): (seq1_id, seq2_id) for seq1_id, seq2_id in seq_pairs}
            
            # Process the results as they complete
            for future in as_completed(future_to_pair):
                # Get the pair of sequences and the alignment score
                seq1_id, seq2_id = future_to_pair[future]
                try:
                    seq1_id, seq2_id, score = future.result()
                    i = seq_ids.index(seq1_id)
                    j = seq_ids.index(seq2_id)
                    # Fill the score matrix
                    score_matrix[i, j] = score
                    score_matrix[j, i] = score  # Symmetric assignment (AB = BA)
                except Exception as exc:
                    logger.error(f"Error aligning {seq1_id} and {seq2_id}: {exc}")
                finally:
                    pbar.update()

    # Construct table for the score matrix visualisation
    logger.debug(f"Score matrix: \n {tabulate(score_matrix, showindex=seq_ids, tablefmt='simple_grid')}")
    logger.success("Pairwise sequence alignment for all non-redundant pairs of sequences completed successfully.\n")
    
    return score_matrix


def score_matrix_to_distance_matrix(score_matrix: np.array, seqs: dict) -> np.array:
    """Transform the score matrix into a distance matrix.
    
    Parameters:
    -----------
    score_matrix: np.array
        A numpy array containing the alignment scores done for all non-redondant pairs of sequences.
    seqs: dict
        A dictionary containing the sequences with the sequence id as key and the sequence as value.
    
    Returns:
    --------
    distance_matrix: np.array
        A numpy array containing the distances between all pairs of sequences.
    """
    logger.info("Transforming the score matrix into a distance matrix...")
    seq_ids = list(seqs.keys())
    # Initialize the distance matrix
    distance_matrix = np.zeros_like(score_matrix)
    # Get the maximum score by ignoring the NaN values
    max_score = np.nanmax(score_matrix)

    # Transform each score into a distance
    # distance = 1 - (score / max_score)
    for i in range(distance_matrix.shape[0]):
        for j in range(distance_matrix.shape[1]):
            score = score_matrix[i, j]
            if np.isnan(score):
                distance_matrix[i, j] = np.nan
            else:
                distance_matrix[i, j] = 1 - (score / max_score)
    
    # Construct table for the distance matrix visualisation
    round_dist_matrix = np.round(distance_matrix, 1)
    logger.debug(f"Distance matrix: \n {tabulate(round_dist_matrix, tablefmt='simple_grid', showindex=seq_ids)}")
    logger.success("Transforming the score matrix into a distance matrix completed successfully.\n")
                   
    return distance_matrix


def tuple_to_newick(tree : Tuple[str, Tuple[str, Tuple]]) -> str:
    """Convert a tuple tree into a Newick string.
    
    Parameters:
    -----------
    tree: Tuple[str, Tuple[str, Tuple]]
        A tuple tree containing the sequences and the nested tuples.
    
    Returns:
    --------
    newick: str
        A string in Newick format.
        Example : "(A,(B,(C,D)))"
    """
    if isinstance(tree, tuple):
        return "(" + ",".join(tuple_to_newick(t) for t in tree) + ")"
    else:
        return tree[0]


def flatten_tree(tree_tuple: Tuple) -> List[str]:
    """Flatten a nested tuple tree into a linear list.
    
    Parameters:
    -----------
    tree_tuple: Tuple
        A tuple tree containing the sequences and the nested tuples.
        Example : ("A", ("B", ("C", "D")))
    
    Returns:
    --------
    flat_tree: List[str]
        A list containing the flattened tree.
        Example : ["A", "B", "C", "D"]
    """
    flat_tree = []
    for x in tree_tuple:
        if isinstance(x, tuple):
            flat_tree.extend(flatten_tree(x))  # Recursive call for nested tuples
        else:
            flat_tree.append(x[0])  # Add element if it's not a tuple
    return flat_tree


def perform_upgma(seqs: Dict[str, str], distance_matrix: Dict[str, str]) -> List[str]:
    """Create a guide tree using the UPGMA algorithm.

    It will be used to guide wich sequences to align first.

    Parameters
    -----------
    seqs: Dict[str, str]
        A dictionary containing the sequences with the sequence id as key and the sequence as value.
    distance_matrix: np.array
        A numpy array containing the distances between all pairs of sequences.

    Returns
    --------
    flat_tree: list
        A list containing the flattened guide tree. 
        Example :
            guide_tree = ((A, (B, C)), D)
            flat_tree = [ A, B, C, D]
    """
    logger.info("Constructing the guide tree using UPGMA algorithm...")

    # Initialize the clusters : each sequence is a cluster            
    clusters = {i: [seq_id] for i, seq_id in enumerate(seqs.keys())}
    # Copy the distance matrix to avoid modifying the original
    dist_matrix = np.copy(distance_matrix)

    while len(clusters) > 1:
        # Find the minimum distance in the distance matrix 
        i, j = np.unravel_index(np.nanargmin(dist_matrix), dist_matrix.shape)
        if i == j :
            # Skip the self-similarities
            continue
        else:
            # Fusionnate the clusters i and j
            clusters[i] = (clusters[i], clusters[j])
            del clusters[j]

            # Update the distance matrix
            for k in range(len(dist_matrix)):
                if k != i and not np.isnan(dist_matrix[k, i]):
                    # calculate the new distance between the new cluster and the other clusters
                    avg_dist = (dist_matrix[i, k] + dist_matrix[j, k]) / 2
                    dist_matrix[i, k] = dist_matrix[k, i] = avg_dist
            
            # Set the distances of the merged cluster to NaN
            dist_matrix[:, j] = np.nan
            dist_matrix[j, :] = np.nan
    
    # Get the guide_tree = last cluster
    guide_tree = list(clusters.values())[0]

    # Flatten the guide tree
    flat_tree = flatten_tree(guide_tree)

    # Visualisation
    newick_tree = tuple_to_newick(guide_tree) + ";"
    tree = Tree(newick_tree)
    logger.debug(f"Tree structure : {tree.get_ascii(show_internal=True)}")

    logger.success("Guide tree construction using UPGMA algorithm completed successfully.\n")

    return flat_tree
        

def calculate_sum_of_pairs_score(aligned_seqs: list, new_seq: str, position_i: str, position_j: str, sub_matrix: dict, gap_open: int) -> int:
    """
    Calculate the sum of pairs score for matching a new sequence to an existing alignment.
    
    Parameters
    ----------
    aligned_seqs: List[str]
        List of already aligned sequences.
    new_seq: str
        The new sequence to be aligned.
    position_i: int
        The position in the aligned sequences.
    position_j: int
        The position in the new sequence.
    sub_matrix: Dict[str, int]
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    gap_open: int
        The gap opening penalty.

    Returns
    -------
    score: int
        The sum of pairs score for the given positions in the alignment and new sequence.
    """
    score = 0
    # Compare each aligned sequence with the new sequence at the given positions
    for aligned_seq in aligned_seqs:
        res_seq_aligned = aligned_seq[position_i]
        res_seq_new = new_seq[position_j]
        if res_seq_new == "-" and res_seq_aligned == "-":
            score += 6
        else:
            score += sub_matrix.get((str(res_seq_aligned + res_seq_new)), gap_open)
    
    # Compare residues within the aligned sequences at the current position (avoiding double comparisons)
    for i in range(len(aligned_seqs)):
        seq1 = aligned_seqs[i]
        for j in range(i + 1, len(aligned_seqs)):
            seq2 = aligned_seqs[j]
            res_seq1 = seq1[position_i]
            res_seq2 = seq2[position_i]
            if res_seq1 == "-" and res_seq2 == "-":
                score += 6
            else:
                score += sub_matrix.get((str(res_seq1 + res_seq2)), gap_open)

    return score


def perform_msa(seqs: Dict[str, str], sub_matrix: Dict[str, Dict[str, int]], gap_open: int, gap_ext: int, guide_tree: List[str]) -> Dict[str, str]:
    """ Perform the multiple sequence alignment using the guide tree.
    
    Parameters
    -----------
    seqs: Dict[str, str]
        A dictionary containing the sequences with the sequence id as key and the sequence as value.
    sub_matrix: Dict[str, Dict[str, int]]
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    gap_open: int
        The gap opening penalty.
    gap_ext : int
        The gap extension penalty.
    guide_tree: List[str]
        A list containing the flattened guide tree. 
        Example : [ A, B, C, D]
    
    Returns
    --------
    aligned_seqs: Dict[str, str]
        A dictionary containing the aligned sequences with the sequence id as key and the aligned sequence as value.
    """
    logger.info("Performing the multiple sequence alignment using the guide tree...")
    aligned_seqs = {}

    # Step 1 : Align the first closest sequences from the guide tree
    id1, seq1 = guide_tree[0], seqs[guide_tree[0]]
    id2, seq2 = guide_tree[1], seqs[guide_tree[1]]
    logger.debug(f"Performing pairwise alignment between {id1} and {id2}...")
    aligned_seqs[id1], aligned_seqs[id2] = pairwise_alignment((id1,seq1), (id2,seq2), sub_matrix, gap_open, gap_ext, True)
    current_alignment = [aligned_seqs[id1], aligned_seqs[id2]]
    
    # Step 2 : Iterate over the guide tree to align the other sequences
    for i in range(2, len(guide_tree)):
        seq_id, seq_to_be_aligned = guide_tree[i], seqs[guide_tree[i]]
        logger.debug(f"Adding sequence '{seq_id}' to the existing alignment...")

        # Initialize the score matrix with zeros
        aligned_seqs_len = len(max(current_alignment, key=len))
        seq_to_be_aligned_len = len(seq_to_be_aligned)
        scores = np.zeros((aligned_seqs_len + 1, seq_to_be_aligned_len + 1), dtype=int)
        # Initialize the first row and column
        scores[:, 0] = np.arange(aligned_seqs_len + 1) * gap_open
        scores[0, :] = np.arange(seq_to_be_aligned_len + 1) * gap_open
        
        # Fill the score matrix
        for i in range(1, aligned_seqs_len + 1):
            for j in range(1, seq_to_be_aligned_len + 1):
                # Calculate the match score
                match = scores[i-1, j-1] + calculate_sum_of_pairs_score(current_alignment, seq_to_be_aligned, i - 1, j - 1, sub_matrix, gap_open)
                # Calculate the gap score for seq1 and seq2
                gap_seq1 = scores[i-1, j] + gap_open
                gap_seq2 = scores[i, j-1] + gap_open
                # Calculate the best score
                scores[i, j] = max(match, gap_seq1, gap_seq2)
        logger.debug(f"Score aligment : \n{scores[aligned_seqs_len, seq_to_be_aligned_len]}")

        # Traceback to get the aligned sequence
        align_m = [''] * len(current_alignment)
        align_n = ''
        # Begin at the bottom right corner of the score matrix
        i, j = aligned_seqs_len, seq_to_be_aligned_len
        # Backtrack to find the aligned sequences
        while i > 0 and j > 0:
            # Define the current score and the scores for the diagonal, left and top cell
            current_score = scores[i, j]
            diagonal_score = scores[i-1, j-1] 
            left_score = scores[i-1, j]
            top_score = scores[i, j-1]
            # Calculate scores possible from the current cell
            match_score = diagonal_score + calculate_sum_of_pairs_score(current_alignment, seq_to_be_aligned, i - 1, j - 1, sub_matrix, gap_open)
            gap_seq1_score = top_score + gap_open
            gap_seq2_score = left_score + gap_open
            # Check which cell the current score came from
            if current_score == match_score:
                # Adding residues in the 2 sequences
                for k in range(len(align_m)):
                    seq_tmp = current_alignment[k]
                    align_m[k] += seq_tmp[i - 1]
                align_n += seq_to_be_aligned[j - 1]
                i -= 1
                j -= 1
            elif current_score == gap_seq1_score:
                # Adding the residue in the new sequence to be aligned and a gap in profil1
                align_n += seq_to_be_aligned[j - 1]
                align_m = [seq + '-' for seq in align_m]
                j -= 1
            elif current_score == gap_seq2_score:
                # Adding the residue in profil1 and a gap in seq2
                align_n += '-'
                for k in range(len(align_m)):
                    seq_tmp = current_alignment[k]
                    align_m[k] += seq_tmp[i - 1]
                i -= 1
            else:
                raise ValueError(f"The score doesn't equal to one of the possible scores. {current_score} is different from gap_seq1_score {gap_seq1_score}, gap_seq2_score {gap_seq2_score} and match_score {match_score}")
        
        # Reverse the aligned sequences (since we backtracked)
        # Join the lists into strings
        # aligned_seqs_from_profil = tuple(''.join(reversed(aligned_seq)) for aligned_seq in aligned_profil1)
        #aligned_seq_str = ''.join(reversed(aligned_seq))
        while j > 0:
            align_n += seq_to_be_aligned[j - 1]
            align_m = [seq + '-' for seq in align_m]
            j -= 1
        while i > 0:
            align_n += '-'
            for k in range(len(align_m)):
                seq_tmp = current_alignment[k]
                align_m[k] += seq_tmp[i - 1]
            i -= 1

        # reverse the strings since we traced back from the last position
        align_m = [k[::-1] for k in align_m]
        align_n = align_n[::-1]

        # Update the aligned_seqs dictionary
        aligned_seqs[seq_id] = align_n
        # Update the existing aligned sequences in the dictionary
        for k, key in enumerate(list(aligned_seqs.keys())):
            if key != seq_id:  # Avoid overwriting the new sequence
                aligned_seqs[key] = align_m[k]  # Update the already aligned sequences


        # Update the current alignment
        del current_alignment
        current_alignment = align_m
        current_alignment.append(align_n)

    return aligned_seqs


def save_aligned_seqs(aligned_seqs: Dict[str, str], job_name: str, input_fasta: str, output_format: str = "clustal") -> None:
    """Save the aligned sequences in either FASTA or CLUSTAL format.

    Parameters
    -----------
    aligned_seqs: Dict[str, str]
        A dictionary containing the aligned sequences with the sequence id as key and the aligned sequence as value.
    job_name: str
        The name of the job to save the aligned sequences.
    input_fasta: str
        The path to the input FASTA file for copying headers when using FASTA output.
    output_format: str
        The output format, either 'fasta' or 'clustal'. Defaults to 'clustal'.
    """
    logger.info(f"Saving the aligned sequences in {output_format} format...")
    
    # Verify if the results directory exists
    if not os.path.exists("results"):
        os.mkdir("results")
    
    if output_format == "fasta":
        # Save in FASTA format, preserving the headers from the input FASTA
        output_file = f"results/{job_name}_aligned.fasta"
        with open(input_fasta, "r") as infile, open(output_file, "w") as outfile:
            seq_id = None
            for line in infile:
                if line.startswith(">"):
                    seq_id = line.split()[0][1:].strip()  # Remove '>' and any trailing spaces
                    aligned_seq = aligned_seqs[seq_id]  # Get the aligned sequence
                    outfile.write(f"{line}")
                    # Break aligned sequence into lines of 60 characters for fasta format
                    for i in range(0, len(aligned_seq), 60):
                        outfile.write(aligned_seq[i:i+60] + "\n")
                else:
                    continue  # Skip original sequences in the input file
        logger.success(f"Aligned sequences saved in {output_file}\n")
    
    elif output_format.lower() == "clustal":
        # Save in CLUSTAL format
        output_file = f"results/{job_name}_aligned.clustal"
        with open(output_file, "w") as f:
            f.write("PYCLUSTAL (1.0.0) multiple sequence alignment\n\n")
    
            # TO DO
        
        logger.success(f"Aligned sequences saved in {output_file}\n")


# MAIN PROGRAM
if __name__ == "__main__":
    logger.info("PyClustal is running...\n")
    
    # Load substitution matrix available from Biopython
    SUB_MATRIXES_AVAILABLE = get_available_sub_matrices()

    # Get the command line arguments
    fasta_file_path, seq_type, sub_matrix_name, gap_open, gap_ext, job_name, tag_log = get_args()

    # Load the sequences
    seqs = parse_fasta_to_dict(fasta_file_path)
   
    # Load the substitution matrix
    sub_matrix = parse_sub_matrix_to_dict(sub_matrix_name)

    # Uncomment to test the pairwise alignment
    # seq1, seq2 = ("seq1", "ACGT"), ("seq2", "TACG")
    # aligned_seq = pairwise_alignment(seq1, seq2, sub_matrix, gap_open, gap_ext, return_alignment=True, tag_log=True)

    # Perform a pairwise sequence alignment for all pairs of sequences in parallel
    # Get the score matrix
    score_matrix = perform_parallel_alignment(seqs, sub_matrix, gap_open, gap_ext)

    # Transform the score matrix into a distance matrix
    distance_matrix = score_matrix_to_distance_matrix(score_matrix, seqs)

    # Create the phylogenetic tree
    tree = perform_upgma(seqs, distance_matrix)

    # Perform multiple sequence alignment
    aligned_seqs = perform_msa(seqs, sub_matrix, gap_open, gap_ext, tree)

    # Save the aligned sequences
    save_aligned_seqs(aligned_seqs, job_name, fasta_file_path, output_format="fasta")

    logger.info("PyClustal has finished running.\n")



    

