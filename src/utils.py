"""This module contains utility functions that are used in the project."""

# METADATA
__authors__ = "Essmay Touami"
__contact__ = "essmay.touami@etu.u-paris.fr"
__date__ = "september 2024"
__version__ = "1.0.0"


# LIBRARY IMPORTS
import os
import sys
import numpy as np
from collections import defaultdict
from typing import Tuple, Dict, List, Union


import re
from loguru import logger
from Bio.Align import substitution_matrices

# LOCAL IMPORTS


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


def validate_args(args):
    """
    Validate the command line arguments.

    Parameters
    ----------
    args : Namespace
        The command line arguments.

    Raises
    ------
    SystemExit
        If any of the validation checks fail, exits the program with an error code.
    """
    # Load substitution matrix available from Biopython
    SUB_MATRIXES_AVAILABLE = get_available_sub_matrices()

    # Fasta file does not exist
    if not os.path.isfile(args.f):
        logger.error(f"The file '{args.f}' does not exist.")
        sys.exit(1)

    # Sequence type is not valid
    if args.seq_type not in ["dna", "protein"]:
        logger.error(f"The sequence type '{args.seq_type}' is not valid. It must be either 'dna' or 'protein'.")
        sys.exit(1)

    # Substitution matrix is not available
    if args.sub_matrix not in SUB_MATRIXES_AVAILABLE.get(args.seq_type, []):
        logger.error(f"The substitution matrix '{args.sub_matrix}' is not available. The available substitution matrices for the sequence type '{args.seq_type}' are: {', '.join(sub_matrices_available.get(args.seq_type, []))}.")
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
        sys.exit(1)


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
    """Wrapper function to perform pairwise alignment between two sequences for parallel processing.
    
    This function calls the `pairwise_alignment` function to compute the alignment score between two sequences.
    It is specifically designed to provide a simplified interface for parallel processing of sequence alignments.
    
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


def get_conservation_groups_from_matrix(sub_matrix: Dict[str, int], strong_threshold: float = 0.5, weak_threshold: float = 0) -> Tuple[List[set], List[set]]:
    """Extract strong and weak conservation groups from a substitution matrix.
    
    Parameters
    ----------
    sub_matrix : Dict[str, int]
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.
        Example : {"AA": 1, "AC": -1, "AD": -2, "AE": -1, ...}
    strong_threshold : float
        The score threshold for strong conservation.
    weak_threshold : float
        The score threshold for weak conservation (between 0 and strong_threshold).
    
    Returns
    -------
    strong_similarity_groups : List[set]
        List of sets of strongly conserved amino acids.
    weak_similarity_groups : List[set]
        List of sets of weakly conserved amino acids.
    """
    strong_similarity_groups = defaultdict(set)
    weak_similarity_groups = defaultdict(set)

    for pair, score in sub_matrix.items():
        aa1, aa2 = pair[0], pair[1]

        if aa1 != aa2:

            if score > strong_threshold:
                strong_similarity_groups[aa1].add(aa2)
                strong_similarity_groups[aa2].add(aa1)
            
            elif weak_threshold < score <= strong_threshold:
                weak_similarity_groups[aa1].add(aa2)
                weak_similarity_groups[aa2].add(aa1)

    strong_similarity_list = [group for group in strong_similarity_groups.values() if len(group) > 1]
    weak_similarity_list = [group for group in weak_similarity_groups.values() if len(group) > 1]

    return strong_similarity_list, weak_similarity_list


def get_consensus_symbol(column_residues: List[str], sub_matrix: Dict[str, int]) -> str:
    """Determine the consensus symbol (*, :, .) for a column of aligned residues.

    Parameters
    -----------
    column_residues: List[str]
        The list of residues in the column.
    sub_matrix: Dict[str, int]
        The substitution matrix as a dictionary with the amino acid pair as key and the substitution score as value.

    Returns
    --------
    symbol: str
        The consensus symbol (*, :, .) based on the conservation of the column.
    """
    # Get the conservation groups from the substitution matrix
    strong_groups, weak_groups = get_conservation_groups_from_matrix(sub_matrix)

    # Filter out gaps ('-') from the column
    column_residues = [res for res in column_residues if res != '-']
    
    if len(set(column_residues)) == 1:
        return '*'
    
    for group in strong_groups:
        if all(res in group for res in column_residues):
            return ':'
    
    for group in weak_groups:
        if all(res in group for res in column_residues):
            return '.'
    
    return ' '


