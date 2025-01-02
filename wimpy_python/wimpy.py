from Bio import SeqIO
from Bio.Seq import Seq
from os import path
from tqdm import tqdm
import glob
import re
import numpy as np
from typing import List, Tuple
from os import PathLike


COMPLEMENT = str.maketrans("ATGC", "TACG")


def rev_comp(seq: str):
    return seq[::-1].translate(COMPLEMENT)


def fastqall(directory: PathLike = "./", prefix: str = ""):
    """read all fastq files in the given directory with given prefix. returns
    all the quality scores, lengths, and sequences as a list

    Args:
        file_directory (str, optional): the location of files Defaults to './'
        prefix (str, optional): the specified prefix. Defaults to ''.

    Returns:
        q_scores: all the quality scores
        lengths: lengths of all reads
        seqs: sequence of all the reads
    """

    # setup directory
    if not path.isdir(directory):
        print("directory not found, using current directory")
        directory = "./"

    q_scores = []
    lengths = []
    seqs = []

    # read all files and store q_scores, lens, and seqs
    for file in tqdm(
        glob.glob(path.join(directory, f"{prefix}*.fastq")), desc="reading fastq files"
    ):
        for read in SeqIO.parse(file, "fastq"):
            q_scores.append(read.letter_annotations["phred_quality"])
            lengths.append(len(read))
            seqs.append(str(read.seq))

    return q_scores, lengths, seqs


def to_tiles(seq: str, tile_len: int = 10):
    """break the sequence into tiles, given the tile length

    Args:
        seq (str): the sequence to break
        tile_len (int, optional): Length of tile. Defaults to 10.

    Returns:
        np.array(str) : a tuple of tiles
    """

    return np.array([seq[i : i + tile_len] for i in range(len(seq) - tile_len + 1)])


def bowtile(seqs, ref, thresh=0.03, tile_len=10, max_len=100):
    """use tiling to determine the occurence of a sequence in a list of reads

    Args:
        seqs (list): a list of reads to be tiled to
        ref (str): the reference sequence
        thresh (float): threshold to determine whether the tiling is valid
        tile_len (int, optional): lrngth of tile. Defaults to 10.
        max_len (int, optional): maximum length of ref seq. Defaults to 100.

    Returns:
        new_seq (list): valid sequences that are aligned to start with ref seq
        right_seq (list): valid sequences
        flip (list): whether the match is on fwd strand (0), rev strand (1), or no match
    """

    # convert to upper case
    seqs = [s.upper() for s in seqs]
    ref = ref.upper()

    # calculate forward and backward reference sequence
    ref_f = ref[:max_len]
    ref_r = rev_comp(ref_f)

    # pre-generate all the tiles
    tiles_f = to_tiles(ref_f, tile_len)
    tiles_r = to_tiles(ref_r, tile_len)

    flip = []
    new_seq = []
    right_seq = []

    # tiling through all sequences
    for seq in tqdm(seqs, desc="bowtile progress"):
        idx_f = []
        idx_r = []

        # find matchs of all tiles
        for tile_f, tile_r in zip(tiles_f, tiles_r):
            idx_f += [x.start() for x in re.finditer(tile_f, seq)]
            idx_r += [x.start() for x in re.finditer(tile_r, seq)]

        # calculate pesudo-median for all indices where match is found
        pos_f = sorted(idx_f)[len(idx_f) // 2] if idx_f else 0
        pos_r = sorted(idx_r)[len(idx_r) // 2] if idx_r else 0

        diff_num_match = (len(idx_f) - len(idx_r)) / (len(ref_f) - tile_len)

        # foward strand
        if diff_num_match > thresh:
            new_seq.append(seq[pos_f:] + seq[:pos_f])
            right_seq.append(seq)
            flip.append(0)
        # reverse strand
        elif diff_num_match < -thresh:
            new_seq.append(rev_comp(seq[pos_r:] + seq[:pos_r]))
            right_seq.append(rev_comp(seq))
            flip.append(1)
        # tiling failed
        else:
            new_seq.append("")
            right_seq.append(seq)
            flip.append(-1)

    return new_seq, right_seq, flip


def tilepin(seqs, ref, thresh=0.03, tile_len=10, verbose=False):
    seqs = [s.upper() for s in seqs]
    ref = ref.upper()

    # pre-generate all the tiles
    tiles = np.array(to_tiles(ref, tile_len))

    found = np.zeros(len(seqs))
    positions = np.zeros(len(seqs)) - 1
    indices = [[None for _ in tiles] for _ in seqs]

    # tiling through all sequences
    if verbose:
        seqs = tqdm(seqs, desc="tilepinning progress")
    for i, seq in enumerate(seqs):

        # find matchs of all tiles
        for j, tile in enumerate(tiles):
            indices[i][j] = [x.start() for x in re.finditer(tile, seq)]

        # calculate pesudo-median for all indices where match is found
        all_matches = np.concatenate(indices[i])
        found[i] = len(all_matches)

        if len(all_matches) / len(ref) > thresh:
            positions[i] = sorted(all_matches)[len(all_matches) // 2]

    return found, positions, indices


def tilepin_v2(seqs, ref, thresh=0.03, tile_len=10, verbose=False):
    """
    Improved tilepin using hashmap search to determine whether
    a tile is in certain region of sequences.

    Args:
        seqs (list[str]): input sequences to search
        ref (str): the reference sequence
        thresh (float, optional): threshold to determine if the
        sequence contains the reference sequence. Defaults to 0.03.
        tile_len (int, optional): length of tile. Defaults to 10.
        verbose (bool, optional): print out progress. Defaults to True.

    Returns:
        num_matches: the number of tiles matched the sequence
        match_index: the index where the reference sequence is located.
        matches: the list of all indices where match is found
    """
    seqs = [s.upper() for s in seqs]
    ref = ref.upper()

    # breaking sequences into tiles and store in hashmap (dictionary)
    ref_tiles = set(to_tiles(ref, tile_len=tile_len))
    matches = np.zeros((len(seqs), max(len(seq) for seq in seqs)), dtype=np.int32)
    if verbose:
        seqs = tqdm(seqs, desc="match sequences to reference")
    for i, seq in enumerate(seqs):
        for j in range(len(seq) - tile_len):
            if seq[j : j + tile_len] in ref_tiles:
                matches[i, j] = 1

    num_matches = np.sum(matches, axis=1)
    match_index = np.zeros_like(num_matches) - 1
    
    is_above_thresh = (num_matches / len(ref)) > thresh
    match_index[is_above_thresh] = np.array(
        [np.median(np.nonzero(match)) for match in matches[is_above_thresh]]
    )

    return num_matches, match_index, matches


def chophat(seqs, positions, end_positions=None, max_length=None, retain=True):
    """
    chop the sequence to keep only the regions of interest.

    Args:
        seqs (list[str]): sequences to be processed.
        positions (np.array[int]): start position for truncation
        end_positions (np.array[int], optional): end position for
        truncation. Defaults to None.
        max_length (int, optional): the max length of truncated sequence. If
        max_length is None, the truncated sequence will be from the start position
        to the end of string. Defaults to None.
        retain (bool, optional): retain the sequence if the start position to the
        end of string is smaller than max_length. Defaults to True.

    Returns:
        list[str]: the truncated sequences
    """
    truncated_sequences = []

    if end_positions is None:
        assert len(seqs) == len(
            positions
        ), "sequences and positions must be the same size"

        for seq, pos in zip(seqs, positions):
            if pos <= 0:
                truncated_sequences.append("")
                continue

            if max_length is None:
                truncated_sequences.append(seq[pos[0] :])
                continue

            if not retain and len(seq) - pos > max_length:
                truncated_sequences.append("")
                continue

            truncated_sequences.append(seq[pos : pos + max_length])

    else:
        assert (
            len(seqs) == len(positions) == len(end_positions)
        ), "sequences, start and end positions must be the same size"

        for seq, start, end in zip(seqs, positions, end_positions):
            if start >= 0 and end > start:
                truncated_sequences.append(seq[start:end])
            else:
                truncated_sequences.append("")

    return truncated_sequences


def viscount(
    seqs, ref_seqs, thresh, return_confusion_matrix=True, tile_len=10, verbose=False
):
    """
    Count the number of occurence of tiles in a list of reference sequences.
    Construct a confusion matrix of index assignments if requested.

    Args:
        seqs (list[str]): sequences to be processed.
        ref_seqs (list[str]): list of reference sequences.
        thresh (int): threshold for a valid occurence in confusion matrix
        return_confusion_matrix (bool, optional): construct the confusion matrix
        or not. Defaults to True.
        tile_len (int, optional): length of tiles. Defaults to 10.
        verbose (bool, optional): printing out progress. Defaults to True.

    Returns:
        match_ratios: the proportion of reference tiles that get assigned to
        the sequences.
        match_counts: the raw count of number of tiles get assigned to the
        sequences.
        conf_matrix (optional): confusion matrix for the assignment
    """
    ref_lengths = np.array([len(ref) for ref in ref_seqs])
    match_counts = np.zeros((len(seqs), len(ref_seqs)))

    if verbose:
        ref_seqs = tqdm(ref_seqs, desc="matching to reference sequences")
    for j, ref_seq in enumerate(ref_seqs):
        tiles = to_tiles(ref_seq.upper(), tile_len=tile_len)
        for tile in tiles:
            for i, seq in enumerate(seqs):
                if tile in seq:
                    match_counts[i, j] += 1

    match_ratios = match_counts / ref_lengths
    if not return_confusion_matrix:
        return match_ratios, match_counts

    conf_matrix = np.zeros((len(ref_seqs), len(ref_seqs)))
    for i in range(len(ref_seqs)):
        for j in range(len(ref_seqs)):
            conf_matrix[i, j] = np.sum(
                (match_ratios[:, i] > thresh) & (match_ratios[:, j] > thresh)
            )

    return match_ratios, match_counts, conf_matrix
