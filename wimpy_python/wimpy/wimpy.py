from Bio import SeqIO, Align
from Bio.Seq import Seq
from os import path
from tqdm.auto import tqdm
from glob import glob
import re
import numpy as np
import warnings
from os import PathLike
from scipy.stats import gaussian_kde
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity


COMPLEMENT = str.maketrans("ATGC", "TACG")


def rev_comp(seq: str):
    """Returns the reverse complement of input sequence."""
    return seq[::-1].translate(COMPLEMENT)


def fastqall(
    directory: PathLike = "./",
    prefix: str = "",
    idx_start: int = None,
    idx_end: int = None,
    return_names: bool = False,
):
    """read all fastq files in the given directory with given prefix. returns
    all the quality scores, lengths, and sequences as a list

    Args:
        file_directory (str, optional): the location of files Defaults to './'
        prefix (str, optional): the specified prefix. Defaults to ''.
        idx_start: start index of the file names. Defaults to None.
        idx_end: end index of the file names. Defaults to None.

    Returns:
        q_scores: all the quality scores
        lengths: lengths of all reads
        seqs: sequence of all the reads

    Example:
    ```
    >>> import numpy as np
    >>> from wimpy import fastqall
    >>> q_scores, lengths, seqs = fastqall(directory="./data", return_names=False)
    >>> print("Number of reads:", len(seqs))
    >>> print("First sequence:", seqs[0])
    ```
    """

    # setup directory
    if not path.isdir(directory):
        print("directory not found, using current directory")
        directory = "./"

    q_scores = []
    lengths = []
    seqs = []
    names = []

    # read all files and store q_scores, lens, and seqs
    file_names = sorted(glob(path.join(directory, f"{prefix}*.fastq")))
    for file in tqdm(file_names[idx_start:idx_end], desc="reading fastq files"):
        for read in SeqIO.parse(file, "fastq"):
            q_scores.append(read.letter_annotations["phred_quality"])
            lengths.append(len(read))
            seqs.append(str(read.seq))
            names.append(read.id)

    if return_names:
        return q_scores, lengths, seqs, names
    return q_scores, lengths, seqs


def to_tiles(seq: str, tile_len: int = 10):
    """break the sequence into tiles, given the tile length

    Args:
        seq (str): the sequence to break
        tile_len (int, optional): Length of tile. Defaults to 10.

    Returns:
        np.array(str) : an array of tiles
    """

    return np.array([seq[i : i + tile_len] for i in range(len(seq) - tile_len + 1)])


def bowtile(seqs, ref, thresh=0.03, tile_len=10, max_len=100, verbose=False):
    """
    use containment search of the fragmented tiles of the reference sequence
    to determine the occurence of a sequence in a list of reads

    Args:
        seqs (list): a list of reads to be tiled to
        ref (str): the reference sequence
        thresh (float): threshold to determine whether the tiling is valid
        tile_len (int, optional): lrngth of tile. Defaults to 10.
        max_len (int, optional): maximum length of ref seq. Defaults to 100.
        verbose (bool, optional): print out progress. Defaults to False.

    Returns:
        new_seq (list): valid sequences that are aligned to start with ref seq
        right_seq (list): valid sequences
        flip (list): whether the match is on fwd strand (0), rev
        strand (1), or no match (-1)

    Example:
    ```
    >>> seqs = [
    ...     "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
    ...     "GCTAGCTAGCTAGCTAGCTAGCTACGTACGTACGTAC",
    ...     "TGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    ... ]
    >>> ref = "ATGCGTACGTAGCTAGCTAGC"

    >>> # Run bowtile
    >>> new_seq, right_seq, flip = wp.bowtile(seqs, ref, thresh=0.03, tile_len=5, max_len=20)
    >>> print("New Sequences:", new_seq)
    >>> print("Original Sequences:", right_seq)
    >>> print("Strand Orientation (0=fwd, 1=rev, -1=fail):", flip)
    ```
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
    if verbose:
        seqs = tqdm(seqs, desc="bowtile progress")
    for seq in seqs:
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
    """Deprecated, use `tilepin_v2` instead"""
    warnings.warn(
        "tilepin will be deprecated soon, please use tilepin_v2 instead. tilepin_v2 is faster and more efficient.",
        DeprecationWarning,
        stacklevel=2,
    )

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
    a reference sequence is present within reads in a given array.

    Args:
        seqs (list[str]): input sequences to search
        ref (str): the reference sequence
        thresh (float, optional): threshold to determine if the
        sequence contains the reference sequence. Defaults to 0.03.
        tile_len (int, optional): length of tile. Defaults to 10.
        verbose (bool, optional): print out progress. Defaults to False.

    Returns:
        num_matches: the number of tiles matched the sequence
        match_index: the index where the reference sequence is located.
        matches: the list of all indices where match is found

    Example:
    ```
    >>> seqs = [
    ...     "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
    ...     "CTAGCCATACTCGGACATGCGTACGTATCTAGCTAGC",
    ...     "TGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    ... ]
    >>> ref = "ATGCGTACGTAGCTAGCTAGC"
    >>> num_matches, match_index, matches = wp.tilepin_v2(seqs, ref, thresh=0.03, tile_len=5)
    >>> print("Number of matches:", num_matches)
    >>> print("Match index:", match_index)
    >>> print("Matches array:", matches)
    ```
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

    Example:
    ```
    >>> import numpy as np
    >>> seqs = [
    ...     "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
    ...     "GCTAGCTAGCTAGCTAGCTAGCTACGTACGTACGTAC",
    ...     "TGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    ... ]
    >>> positions = np.array([0, 5, 10])
    >>> end_positions = np.array([10, 20, 30])
    >>> truncated_seqs = wp.chophat(seqs, positions, end_positions, max_length=15)
    >>> print(truncated_seqs)
    ```
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

    Example:
    ```
    >>> seqs = [
    ...     "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
    ...     "GCTAGCTAGCTAGCTAGCTAGCTACGTACGTACGTAC",
    ...     "TGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    ... ]
    >>> ref_seqs = ["ATGCGTACGTAGCTAGCTAGC", "CATGCATGCATGCATGCA"]
    >>> match_ratios, match_counts, conf_matrix = wp.viscount(seqs, ref_seqs, thresh=0.03, tile_len=5)
    >>> print("Match Ratios:\n", match_ratios)
    >>> print("Match Counts:\n", match_counts)
    >>> print("Confusion Matrix:\n", conf_matrix)
    ```
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


def fastar(seqs, ref, tile_len=6, bw=None):
    #! not final, to be updated after matlab script update
    """
    Uses kernel density estimation to determine the occurence, or 
    multiple occurences of a motif in the sequence. 

    Args:
        seqs (list[str]): list of sequence to search within.
        ref (str): reference sequence.
        tile_len (int): size of tile for sliding window.
        bw (float): bandwidth for kernel density estimation.

    Returns:
        nums (np.array): number of local maxima found in each sequence.
        locs (list): list of locations of local maxima for each sequence.

    Example:
    ```
    >>> seqs = [
    ...     "ACTGATCGACTGATCGACTGATGGACTGATCGACTG",
    ...     "GCTAGCTAGCTAGACTGATCGAGCTACGTAACGTAC",
    ...     "TGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    ... ]
    >>> ref = "ACTGAT"
    >>> nums, locs = wp.fastar(seqs, ref, tile_len=3)
    >>> print("Number of local maxima:", nums)
    >>> print("Locations of local maxima:", locs)
    ```
    """

    if bw is None:
        bw = 0.5 * len(ref)

    nums = np.zeros(len(seqs))
    locs = [None] * len(seqs)

    for j, x in enumerate(seqs):
        a1 = []
        if x == "":
            continue

        tiles = to_tiles(ref, tile_len=tile_len)
        for tile in tiles:
            a1 += [m.start() for m in re.finditer(tile, x)]

        if a1:
            # kde = gaussian_kde(a1, bw_method=0.005*bw)
            # density = kde.evaluate(np.linspace(min(a1)-100, max(a1)+100, max(a1)+100))
            kde_skl = KernelDensity(bandwidth=bw)
            kde_skl.fit(np.array(a1).reshape(-1, 1))
            density = np.exp(
                kde_skl.score_samples(
                    np.linspace(min(a1) - 100, max(a1) + 100, max(a1) + 100).reshape(
                        -1, 1
                    )
                )
            )
            dmax = argrelextrema(density, np.greater)
            dval = density[dmax[0]]
            dmax = dmax[0][dval > 0.1 * max(dval)]
            locs[j] = dmax
            nums[j] = len(dmax)

    return nums, locs


def barcoat(seqs, preset="BBA", barcode_construct=None, alinger=None):
    """
    Aligns sequences to a barcode construct using a preset or custom aligner.
    Args:
        seqs (list of str): List of sequences to be aligned.
        preset (str, optional): Preset configuration for barcode construct and aligner.
        Options are "BBA", "DDC", or None. Default is "BBA".
        barcode_construct (Bio.Seq.Seq, optional): Custom barcode construct sequence.
        Required if preset is None.
        aligner (Bio.Align.PairwiseAligner, optional): Custom aligner object.
        Required if preset is None.

    Returns:
        barcode (numpy.ndarray): Array of aligned barcode sequences.
        position (numpy.ndarray): Array of starting positions of alignments.
        length (numpy.ndarray): Array of lengths of alignments.
        score (numpy.ndarray): Array of alignment scores.

    Raises:
    ValueError: If preset is None and either barcode_construct or aligner is not provided.

    Example:
    ```
    >>> from Bio.Seq import Seq
    >>> from Bio import Align
    >>> import numpy as np
    >>> seqs = ["ATTATTATTATTATTATTA", "CTTCTTCTTCTTCTTCTTC"]
    >>> barcode_construct = Seq("ATTATTATTATTATTATTA")
    >>> aligner = Align.PairwiseAligner()
    >>> aligner.mode = "local"
    >>> aligner.substitution_matrix = Align.substitution_matrices.Array(
    ...     data=np.array(
    ...         [[5, -5, -5, -5], [-5, 5, 5, 5], [-5, 5, 5, 5], [-5, 5, 5, 5]]
    ...     ),
    ...     alphabet="ATGC",
    ... )
    >>> aligner.open_gap_score = -8
    >>> aligner.extend_gap_score = -8
    >>> barcoat(seqs, preset=None, barcode_construct=barcode_construct, aligner=aligner)
    ```
    """

    if preset not in ["BBA", "DDC", None]:
        raise ValueError(
            f"preset '{preset}' is not recognized. must be either 'BBA', 'DDC', or None."
        )

    if preset is None and (barcode_construct is None or aligner is None):
        raise ValueError(
            "barcode_construct and aligner must be provided if preset is None"
        )

    # load barcode structure and substitution matrix
    if preset == "BBA":
        barcode_construct = Seq("ATTATTATTATTATTATTA")

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = Align.substitution_matrices.Array(
            data=np.array(
                [[5, -5, -5, -5], [-5, 5, 5, 5], [-5, 5, 5, 5], [-5, 5, 5, 5]]
            ),
            alphabet="ATGC",
        )
        aligner.open_gap_score = -8
        aligner.extend_gap_score = -8

    elif preset == "DDC":
        barcode_construct = Seq("CTTCTTCTTCTTCTTCTTC")

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = Align.substitution_matrices.Array(
            data=np.array(
                [[5, 5, 5, -5], [5, 5, 5, -5], [5, 5, 5, -5], [-5, -5, -5, 5]]
            ),
            alphabet="ATGC",
        )
        aligner.open_gap_score = -8
        aligner.extend_gap_score = -8

    # else: use the provided barcode_construct and aligner

    num_seqs = len(seqs)
    barcode = np.empty(num_seqs, dtype="<U32")
    position = np.zeros(num_seqs, dtype=np.int32) - 1
    length = np.zeros(num_seqs, dtype=np.int32) - 1
    score = np.zeros(num_seqs, dtype=np.int32) - 1

    for i, seq in enumerate(tqdm(seqs, desc="searching barcode using alignment")):
        if len(seq) <= 0:
            continue

        alignment = aligner.align(seq, barcode_construct)[0]
        barcode[i] = alignment[0]
        position[i] = int(alignment.indices[0, 0])
        length[i] = alignment.length
        score[i] = alignment.score

    return barcode, position, length, score
