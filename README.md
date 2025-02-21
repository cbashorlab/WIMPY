# WIMPY: A software package for nanopore sequencing analysis of combinatorial genetic libraries of arbitrary length scales

<!-- TODO: add link to preprint/journal -->
Checkout the full article at [OUP Bioinformatics]()!

![abstract](./pipeline.jpg)

`wimpy` (What’s In My Pot, Y’all) is a software package with implementations in both python and MATLAB, that can analyze large-scale combinatorial libraries of synthetic DNA. WIMPY leverages features that are unique to libraries of synthetic DNA, such as known locations of expected diversity flanked by constant sequences, and uses localized containment search algorithms for rapid and accurate variant assignment on a single-read basis, without relying on any consensus-based sequence aggregation methods.

## Getting Started - Python

We recommend running `wimpy` with UNIX-based operating system (Linux/MacOS). For Windows users we recommend to use [Windows Subsystem for Linux (WSL 2)](https://learn.microsoft.com/en-us/windows/wsl/install) for better compatibility (although it should also be compatible with Windows installation of Python and Anaconda, but it's not officially supported yet).

### Download Conda and Clone Repository

Install [Miniconda](https://docs.anaconda.com/miniconda/install/#quick-command-line-install) (recommended) or [Anaconda](https://docs.anaconda.com/anaconda/install/). For Miniconda you can install it with the following commands:

- Linux/WSL:

    ```bash
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    ```

- MacOS (M1 or later):

    ```bash
    mkdir -p ~/miniconda3
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    ```

Clone the repository into your local directory:

```bash
git clone git@github.com:cbashorlab/WIMPY.git
```

### Setting Up Environment and Package

Open folder `wimpy_python`:

```bash
cd ./WIMPY/wimpy_python
```

create `wimpy` virtual environment from the `environment.yml` file:

```bash
conda env create -f environment.yml
```

activate the environment:

```bash
conda activate wimpy
```

To use `wimpy` as a package, install it with the following command:

```bash
pip install -e .
```

### Using WIMPY

Checkout [`example_script_python.ipynb`](./wimpy_python/example_script_python.ipynb) for an exmple of how to use `wimpy` to process sequencing files in a pipeline.

## Getting Started - MATLAB

- Install the latest version of [MATLAB](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html)
- Clone the repository into your local directory:

    ```bash
    git clone git@github.com:cbashorlab/WIMPY.git
    ```

- Open folder `wimpy_matlab` in MATLAB
- Add `wimpy_helper_functions` to path by running the following command on your MATLAB console:

    ```MATLAB
    addpath('./wimpy_helper_functions')
    ```

  - Alternatively, you can add `wimpy_helper_functions` to path by going to "Current Folder" -> rightclick `wimpy_helper_functions` -> "Add to Path" -> "Selected Folder"

- Checkout [`example_script_matlab.m`](./wimpy_matlab/example_script_matlab.m) for an exmple of how to use `wimpy` to process sequencing files in a pipeline.

## WIMPY Functions Documentation

***⚠️ Work in Progress - making edits to the functions at the moment, documentation may be inconsistent. Documentation will be updated once edit on functions are finalized***

### `rev_comp(seq: str)`

Returns the reverse complement of input sequence.

### `fastqall(directory: PathLike = "./", prefix: str = "", idx_start: int = None, idx_end: int = None)`

Reads all fastq files in the given directory with the given prefix. Returns all the quality scores, lengths, and sequences as a list.

**Args:**

- `directory` (str, optional): The location of files. Defaults to './'.
- `prefix` (str, optional): The specified prefix. Defaults to ''.
- `idx_start`: Start index of the file names. Defaults to None.
- `idx_end`: End index of the file names. Defaults to None.

**Returns:**

- `q_scores`: All the quality scores.
- `lengths`: Lengths of all reads.
- `seqs`: Sequence of all the reads.

### `to_tiles(seq: str, tile_len: int = 10)`

Breaks the sequence into tiles, given the tile length.

**Args:**

- `seq` (str): The sequence to break.
- `tile_len` (int, optional): Length of tile. Defaults to 10.

**Returns:**

- `np.array(str)`: A tuple of tiles.

### `bowtile(seqs, ref, thresh=0.03, tile_len=10, max_len=100)`

Uses tiling to determine the occurrence of a sequence in a list of reads.

**Args:**

- `seqs` (list): A list of reads to be tiled to.
- `ref` (str): The reference sequence.
- `thresh` (float): Threshold to determine whether the tiling is valid.
- `tile_len` (int, optional): Length of tile. Defaults to 10.
- `max_len` (int, optional): Maximum length of ref seq. Defaults to 100.

**Returns:**

- `new_seq` (list): Valid sequences that are aligned to start with ref seq.
- `right_seq` (list): Valid sequences.
- `flip` (list): Whether the match is on fwd strand (0), rev strand (1), or no match (-1).

### `tilepin(seqs, ref, thresh=0.03, tile_len=10, verbose=False)`

Tiles the sequences and finds matches.

**Args:**

- `seqs` (list): Input sequences.
- `ref` (str): Reference sequence.
- `thresh` (float, optional): Threshold for valid tiling. Defaults to 0.03.
- `tile_len` (int, optional): Length of tile. Defaults to 10.
- `verbose` (bool, optional): Print progress. Defaults to False.

**Returns:**

- `found`: Number of tiles found.
- `positions`: Positions of tiles.
- `indices`: Indices of tiles.

### `tilepin_v2(seqs, ref, thresh=0.03, tile_len=10, verbose=False)`

Improved tilepin using hashmap search to determine whether a tile is in certain region of sequences.

**Args:**

- `seqs` (list[str]): Input sequences to search.
- `ref` (str): The reference sequence.
- `thresh` (float, optional): Threshold to determine if the sequence contains the reference sequence. Defaults to 0.03.
- `tile_len` (int, optional): Length of tile. Defaults to 10.
- `verbose` (bool, optional): Print out progress. Defaults to True.

**Returns:**

- `num_matches`: The number of tiles matched the sequence.
- `match_index`: The index where the reference sequence is located.
- `matches`: The list of all indices where match is found.

### `chophat(seqs, positions, end_positions=None, max_length=None, retain=True)`

Chops the sequence to keep only the regions of interest.

**Args:**

- `seqs` (list[str]): Sequences to be processed.
- `positions` (np.array[int]): Start position for truncation.
- `end_positions` (np.array[int], optional): End position for truncation. Defaults to None.
- `max_length` (int, optional): The max length of truncated sequence. Defaults to None.
- `retain` (bool, optional): Retain the sequence if the start position to the end of string is smaller than max_length. Defaults to True.

**Returns:**

- `list[str]`: The truncated sequences.

### `viscount(seqs, ref_seqs, thresh, return_confusion_matrix=True, tile_len=10, verbose=False)`

Counts the number of occurrences of tiles in a list of reference sequences. Constructs a confusion matrix of index assignments if requested.

**Args:**

- `seqs` (list[str]): Sequences to be processed.
- `ref_seqs` (list[str]): List of reference sequences.
- `thresh` (int): Threshold for a valid occurrence in confusion matrix.
- `return_confusion_matrix` (bool, optional): Construct the confusion matrix or not. Defaults to True.
- `tile_len` (int, optional): Length of tiles. Defaults to 10.
- `verbose` (bool, optional): Print out progress. Defaults to True.

**Returns:**

- `match_ratios`: The proportion of reference tiles that get assigned to the sequences.
- `match_counts`: The raw count of number of tiles get assigned to the sequences.
- `conf_matrix` (optional): Confusion matrix for the assignment.

### `FASTar(pregions, ref, step, bw)`

Python equivalent of the FASTar function in MATLAB.

**Args:**

- `pregions` (list[str]): List of pre-regions to search within.
- `ref` (str): Reference sequence.
- `step` (int): Step size for sliding window.
- `bw` (float): Bandwidth for kernel density estimation.

**Returns:**

- `nums` (np.array): Number of local maxima found in each pre-region.
- `locs` (list): List of locations of local maxima for each pre-region.
