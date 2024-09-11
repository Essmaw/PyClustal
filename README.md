<h1 align="center">
  <img style="vertical-align:middle; width:70%; position:fixed;"
  src="/data/img/banner.png">
</h1>
<p align="center" style="width: 500px;">
  <i> Reimplementation of Clustal Software in Python
  </i>
</p>

<p align="center">
    <img alt="Made with Python" src="http://ForTheBadge.com/images/badges/made-with-python.svg">
    <img alt="Made with heart" src="http://ForTheBadge.com/images/badges/built-with-love.svg">
</p>

## Introduction 📚
[Clustal](https://en.wikipedia.org/wiki/Clustal) is a widely-used tool for multiple sequence alignment, originally created by Des Higgins in 1988. It aligns biological sequences like DNA or proteins to help identify similarities and evolutionary relationships. In this project, we've reimplemented Clustal's core algorithm in Python, making it lightweight and compatible with smaller computers. Instead of aiming for the perfect solution, our approach is more heuristic, focusing on speed and efficiency, which works well for limited computational resources.

The process involves three main steps:
1. **Pairwise Alignment**: Align sequences in pairs using the Needleman-Wunsch algorithm.
2. **Guide Tree Creation**: Generate a guide tree based on the pairwise alignment scores.
3. **Multiple Sequence Alignment**: Utilize the guide tree to align all sequences progressively.


## Setup ⚙️
To install PyClustal and its dependencies, you need to perform the following steps:

### Clone the repository

```bash
git clone https://github.com/Essmaw/PyClustal.git
cd pyclustal
```

### Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Create a Conda environment

```bash
conda env create -f pyclustalenv.yml
```

### Activate the Conda environment

```bash
conda activate pyclustalenv
```

## Usage 🚀

### Command Line Interface (CLI) 🖥️

To run PyClustal, you can use the following command:

```bash
python src/pyclustal.py --f [fasta_file_path] --seq-type [seq_type] --sub-matrix [sub_matrix] --gap-open [gap_open] --gap-ext [gap_ext] --job-name [job_name] --tag-log [tag_log]
```

Here's a brief explanation of the arguments:
- `--f` or `--fasta_file_path`: Path to the input FASTA file.
- `--seq-type`: Type of sequences (dna or protein).
- `--sub-matrix`: Substitution matrix to use.
- `--gap-open`: Gap opening penalty.
- `--gap-ext`: Gap extension penalty.
- `--job-name`: Name of the job.
- `--tag-log`: The flag to enable the logging of the pairwise alignment process in precision or not.

Example :

```bash
python src/pyclustal.py --f data/insulins.fasta --seq-type protein --sub-matrix BLOSOM62 --gap-open -5 --gap-ext -1 --job-name aligned_insulins.fasta --tag-log False
```
> This command will align the sequences contained in the file "data/insulins.fasta" using the blosum62 substitution matrix, a gap opening penalty of -5 and a gap extension penalty of -1. The aligned sequences will be saved in the file "aligned_insulins.fasta" in the results folder. And the alignment process of each pair of sequences will not be logged.


### Web Interface 🌐

To run the web interface, you can use the following command:

```bash
gradio_app.py
```

This will run the Gradio app in your web browser.


Enjoy aligning your sequences! 🎉