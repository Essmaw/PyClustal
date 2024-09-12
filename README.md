<h1 align="center">
  <img style="vertical-align:middle; width:90%; position:fixed;"
  src="/data/img/banner.png">
</h1>
<p  align="center">
  <img alt="Static Badge" src="https://img.shields.io/badge/Built_with_science_and_%E2%9D%A4%EF%B8%8F-%23c4e4ff?style=flat-square&logoColor=%23ffcbeb&color=%23c4e4ff">
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
python src/pyclustal.py --f [fasta_file_path] --seq-type [seq_type] --sub-matrix [sub_matrix] --gap-open [gap_open] --gap-ext [gap_ext] --job-name [job_name] --tag-log [tag_log] --output-format [output_format]
```

Here's a brief explanation of the arguments:
- `--f`: Path to the input FASTA file.
- `--seq-type` (optionnal): Type of sequences (dna or protein).  Default is protein.
- `--sub-matrix` (optionnal): Substitution matrix to use. Default is BLOSUM62.
- `--gap-open` (optionnal): Gap opening penalty. Default is -5.
- `--gap-ext` (optionnal): Gap extension penalty. Default is -1.
- `--job-name` (optionnal): Name of the job. Default is name of the input file with '_aligned' appended.
- `--tag-log` (optionnal): The flag to enable the logging of the pairwise alignment process in precision or not. Default is False.
- `--output-format` (optionnal): The format of the output file (fasta or clustal). Default is clustal.

Example :

```bash
python src/pyclustal.py --f data/fasta_files/insulins.fasta --seq-type protein --sub-matrix BLOSOM62 --gap-open -5 --gap-ext -1 --job-name aligned_insulins.fasta --tag-log False  --output-format clustal
```

💡 This command will align the sequences contained in the file "data/fasta_files/insulins.fasta" using the blosum62 substitution matrix, a gap opening penalty of -5 and a gap extension penalty of -1. The aligned sequences will be saved in the file "aligned_insulins.clustal" in the results folder. And the alignment process of each pair of sequences will not be logged.


### Web Interface 🌐

To run the web interface, you can use the following command:

```bash
gradio_app.py
```

This will run the Gradio app in your web browser.


## Results 📊

The alignment results are provided in both CLUSTAL and FASTA formats. You can find the results in the [results folder](https://github.com/Essmaw/PyClustal/tree/main/results). These alignments are derived from the FASTA files located in the [data folder](https://github.com/Essmaw/PyClustal/tree/main/data/fasta_files). 

Below are the commands used to obtain the results for each file. Each command aligns the sequences contained in the specified FASTA file, using the appropriate substitution matrix and gap penalties. The aligned sequences will be saved in the results folder. If you prefer FASTA format output, simply add the `--output-format fasta` option to the command.

### Commands for Alignment

- **For `dna.fasta`**:
  ```bash
  python src/pyclustal.py --f data/fasta_files/dna.fasta --seq-type dna --sub-matrix NUC.4.4
  ```

- **For `insulins.fasta`**:
  ```bash
  python src/pyclustal.py --f data/fasta_files/insulins.fasta
  ```

- **For `p53.fasta`**:
  ```bash
  python src/pyclustal.py --f data/fasta_files/p53.fasta
  ```

- **For `zinc_finger.fasta`**:
  ```bash
  python src/pyclustal.py --f data/fasta_files/zinc_finger.fasta
  ```



Enjoy aligning your sequences! 🎉