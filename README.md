[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/djmaity/allobench/HEAD?labpath=AlloBench.ipynb)

# AlloBench

The AlloBench is a pipeline to create dataset of proteins with known allosteric and active sites.

## Installation
Install Python 3 and the required dependencies using either of the two methods below. It is highly recommended to create a virtual environment to avoid conflicts with existing packages on the computer.

### Using Anaconda or Miniconda (Recommended)
If not already intalled, install Anaconda or Miniconda following the instructions at https://docs.anaconda.com/anaconda/install/ and https://docs.anaconda.com/miniconda/miniconda-install/, respectively.

Create the allostery virtual environment using the environment.yml file in the repository:

```
conda env create --file=environment.yml
```

### Using Python pip

Python 3 comes preinstalled on Linux and macOS. If pip is not installed please refer the [PIP Documentation](https://pip.pypa.io/en/stable/installation/)
```
python3 -m venv allostery
```
```
conda install pip
pip install -r requirements.txt
```
The nootbooks require an internet connection to automatically obtain the required data in many of its cells. The output is written to the output directory of the repository. The input and output locations can be customized in the notebooks. Please refer to the text cells in each notebook for further instructions. 

## Running AlloBench
The AlloBench.ipynb Jupyter notebook contains the entire pipeline. An internet connection is required to run the notebook as it fetches data from UniProt and RCSB Protein Data Bank (PDB).

### 1. Download AlloSteric Database (ASD) Data
The data from ASD cannot be automatically fetched and must be manually downloaded after requesting access. Download the ASD_Release_202306_XF.tar.gz archive from from https://mdl.shsmu.edu.cn/ASD/module/download/download.jsp?tabIndex=1 and place it in the directory with the notebook. If using mybinder.org upload the file using 

### 2. Launch Jupyter Lab and open AlloBench.ipynb
```
conda activate allostery
jupyter lab AlloBench.ipynb
```

### 3. Run the AlloBench
The easiest way is to run the entire notebook is to select:

Kernel > Restart Kernel and Run All Cells...

Please refer the [JupyterLab User Guide](https://jupyterlab.readthedocs.io/en/latest/getting_started/overview.html) if you encounter any issues

## Output

### Files
* **AlloBench.csv** - The unfiltered AlloBench dataset with allosteic and active site residues
* **ASD_Updated.csv** - The ASD data after updating obsolete PDB and UniProt IDs
* **ASD_Enriched.csv** - The ASD data complemented by additional information from UniProt and PDB

### Directories
* **pdb_structures** - Directory containing downloaded PDB structures
