# RLSRP_2025

This is the GitHub repository containing all code for the preprint  
**Creating Pan-Cancer Transcriptomic Mutational Signatures Using Varseek**  
by *Joseph Matthew Rich, Laura Luebbert, Delaney Sullivan, Reginaldo Rosa, and Lior Pachter*.

📄 **Full preprint:** [LINK]  

🔗 **Learn more about Varseek:** [GitHub - pachterlab/varseek](https://github.com/pachterlab/varseek.git)

## Getting Started

To run the code in this repository, follow these steps:

```sh
git clone https://github.com/pachterlab/RLSRP_2025.git
cd RLSRP_2025
```

We recommend using an environment manager such as conda. Some additional non-python packages must be installed for full functionality. If using conda (recommended), simply run the following:

```sh
conda env create -f environment.yml
conda activate RLSRP_2025
```

Otherwise, install these packages manually as-needed (see environment.yml for the list of packages and recommended versions).

Once the environment is set up, install the repository as a package.

```sh
pip install .
```

---

## Repository Contents

### 📓 Notebooks  
Jupyter notebooks to reproduce each main and supplemental figure, named according to the figure number. Each notebook can also be run independently in a Google Colab environment.
- `notebooks/`
  - `1.ipynb`
  - `2.ipynb`
  - `3.ipynb`
  - `4.ipynb`
  - `5.ipynb`
  - `S1.ipynb`
  - `S2.ipynb`
  - `S3.ipynb`
  - `S4.ipynb`
  - `S5.ipynb`

### 🛠 RLSRP_2025
Core functions used within the notebooks
- `RLSRP_2025/`
    - `logger_utils.py`
    - `seq_utils.py`
    - `visualization_utils.py`

### 🛠 Scripts  
Python scripts called within some notebooks:
- `scripts/`

### ⚙️ Environment Setup  
Files to reproduce the coding environment:
- `environment.yml` – Conda and pip package dependencies  
- `requirements.txt` – Pip package dependencies (used within `environment.yml`)  

### License  
This project is licensed under the **BSD 2-Clause License**. See the [LICENSE](LICENSE) file for details.

---

For any issues or contributions, feel free to open a pull request or issue in this repository.
