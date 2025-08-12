# RLSRP_2025

This is the GitHub repository containing all code for the preprint  
**Reference-based variant detection with varseek**  
by *Joseph Matthew Rich, Laura Luebbert, Delaney Sullivan, Reginaldo Rosa, and Lior Pachter*.

**varseek:** [GitHub - pachterlab/varseek](https://github.com/pachterlab/varseek.git)
**Examples for getting started:** [GitHub - pachterlab/varseek-examples](https://github.com/pachterlab/varseek-examples.git)

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

`notebook/`: Jupyter notebooks to reproduce each main and supplemental figure, named according to the figure number.
`RLSRP_2025/`: Core functions used within notebooks
`scripts/`: Long scripts for generating variant indices or running variant calling with varseek and other tools

## License  
This project is licensed under the **BSD 2-Clause License**. See the [LICENSE](LICENSE) file for details.

---

For any issues or contributions, feel free to open a pull request or issue in this repository.
