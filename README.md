# NanoZoomer Connectomics Pipeline: standalone version of the Gauss-Seidel Iteration Scheme
Python standalone version of the Gauss-Seidel Iteration Scheme by [Gaffling et al. (2015)](https://www.ncbi.nlm.nih.gov/pubmed/25312918) to smooth high-frequency distortions across a stack of brain section images. This is a standalone version based on code from our [NanoZoomer Connectomics Pipeline](https://doi.org/10.1101/748376). The essential 'run' function is Nipype ready code; contains all necessary calls as inner functions.

Author: Alexander Woodward, Connectome Analysis Unit, RIKEN CBS, Wako, Japan. Email: alexander.woodward at riken dot jp

## Usage

The code was tested using Mac OSX, Python 3.7 using [Anaconda](https://www.anaconda.com/distribution/), and [ANTs Advanced Normalization Tools 2.1.0](https://github.com/ANTsX/ANTs/releases/tag/v2.1.0). 
1. Setup ANTs and make sure *antsRegistration* is accessible from your shell path.
1. Use `conda create --name <env> --file requirements.txt`, where `<env>` is the name of the environemt, to create a suitable conda environment.
2. Run `python main.py <input_directory>`, where `<input_directory>` is a folder of brain section images in sequence, starting with index 1, e.g. image_001.tif, image_002.tif, image_003.tif, ...
3. Working and output directories will be generated in the same directory as main.py.
## Citation

If you use this code please cite the paper that describes the computational pipeline that it is a part of:

```
@article {Woodward748376,
	author = {Woodward, Alexander and Gong, Rui and Abe, Hiroshi and Nakae, Ken and Hata, Junichi and Skibbe, Henrik and Yamaguchi, Yoko and Ishii, Shin and Okano, Hideyuki and Yamamori, Tetsuo and Ichinohe, Noritaka},
	title = {The NanoZoomer Connectomics Pipeline for Tracer Injection Studies of the Marmoset Brain},
	elocation-id = {748376},
	year = {2019},
	doi = {10.1101/748376},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2019/08/28/748376},
	eprint = {https://www.biorxiv.org/content/early/2019/08/28/748376.full.pdf},
	journal = {bioRxiv}
}
```

