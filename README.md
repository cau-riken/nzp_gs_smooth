# NanoZoomer Connectomics Pipeline: standalone version of the Gauss-Seidel Iteration Scheme
Python standalone version of the Gauss-Seidel Iteration Scheme by [Gaffling et al. (2015)](https://www.ncbi.nlm.nih.gov/pubmed/25312918) to smooth high-frequency distortions across a stack of brain section images. This is a standalone version based on code from our [NanoZoomer Connectomics Pipeline](https://doi.org/10.1101/748376). The essential 'run' function is Nipype ready code (all of the inner functions that are used are placed inside it).

Author: Alexander Woodward, Connectome Analysis Unit, RIKEN CBS, Wako, Japan. Email: alexander.woodward at riken dot jp

If you use this code for your exciting research, please cite the paper that describes the computational pipeline that this code is a part of:

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

