# NanoZoomer Connectomics Pipeline: standalone version of the Gauss-Seidel Iteration Scheme
Python standalone version of the Gauss-Seidel Iteration Scheme by [Gaffling et al. (2015)](https://www.ncbi.nlm.nih.gov/pubmed/25312918) to smooth high-frequency distortions across a stack of brain section images. This is a standalone version based on code from our [NanoZoomer Connectomics Pipeline](https://doi.org/10.1101/748376). The essential 'run' function is Nipype ready code; contains all necessary calls as inner functions.

Author: Alexander Woodward, Connectome Analysis Unit, RIKEN CBS, Wako, Japan. Email: alexander.woodward at riken dot jp

## Usage

The code was tested using Mac OSX, Python 3.7 using [Anaconda](https://www.anaconda.com/distribution/), and [ANTs Advanced Normalization Tools 2.1.0](https://github.com/ANTsX/ANTs/releases/tag/v2.1.0). 
1. Setup ANTs and make sure *antsRegistration* is accessible from your shell path.
1. Use `conda create --name <env> --file requirements.txt`, where `<env>` is the name of the environment, to create a suitable conda environment.
2. Run `python main.py <input_directory>`, where `<input_directory>` is a folder of brain section images in sequence, starting with index 1, e.g. image_001.tif, image_002.tif, image_003.tif, ...
3. Working and output directories will be generated in the same directory as main.py.
## Citation

If you use this code please cite the paper that describes the computational pipeline that it is a part of:

```
@Article{Woodward2020,
author={Woodward, Alexander
and Gong, Rui
and Abe, Hiroshi
and Nakae, Ken
and Hata, Junichi
and Skibbe, Henrik
and Yamaguchi, Yoko
and Ishii, Shin
and Okano, Hideyuki
and Yamamori, Tetsuo
and Ichinohe, Noritaka},
title={The NanoZoomer artificial intelligence connectomics pipeline for tracer injection studies of the marmoset brain},
journal={Brain Structure and Function},
year={2020},
month={May},
day={04},
abstract={We describe our connectomics pipeline for processing anterograde tracer injection data for the brain of the common marmoset (Callithrix jacchus). Brain sections were imaged using a batch slide scanner (NanoZoomer 2.0-HT) and we used artificial intelligence to precisely segment the tracer signal from the background in the fluorescence images. The shape of each brain was reconstructed by reference to a block-face and all data were mapped into a common 3D brain space with atlas and 2D cortical flat map. To overcome the effect of using a single template atlas to specify cortical boundaries, brains were cyto- and myelo-architectonically annotated to create individual 3D atlases. Registration between the individual and common brain cortical boundaries in the flat map space was done to absorb the variation of each brain and precisely map all tracer injection data into one cortical brain space. We describe the methodology of our pipeline and analyze the accuracy of our tracer segmentation and brain registration approaches. Results show our pipeline can successfully process and normalize tracer injection experiments into a common space, making it suitable for large-scale connectomics studies with a focus on the cerebral cortex.},
issn={1863-2661},
doi={10.1007/s00429-020-02073-y},
url={https://doi.org/10.1007/s00429-020-02073-y}
}


```

