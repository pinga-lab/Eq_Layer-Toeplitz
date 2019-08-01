# Fast equivalent layer for gravity data processing with BTTB systems

by
[Diego Takahashi](http://www.pinga-lab.org/people/tomazella.html)<sup>1</sup>,
[Vanderlei C. Oliveira Jr.](http://www.pinga-lab.org/people/oliveira-jr.html)<sup>1</sup> and
[Valéria C. F. Barbosa](http://www.pinga-lab.org/people/oliveira-jr.html)<sup>1</sup>

<sup>1</sup>[Observatório Nacional](http://www.on.br/index.php/pt-br/)

This work will be submitted for publication in
[*Geophysics*](https://seg.org/Publications/Journals/Geophysics).


## Abstract

We have developed an efficient and very fast equivalent layer technique for gravity data 
processing by modifying an iterative method grounded on excess mass constraint that does not 
require the solution of linear systems. Taking advantage of the properties related to the 
symmetric Block-Toeplitz Toeplitz-block (BTTB) and Block-Circulant Circulant-Block (BCCB) 
matrices, that raises when regular grids of observation points and equivalent sources 
(point masses) are used to set up a fictitious equivalent layer, we have developed an 
algorithm which greatly reduces the number of flops and memory RAM necessary to estimate 
of a 2D mass distribution over the equivalent layer. The algorithm is based on the structure 
of symmetric BTTB matrix, where all its elements consist of the elements of the first row, 
which in turn can be embedded into a symmetric BCCB matrix. Likewise, only the first row of 
the BCCB matrix is needed to reconstruct the full matrix completely. From the first column 
of BCCB matrix, its eigenvalues can be calculated using the fast Fourier transform, which 
can be used to readily compute the matrix-vector product. As a result, our method is 
efficient to process very large datasets using either fine- or mid-grid meshes. The larger 
the dataset, the faster and more efficient our method becomes compared to the available 
equivalent-layer techniques. Synthetic tests demonstrate the ability of our method to 
satisfactorily upward- and downward-continuing the gravity data.

![](manuscript/Fig/sensibility_grav_mod.png)

**Figure 1:** *symmetric block-toeplitz toeplitz-block structure of the gravimetric sensibility 
matrix that arises when using the equivalent-layer technique for regular grids of data and 
equivalent sources.*

![](manuscript/Fig/float.png)

**Figure 2:** *floating points to estimate the parameter vector using the fast equivalent 
layer with Siqueira et al.'s method and our approach versus the numbers of observation 
points varyig from $N = 5000$ to $N = 1000000$ with $50$ iterations. The number of operations 
is drastically decreased.*


![](manuscript/Fig/time_comparison.png)

**Figure 3:** *time necessary to run 50 iterations of the Siqueira et al.'s method and the 
one presented in this work. With the limitation of $16$ Gb of memory RAM in our system, we 
chose to test only up to $22500$ obervation points.*


## Reproducing the results

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/DiegoTaka/Eq_Layer-Toeplitz.git


All source code used to generate the results and figures in the paper are in
the `code` folder. The sources for the manuscript text and figures are in `manuscript`.
See the `README.md` files in each directory for a full description.

The calculations and figure generation are all run inside
[Jupyter notebooks](http://jupyter.org/).
You can view a static (non-executable) version of the notebooks in the
[nbviewer](https://nbviewer.jupyter.org/) webservice:

http://nbviewer.jupyter.org/github/DiegoTaka/Eq_Layer-Toeplitz

See sections below for instructions on executing the code.


### Setting up your environment

You'll need a working Python **2.7** environment with all the standard
scientific packages installed (numpy, scipy, matplotlib, etc).  The easiest
(and recommended) way to get this is to download and install the
[Anaconda Python distribution](http://continuum.io/downloads#all).
Make sure you get the **Python 2.7** version.

Use `conda` package manager (included in Anaconda) to create a
[virtual environment](https://conda.io/docs/using/envs.html) with
all the required packages installed.
Run the following command in this folder (where `environment.yml`
is located):

    conda env create

To activate the conda environment, run

    source activate ellipsoids

or, if you're on Windows,

    activate ellipsoids

This will enable the environment for your current terminal session.
After running the code, deactivate the environment with the following
commands:

    source deactivate

or, if you're on Windows,

    deactivate


**Windows users:** We recommend having a bash shell and the `make` installed
to run the code, produce the results and check the code. You may download the
[*Git for Windows*](https://git-for-windows.github.io/) and the
[*Software Carpentry Windows Installer*](https://github.com/swcarpentry/windows-installer/releases).


### Running the code

To execute the code in the Jupyter notebooks, you must first start the
notebook server by going into the repository folder and running:

    jupyter notebook

Make sure you have the `conda` environment enabled first.

This will start the server and open your default web browser to the Jupyter
interface. In the page, go into the `code` folder and select the
notebook that you wish to view/run.

The notebook is divided into cells (some have text while other have code).
Each cell can be executed using `Shift + Enter`.
Executing text cells does nothing while executing code cells runs the code
and produces it's output.
To execute the whole notebook, run all cells in order or use "Cell -> Run All"
from the menu bar.

## License

All source code is made available under a BSD 3-clause license.  You can freely
use and modify the code, without warranty, so long as you provide attribution
to the authors.  See `LICENSE.md` for the full license text.

The manuscript text is not open source. The authors reserve the rights to the
article content.