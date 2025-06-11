## Petrophysical Joint Inversion Framework
This repository contains a collection of inversion frameworks developed to estimate hydrogeological subsurface properties by jointly inverting seismic refraction and electrical resistivity data. The core concept is based on **petrophysical joint inversion (PJI)**, where geophysical responses are linked through petrophysical laws to shared state variables such as porosity, saturation, temperature, and others.

The framework was originally developed for applications in **partially frozen porous media**, but the modifications available in this repository make it broadly applicable to various environments where seismic and electrical properties are coupled through pore-scale physics.

---

[![License](https://img.shields.io/badge/license-BSD-green)](LICENSE.md)

## Motivation

Standard inversion methods for geophysical data often treat different data types independently and use overly simplified petrophysical models. The PJI framework addresses these limitations by:

- Integrating multiple geophysical data types (e.g., seismic and electrical)
- Linking them through **shared petrophysical parameters** (e.g., porosity, saturation, temperature, cementation exponent)
- Accounting for **nonlinear couplings** and **state dependencies** of physical properties
- Allowing **uncertainty-aware** and **data-driven** estimates of hydrogeological properties

The modular design of the framework supports both well-established petrophysical laws (e.g., Archie's law) and extended formulations (e.g., Waxman–Smits, temperature dependence, CEC estimation), facilitating tailored model development.

## Branch Overview

The repository is organized into several branches, each representing a different stage or variant of the PJI framework. Below is a description of the main branches and their associated research contexts:

- **`wagner2019`**  
  This branch contains the implementation published by Wagner et al. (2019), which uses Archie's law to model the electrical response. It represents the original PJI formulation and is based on a four-phase petrophysical model of unfrozen water, ice, air, and solid matrix.  
  Repository: [four-phase-inversion by Florian Wagner](https://github.com/florian-wagner/four-phase-inversion)

- **`dev`**  
  This branch corresponds to the implementation presented in Steiner et al. (2022), which extends the framework by incorporating the Waxman–Smits (1968) conduction model. This formulation includes **surface conductivity effects** via exchangeable cations and uses **dual-frequency electrical resistivity data** to better constrain these processes. It introduces a frequency-dependent treatment of conductivity, addressing challenges in interpreting electrical responses in clay-rich or saline environments.

- **`dev-estm`**  
  This experimental branch implements a method proposed in Steiner's PhD thesis (2023), where the **cementation exponent \( m \)** in Archie's law is not treated as fixed but is instead **estimated as part of the inversion**. This allows site-specific calibration and adapts the model to variations in pore structure or saturation history, potentially improving model realism in heterogeneous environments.

- **`dev-tempdep`**  
  Focuses on **temperature-dependent electrical conductivity**. This implementation accounts for the influence of temperature on both bulk and surface conduction, which is especially relevant for permafrost or seasonal thaw environments. While it does not estimate the cementation exponent, it introduces temperature-corrected conductivity models that adjust dynamically during the inversion.

- **`dev-estcec`**  
  This branch combines **temperature modeling** with an advanced inversion strategy that directly estimates the **cation exchange capacity (CEC)** from the data. Unlike previous approaches that updated CEC externally, here it is included in the parameter vector, allowing full coupling with other state variables and improving inversion consistency for systems with significant surface conduction (e.g., fine-grained sediments, weathered profiles).


## Structure of this repository

All source code used to generate the results and figures in the paper are in the
`code` folder. A Python library holds the important bits and pieces, which are
resued for calculations and figure generation run in Python scripts. The field
data used in this study are provided in the `data` folder and the sources for
the manuscript text and figures  (LaTeX) are in `manuscript`. See the
`README.md` files in each directory for a full description.

## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/florian-wagner/four-phase-inversion.git

or [download a zip archive](https://github.com/florian-wagner/four-phase-inversion/archive/master.zip).

## Dependencies

You'll need a working Python environment on a Linux machine to run the code.
Other operating systems are generally possible, but have not been tested. The
recommended way to set up your environment is through the [Anaconda Python
distribution](https://www.anaconda.com/download/) which provides the `conda`
package manager. Anaconda can be installed in your user directory and does not
interfere with the system Python installation. The required dependencies are
specified in the file `environment.yml`.

We use `conda` virtual environments to manage the project dependencies in
isolation. Thus, you can install our dependencies without causing conflicts with
your setup (even with different Python versions).

Run the following command in the repository folder (where `environment.yml` is
located) to create a separate environment and install all required dependencies
in it:

    conda env create


## Reproducing the results

Before running any code you must activate the conda environment:

    conda activate four-phase-inversion

This will enable the environment for your current terminal session.
Any subsequent commands will use software that is installed in the environment.

To build the software, produce all results and figures, and compile
the manuscript PDF, run this in the top level of the repository:

    make all

If all goes well, the manuscript PDF will be placed in `manuscript/output`.

You can also run individual steps in the process using the `Makefile`s from the
`code` and `manuscript` folders. See the respective `README.md` files for
instructions.

## Quick example combining the steps above

    unset PYTHONPATH # to avoid conflicts with packages outside the conda env
    git clone https://github.com/florian-wagner/four-phase-inversion.git
    cd four-phase-inversion
    conda env create
    conda activate four-phase-inversion
    cd code
    make build
    make all

## License

All source code is made available under a BSD 3-clause license. You can freely
use and modify the code, without warranty, so long as you provide attribution
to the authors. See `LICENSE.md` for the full license text.

The manuscript itself is an Open Access article distributed under the terms of
the Creative Commons Attribution License
(http://creativecommons.org/licenses/by/4.0/), which permits unrestricted reuse,
distribution, and reproduction in any medium, provided the original work is
properly cited.

## Credits

The software implementation is based on [pyGIMLi](https://www.pygimli.org) (and
its dependencies), which would not exist without the dedication of Carsten
Rücker. This repository is heavily inspired by a [template for
reproducible research papers](https://www.leouieda.com/blog/paper-template.html)
by Leonardo Uieda.
