# Spin-chain models simulator

Worm-like spin chain models — a 2020 summer project with Abhishek Chaudhuri.

> [!IMPORTANT]
> **Status: placeholder.** The original 2020 code and results have not been added to this repository yet. For now this README sets out the project's scope, and the sections marked **TODO** will be filled in once that material is uploaded. No results from the project are reported here yet.

## Project scope

- **Project:** summer project, 2020
- **With:** Abhishek Chaudhuri
- **Topic:** worm-like spin chain models
- **Institution / programme:** **TODO**
- **Duration:** **TODO**
- **Language(s) and tools:** **TODO**

> **TODO:** Summarise the goals of the project and the questions it set out to answer.

## Background

*General textbook background for readers new to the topic. It does not describe what this project's code implements; that belongs under [Models and methods](#models-and-methods).*

The worm-like chain (WLC), or Kratky–Porod model [1], describes a semiflexible polymer as an inextensible curve of contour length $`L`$. Writing $`\hat{\mathbf t}(s)`$ for the unit tangent at arc length $`s`$, the energy penalises bending:

```math
E_{\mathrm{bend}} = \frac{\kappa}{2} \int_0^L \left| \frac{\mathrm{d}\hat{\mathbf t}}{\mathrm{d}s} \right|^2 \mathrm{d}s ,
```

where $`\kappa`$ is the bending rigidity. Tangent directions decorrelate over the persistence length $`\ell_p`$, which in three dimensions is $`\ell_p = \kappa / k_B T`$.

Splitting the chain into $`N`$ segments of length $`b`$ with unit tangents $`\hat{\mathbf t}_i`$, and replacing the derivative by a finite difference, gives, up to an additive constant,

```math
E = -J \sum_{i=1}^{N-1} \hat{\mathbf t}_i \cdot \hat{\mathbf t}_{i+1} , \qquad J = \frac{\kappa}{b} .
```

This has the form of a one-dimensional classical ferromagnetic spin chain, with the tangents as unit spins: a Heisenberg chain in three dimensions [2] and an XY (planar-rotator) chain in two. A force $`\mathbf f`$ pulling on the chain ends adds the term $`-\mathbf f \cdot \mathbf R = -b \, \mathbf f \cdot \sum_i \hat{\mathbf t}_i`$, where $`\mathbf R`$ is the end-to-end vector, and this term acts on the spins like an external magnetic field. The chain's response to such a force is the WLC force–extension relation, widely used for stretched DNA [3].

## Models and methods

> **TODO:** Describe the specific model(s) studied in the project (definitions, dimensionality, boundary conditions, parameter values), the numerical or analytical methods used, and the quantities computed.

## Repository contents

> **TODO:** Add the original 2020 code and list its main files and folders here.

## Requirements and installation

> **TODO:** List the language version and dependencies the code needs (for example a `requirements.txt` if it is written in Python) and how to install them.

## Reproducing the key figures

> **TODO:** Add a script or notebook that regenerates each key figure from the original code, with the exact command and parameters for each one.

## Results

> **TODO:** Add the key figures with a short description of each. This section will only contain results produced by the original project code.

## Tests

> **TODO:** Add tests for the core functions and the command that runs them.

## References

Background references for the general physics above. The references used in the project itself are still to be added (**TODO**).

1. O. Kratky and G. Porod, "Röntgenuntersuchung gelöster Fadenmoleküle," *Recl. Trav. Chim. Pays-Bas* **68**, 1106–1122 (1949). [doi:10.1002/recl.19490681203](https://doi.org/10.1002/recl.19490681203)
2. M. E. Fisher, "Magnetism in one-dimensional systems—the Heisenberg model for infinite spin," *Am. J. Phys.* **32**, 343–346 (1964). [doi:10.1119/1.1970340](https://doi.org/10.1119/1.1970340)
3. J. F. Marko and E. D. Siggia, "Stretching DNA," *Macromolecules* **28**, 8759–8770 (1995). [doi:10.1021/ma00130a008](https://doi.org/10.1021/ma00130a008)

## Acknowledgements

This project was carried out in summer 2020 with Abhishek Chaudhuri.

> **TODO:** Add affiliations, the programme, and any funding.

## License

Released under the [MIT License](LICENSE).
