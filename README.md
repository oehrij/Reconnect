# Reconnect
Rapid Evaluation of Multispecies Connectivity (Reconnect) R-tool.

This tool enables the rapid computation of connectivity indicators simultaneously for multiple species dispersal and habitat needs, across large regions of interest. Input data need to be in raster or shape file format and describe species habitat, land-cover, or protected area distribution. Cost distances can be computed in case a resistance raster layer is given.
The Reconnect R-tool application is demonstrated in a recent study: [Rapid evaluation of habitat connectivity change to safeguard multispecies persistence in human-transformed landscapes](https://www.biorxiv.org/content/10.1101/2023.11.23.568419v1). This work was developed as part of a postdoc project collaboration of McGill University and the environmental research firm [Habitat](https://www.habitat-nature.com/) co-supervised by [Dr. Andrew Gonzalez](https://www.thegonzalezlab.org/) and [Dr. Brian Leung](https://leung-lab.github.io/leunglab/). Have a look at the [Reconnect documentation](https://github.com/oehrij/Reconnect/blob/main/doc/Reconnect_approach.pdf) for a more detailed description.

![Slide2](https://github.com/oehrij/Reconnect/assets/78751500/aff1c93e-cf18-4d70-8393-c9edc3354cf3)


## Installation 
```r
#install_packages("devtools")
library(devtools)
install_github("oehrij/Reconnect",build_vignettes = TRUE)
library(Reconnect)
```


## License

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0). See the [LICENSE](LICENSE) file for details.

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)


## Feedback
Please write me an email at: jacqueline.oehri@gmail.com


## Reconnect_functions
Reconnect - core, wrapper and summary functions, as well as general functions for processing geographic data.

![Slide3](https://github.com/oehrij/Reconnect/assets/78751500/3745bdf7-de6c-4b7a-8b4f-4e149da0f675)


## CONN_functions
Functions to compute a set of complementary connectivity-indicators at pixel-, patch- and landscape-level that were selected based on a literature review by Eluna Touratier and Jacqueline Oehri, within the framework of the Postdoc project: "Linking multispecies connectivity modelling and ecosystem services in the context of landscape urbanization" (see above). Type ```help(package = 'Reconnect', help_type = 'html')``` for a demonstration of the functions in the 'CONN_MPC_vignette'.

![Slide4](https://github.com/oehrij/Reconnect/assets/78751500/16dd39f4-f783-409a-92ae-284b2bfb5480)


## MPC_functions
The MPC_functions can derive metapopulation capacity (MPC) and -related indices from binary habitat distribution maps (shapefile, raster file or patch-area vector/patch-distance matrix). They additionally enable the individual setting of species-specific parameters, such as dispersal capacity and distance-decay function. 
The implemented MPC formula is based on [Hanski and Ovaskainen 2000](https://doi.org/10.1038/35008063), [Hanski 1994](151-162.https://doi.org/10.2307/5591), as well as the modifications suggested by [Schnell et al. 2013](https://doi.org/10.1111/cobi.12047). The R-code is largely based on the R-codes of [Huang et al. 2020](https://doi.org/10.1111/cobi.13364) and [Strimas-Mackey & Brodie 2018](https://doi.org/10.1002/eap.1739). Have a look at the [Full MPC report](https://oehrij.shinyapps.io/MPC_report/) for more details.

![Slide6](https://github.com/oehrij/Reconnect/assets/78751500/c2deca0f-f0bc-410f-81e3-094ceb80544d)


## NLMC_functions
Functions for generating simulated and neutral landscapes.

![Slide7](https://github.com/oehrij/Reconnect/assets/78751500/91203fd8-d2e2-4348-9317-09056ccaa183)


## References
Albert, C. H., Rayfield, B., Dumitru, M., & Gonzalez, A. (2017). Applying network theory to prioritize multispecies habitat networks that are robust to climate and land‐use change. Conservation Biology, 31(6), 1383-1396.

Brandes, U. A faster algorithm for betweenness centrality*. J. Math. Sociol. 25, 163–177 (2001).

Hanski, Ilkka, and Otso Ovaskainen. 2000. The Metapopulation Capacity of a Fragmented Landscape. Nature 404 (6779): 755-58.  https://doi.org/10.1038/35008063

Hanski, I. (1994). A practical model of metapopulation dynamics. Journal of animal ecology, 151-162.https://doi.org/10.2307/5591

Huang, R., Pimm, S. L., & Giri, C. (2020). Using metapopulation theory for practical conservation of mangrove endemic birds. Conservation Biology, 34(1), 266-275. https://doi.org/10.1111/cobi.13364

McGarigal, K., Cushman, S.A., and Ene E. 2012. FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical and Continuous Maps. Computer software program produced by the authors at the University of Massachusetts, Amherst. Available at the following website: https://www.umass.edu/landeco/

Minor, E. S. & Urban, D. L. A graph-theory framework for evaluating landscape connectivity and conservation planning. Conserv. Biol. 22, 297–307 (2008).

Saura, S. & Martínez-Millán, J. Landscape patterns simulation with a modified random clusters method. Springer Science and Business Media LLC 15, 661–678 (2000).

Saura, S., & Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and urban planning, 83(2-3), 91-103.

Saura, S., & Rubio, L. (2010). A common currency for the different ways in which patches and links can contribute to habitat availability and connectivity in the landscape. Ecography, 33(3), 523-537.

Saura, S., Estreguil, C., Mouton, C., & Rodríguez-Freire, M. (2011). Network analysis to assess landscape connectivity trends: application to European forests (1990–2000). Ecological Indicators, 11(2), 407-416.

Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they?. Ecological indicators, 76, 144-158.

Schnell, Jessica K., Grant M. Harris, Stuart L. Pimm, and Gareth J. Russell. 2013. Estimating Extinction Risk with Metapopulation Models of Large-Scale Fragmentation. Conservation Biology 27(3): 520-30  https://doi.org/10.1111/cobi.12047

Sciaini, M., Fritsch, M., Scherer, C. & Simpkins, C. E. NLMR andlandscapetools : An integrated environment for simulating and modifying neutral landscape models in R. Methods Ecol. Evol. 9, 2240–2248 (2018).

Stott, I., Townley, S., Carslake, D., & Hodgson, D. J. (2010). On reducibility and ergodicity of population projection matrix models. Methods in Ecology and Evolution, 1(3), 242-252.

Strimas-Mackey, M., & Brodie, J. F. (2018). Reserve design to optimize the long-term persistence of multiple species. Ecological Applications, 28(5), 1354-1361. https://doi.org/10.1002/eap.1739; https://github.com/mstrimas/metacapa  
