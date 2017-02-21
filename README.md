# flexsimple
-Description-
The FLEX-based model (Fenicia et al., 2008) consists of five storage components. First, a snow routine has to be run, which is a simple degree-day module, similar as used in, for example, HBV (Bergström, 1976). After the snow routine, the precipitation enters the interception reservoir. Here, water evaporates at potential rates or, when exceeding a threshold, directly reaches the soil moisture reservoir. The soil moisture routine is modelled in a similar way as the Xinanjiang model (Zhao, 1992). Briefly, it contains a distribution function that determines the fraction of the catchment where the storage deficit in the root zone is satisfied and that is therefore hydrologically connected to the stream and generating storm runoff. From the soil moisture reservoir, water can further vertically percolate down to recharge the groundwater or leave the reservoir through transpiration. Transpiration is a function of maximum root zone storage Su,max and the actual root zone storage, similar to the functions described by Feddes et al. (1978). Water that cannot be stored in the soil moisture storage then is split into preferential percolation to the groundwater and runoff generating fluxes that enter a fast reservoir, which represents fast responding system components such as shallow subsurface and overland flow.

snow
     ||
     \/
interception
     ||
     \/
     Su  -> Sf
     ||
     \/
     Ss

For optimizing the following algorithms are implemented: 



-Requirements-
Compiling:
The code has been tested and compiled previously with gfortran. A makefile has been provided, to compile:
-Type "Make" in a terminal window
-Type "Make clean" to remove the old compiled files
Windows binaries:


-Usage-
Namelist for parameters (param.nml) and configuration (config.nml) are provided, where setting can be adjusted. Input files are required with the names Etp.txt, Prec.txt, Qobs.txt, Temp.txt, Date.txt (DD-MM-YYYY). An ascii-file in esri-format containing a DEM is needed for snow modelling. 

-References-
Nijzink, R., Hutton, C., Pechlivanidis, I., Capell, R., Arheimer, B., Freer, J., Han, D., Wagener, T., McGuire, K., Savenije, H., and Hrachowitz, M.: The evolution of root-zone moisture capacities after deforestation: a step towards hydrological predictions under change?, Hydrol. Earth Syst. Sci., 20, 4775-4799, doi:10.5194/hess-20-4775-2016, 2016. 


