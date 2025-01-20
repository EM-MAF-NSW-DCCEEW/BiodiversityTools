# BiodiversityTools
A suite of spatial biodiversity modelling tools including software and code developed and maintained by NSW DCCEEW's Metrics and Forecasting Ecological Modelling Team (EMT). 

The EMT are part of the Science Impact and Assessment Branch in the New South Wales (NSW) Department of Climate Change, Energy, the Environment and Water’s (DCCEEW’s) Science and Insights Division. EMT deliver spatial biodiversity modelling methods, tools and information used to assess and improve biodiversity outcomes in NSW, Australia. 

Biodiversity modelling tools including software and code developed and maintained by the EMT and used by DCCEEW in NSW are being published here in the hope that they may have broader application within and beyond NSW. 

All software and code are Copyright(C) 2025 State of New South Wales and Department of Climate Change, Energy, the Environment and Water (DCCEEW) and are published here under version 3 of the GNU General Public Licence.

All software and code are provided in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See version 3 of the GNU General Public License https://www.gnu.org/licenses/ for more details.

## Methods

The software and code implement the following published spatial biodiversity modelling methods and tools:

* Spatial Links Tool (1).

* Spatial Cost Benefit Approach (2).

* Generalised Landscape Connectivity Model (3).

* Biodiversity Forecasting Tool (4).

* Rapid Evaluation of Metapopulation Persistence (5, 6).

Refer to the papers below for technical descriptions of these methods and their application.

These methods have supported state and regional biodiversity assessments including regional conservation initiatives, habitat connectivity, native vegetation management, threatened species and climate impact assessments. Most recently, they have underpinned biodiversity and ecological integrity indicators published under the Biodiversity Indicator Program (BIP) and are currently supporting Saving our Species investment in landscape managed species. 

Adoption of these methods extend beyond NSW with use in external government and research programs. Examples include the Queensland Government’s Spatial Biocondition method, and the Australian Department of Climate Change, Energy, the Environment and Water’s National Connectivity Index; as well as collaborative research with CSIRO, the University of New England, Macquarie University and the University of Queensland.


## Software and code

These related and interdependent tools are packaged together in a single ‘Biodiversity Tools’ software package. The software is written in C++ and makes use of Nvidia’s CUDA graphics processing capabilities. When built, the software package contains binary executable console applications for the above methods and tools, along with static and dynamic software libraries required to implement the methods. 

The tools operate on spatial raster data stored in ESRI's floating point file format (.flt, .hdr) however GeoTiff support is gradually being implemented across tools along with other improvements as time permits. Each of the console applications process plain text parameter (.par) files defining inputs, outputs and parameters. Template parameter files containing parameter descriptions can be created by passing the applications a parameter filename that doesn’t already exist. 

Whilst originally developed and maintained for the Windows OS environment only, a cmake build system (https://cmake.org/) has recently been adopted and the software has been modified to build and run in both Windows and Linux environments. Development and testing in Linux is still ongoing. 
The software is written using features from the ISO C++17 standard and relies on the following 3rd party software libraries for various functionality.

* GDAL https://gdal.org/index.html https://gdal.org/license.html

* Nvidia https://developer.nvidia.com/cuda-toolkit https://docs.nvidia.com/cuda/eula/index.html

All tools except Spatial Links require a CUDA compatible Nvidia graphics processing (GPU) device (graphics card). A list of compatible devices can be found at:

https://developer.nvidia.com/cuda-gpus

## Copyright

All software and code are Copyright(C) 2025 State of New South Wales and Department of Climate Change, Energy, the Environment and Water (DCCEEW) and are published here under version 3 of the GNU General Public Licence.


## Licence

All software and code are made available under the GNU General Public Licence V3.0

https://www.gnu.org/licenses/gpl-3.0.txt


## References

(1) The spatial links tool: Automated mapping of habitat linkages in variegated landscapes. 
Michael Drielsma, Glenn Manion, Simon Ferrier, 2007. 
Ecological Modelling, Volume 200, Issues 3–4, Pages 403-411, ISSN 0304-3800. 
https://doi.org/10.1016/j.ecolmodel.2006.08.017

(2) A raster-based technique for analysing habitat configuration: The cost–benefit approach. 
Michael Drielsma, Simon Ferrier, Glenn Manion, 2007. 
Ecological Modelling, Volume 202, Issues 3–4, Pages 324-332, ISSN 0304-3800. 
https://doi.org/10.1016/j.ecolmodel.2006.10.016

(3) General Landscape Connectivity Model (GLCM): a new way to map whole of landscape biodiversity functional connectivity for operational planning and reporting. 
Michael J. Drielsma, Jamie Love, Subhashni Taylor, Rajesh Thapa, Kristen J. Williams, 2022. 
Ecological Modelling, Volume 465, 109858, ISSN 0304-3800. 
https://doi.org/10.1016/j.ecolmodel.2021.109858

(4) The Biodiversity Forecasting Toolkit: Answering the ‘how much’, ‘what’, and ‘where’ of planning for biodiversity persistence. 
Michael Drielsma, Simon Ferrier, Gary Howling, Glenn Manion, Subhashni Taylor, Jamie Love, 2014. 
Ecological Modelling, Volume 274, Pages 80-91, ISSN 0304-3800. 
https://doi.org/10.1016/j.ecolmodel.2013.11.028

(5) Rapid evaluation of metapopulation persistence in highly variegated landscapes. 
Michael Drielsma, Simon Ferrier, 2009. 
Biological Conservation, Volume 142, Issue 3, Pages 529-540, ISSN 0006-3207. 
https://doi.org/10.1016/j.biocon.2008.11.018

(6) An equitable method for evaluating habitat amount and potential occupancy. 
Drielsma M, & Love J, 2021. 
Ecological Modelling, Vol.440, pp.109388. 
https://doi.org/10.1016/j.ecolmodel.2020.109388
