<!-- # Overview
Here we share code and interactive notebooks that introduce users to the application
of Density-functional Fluctuation Theory (DFFT) to the analysis of residential segregation
through the use of the classic Schelling simulation of residential segregation.

A useful introduction to the application of DFFT onto segregation data is the analysis of simulations with only **two distinct types of agents** shown in the notebook, DFFT_demonstration_Schelling_binary.ipynb, which can be viewed either in a [static online](insert link) version or in an [interactive online](insert link) version. The interactive version allows users to change the parameters of the Schelling simulation to perform their own investigations.

To investigate Schelling models with **three distinct types of agents** (trinary) instead of only two (binary), users can view the "DFFT_demonstration_Schelling_trinary.ipynb" notebook also in a [static](insert link) or [interactive version](insert link). Lastly, if users wish to conduct analysis on the more classic Schelling model where agents move to **vacant sites** instead of switching locations, as was studied in our work [(arXiv:2008.09663)](https://arxiv.org/abs/2008.09663), users can run the "DFFT_demonstration_Schelling_vacancies.ipynb" notebook in either a [static]() or [interactive]() version.

To run simulations on a local machine using Julia, a free open-source language, on a Jupyter notebook users should follow these steps.
1. Download and install [Julia](https://julialang.org/downloads/). This is the language that simulations are written in.
2. Install the necessary additional packages by running the following line of code within a Julia prompt.
3. Steps 1-2 install software needed to run code. To view the interactive notebooks linked above, you must additionally download and install [Jupyter notebook](https://jupyter.readthedocs.io/en/latest/install/notebook-clas  sic.html).
4. Download this repository, open Jupyter (for instance by running `jupyter notebook` from the command prompt), navigate to the downloaded code, and open any of the interactive notebooks. -->

### Schelling agent-based model of segregation implementation for [arXiv:2008.09663](https://arxiv.org/abs/2008.09663)
For the report published in [arXiv:2008.09663](https://arxiv.org/abs/2008.09663), we use code contained in the folder
labeled matlab-implementation. Stay tuned (or contact us) for a faster implementation using open-source software and an interactive walkthrough.

#### To collect statistics from the Schelling simulation as shown in [arXiv:2008.09663](https://arxiv.org/abs/2008.09663) using MATLAB R2020a (tested on both Ubuntu 16.04 and Windows 10), follow these steps below.

1. Download the folder "matlab-implementation" and add the folder to your MATLAB path.
2. Open the script titled "Schelling1.m". To collect the number of statistics used in [arXiv:2008.09663](https://arxiv.org/abs/2008.09663),
set `SampleSize=50000;` at the top of the script and click "Run" from the top menu, This will take hours to run on a standard laptop PC. To skip this step, load the results by running `load("results_schelling1_SampleSize_2500.mat")`. To minimize the size of the data file, this is run with `SampleSize=2500;` One can also download datafile with `SampleSize=50000;` here: https://drive.google.com/file/d/1xg9_GA5avQIbr61nshIdldsq8YBZs5W7/view?usp=sharing 

    "Schelling1.m" generates the statistics of a steady-state ensemble and the initial conditions for an ensemble where the total population of red and blue agents have changed. A snapshot of the simulation look like the following. ![Schelling smapshot](/matlab-implementation/images/Schelling1_1.png)

3. To visualize the DFFT functions for the densities observed in this data, run the script `Visualize_DFFT_functions`. This script takes the steady state statistics from the simulation and uses a least-squares
algorithm to find the best approximation for a global frustration function and spatial Vexation functions according to equation (3) in the main text. It will generate plots like the following. (Note that, due the invariance of the frustration to changes in its overall slope, the plot of the frustration below has an overall slope that, in the text, we move into a constant shift in the Vexation functions. They are equivalent)
![V_blue](/matlab-implementation/images/V_blue.png)
![V_blue](/matlab-implementation/images/V_blue.png)
![frustration_function](/matlab-implementation/images/frustration_function.png)

#### To perturb the steady-state distributions generated above by changing the number of red and blue agents and letting them adjust to a new final steady-state, follow the steps below

4. Open the script title "Schelling2.m" and click "Run" from the top menu. Again, to save time, you can skip this step for a smaller `SampleSize=2500;` by running `load("results_schelling2_SampleSize_2500.mat")`. One can also download datafile with `SampleSize=50000;` here: https://drive.google.com/file/d/1xm3TDmserPSt5wXZbWAI9NyHAnJ-xLTO/view?usp=sharing

    "Schelling2.m" takes the initial conditions from "Schelling1.m" and simulates how the population
    of each modified city will evolve into the future, thus collecting the statistics for the out-of-equilibrium
    time-dependent evolution of the Schelling simulation. This will also plot the observed response of the probability of finding a given density over time (equivalent to figure 4d and 4e).

    ![Schelling2_1](/matlab-implementation/images/Schelling2_1.png)
    ![Schelling2_2](/matlab-implementation/images/Schelling2_2.png)


#### To generate forecasts of the dynamics from step 4, follow these steps below.

5. If you have not run steps (1, 2, and 4), download the folder "matlab-implementation", navigate to that folder in MATLAB and load the results using the command `load(results_schelling2.mat)`

6. Open the script "TPreds.m" and click "Run". This will generate the variable `Omove` that characterizes the time scale between the the Schelling simulation and the scale of the DFFT forecasts.

7. Open the script "TPred.m" and click "Run". This will use the time scaling from step 5 to generate forecasts of the time evolution using the model described in Section 4.1 of the paper. It will produce forecasts like those shown in figure 4b and 4c.

![Schelling2_1](/matlab-implementation/images/Schelling_DFFT_pred1.png)
![Schelling2_2](/matlab-implementation/images/Schelling_DFFT_pred2.png)

8. To generates the Mean Value Estimate (MVE) forecast that quickly forecasts the evolution of the mean of the distribution, open the script "MVE.m" and click "Run". This will generate the time evolution of the mean number of each type of agent across each block as shown below.

![Schelling2_1](/matlab-implementation/images/MVE_blue.png)
![Schelling2_2](/matlab-implementation/images/MVE_red.png)

9. Run the function `[pN1_b,pN2_b,P]=SDFFTPred2(DataSheet,fN1,fN2,s)` with `fN1` and `fN2` set to the desired new total number of blue and red agents respectively in order to find the new steady state probability distributions as described in section 5 of the paper.
