# FoodWebFiveFundamentalWays

When using this code for further research, please cite the original publication: Lever, J.J., Van Nes, E. H., Scheffer, M., & Bascompte, J.  (in press). Five fundamental ways in which complex food webs may spiral out of control. *Ecology Letters*.

## PURPOSE AND SETUP

The purpose of this code is to provide insight in the analysis done for the above-mentioned publication. Adjustments will likely be necessary when using this code for your own research. Please contact the authors for more information.

The file 'MAIN_makeAnalyseFoodWebs.m' is the main function executing three steps from generating and parameterizing food webs (step 1), analysis of the properties of these food web (step 2), and the application of a change in environmental conditions (step 3).
The files 'STEP0_*.m' are used to set parameter ranges and need to be run before execution of the main function.
The folder '+code' contains the Matlab functions supporting the main function (MAIN_makeAnalyseFoodWebs).  
The folder 'parametersets' contains the settings needed to make simulations, i.e. the ranges from which parameters are sampled.
Model-generated networks and the outcome of analysis are stored under 'results'.

## NOTES ON MAIN_makeAnalyseFoodWebs.m

To run the main function: open Matlab, go to the folder where the main function is stored, and type 'MAIN_makeAnalyseFoodWebs(NETnr_min,NETnr_max)' in the command line. This will generate and analyse food webs with the numbers starting with number 'NETnr_min' and ending with 'NETnr_max', e.g. to generate 100 networks with the number 1-100, type 'MAIN_makeAnalyseFoodWebs(1,100)'. Multiple networks (i.e. with the same number) may be generated and stored for different parameter ranges (i.e. specified under 'paramSetNameSeries'). Each network may be subjected to different scenarios of environmental change specified under 'changeSetNameSeries'.

Most of the functions that are called within the main function (MAIN_makeAnalyseFoodWebs) have a similar setup. Model-generated food webs and analysis are loaded from a file when 'replace*=false' (e.g. 'replaceNETWORK=false'). New food webs are generated and analysis are done when 'replace*=true' (e.g. 'replaceNETWORK=true'), or when there is no datafile available. The name of the file under which data are stored is returned by '*File' (e.g. 'networkFile'). Data of three example networks are provided.

For each food web and environmental change scenario, the impact of changes towards multiple sets of random 'final' parameter values (at E=1) is explored of which the number is defined by 'NRchange'. The number of steps taken to increase environmental condition E from 0 to 1 (i.e. towards each of these final parameter sets) by the function 'changeNETWORK' is defined by NRsteps (e.g. when NRsteps=10000, E is increased with 0.0001 at each step. From this we may know that a threshold is passed between, e.g., a value of E=0.0532-0.0533). In addition to this, this code provides a function (i.e. 'changeNETWORKPrescision') which divides the step during which a threshold is passed into a further number of sub-steps defined by NRprecisionSteps (e.g. when NRsteps=10000 and NRsubStepsPrecision=100 and a threshold is passed between a value of E=0.0532-0.0533, E is increased with sub-steps of 0.000001 from the value of E=0.0532 onwards).

## NOTES ON STEP0_parameterSettings.m

Parameter ranges from which parameters are sampled can be specified with the script 'STEP0_parameterSettings.m'. The chosen ranges are stored using the name specified under 'SETNAME' in the folder 'parametersets' and can be loaded by the main function (MAIN_makeAnalyseFoodWebs.m) using the same name. Default parameter settings are already stored and made available under the name 'default_v1'.

## NOTES ON STEP0_envChangeSettings.m

The nature of different environmental change scenarios can be specified using 'STEP0_envChangeSettings.m'. A default scenario is provided under the name 'default_v1'. Please note that, e.g. 'CHANGE_totFeedRate=0' means that, as a part of this scenario, no change is applied to the consumers' total feeding rates.
