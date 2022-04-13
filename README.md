# Global vs Local
Repository of code for establishing how the global method compares to the local method for the paper "A global method for simulating intracellular signaling reduces computational time in multiscale agent-based models with systems biology applications"

## Reproducing Figures

### Figure 3 (Global method is comparable to local method but with substantially reduced wall time)

1. Run_Figure3.m: produces panels A-E. (this script also produces some panels for Figure 4)
2. Run_Figure3FG.m: produces panels F-G (data also used to create some SI figs)

The files beginning with "Test_" can be run to create similar data to what is in the figures.

### Figure 4 (Noise in the local method is negligible and follows predictable patterns)

1. Run_Figure3.m: produces panels A-B.
2. Run_Figure4CD.m: produces panels C-D. (also makes the same plot at other time points)

### Figure 5 (Global method agrees with local method in growth dynamics and spatial distributions without aIL-6R)

1. Run_Figure5.m: produces all panels.

### Figure 6 (Global method agrees with local method in growth dynamics and spatial distributions with aIL-6R)

1. Run_Figure6.m: produces all panels.

### SI Figures

1. Run_varyPK.m: run simulations of FGFR3 model with varying PK parameters.
2. Run_varyReactionRates.m: run simulations of FGFR3 model with varying reaction rate parameters.
3. Run_varyInflux.m: run simulations of FGFR3 model with varying vasculature assumptions.
4. Run_AbDynamics.m: run simulations of FGFR3 model assuming aFGFR3 is an antibody instead of an SMI

## All files

### Scripts

| Name | Description |
| - | - |
|  Run_AbDynamics.m | Runs FGFR3 model treating the SMI as an antibody. |
| Run_Figure3.m | Runs FGFR3 model to compare the local and global methods at the base parameter values. Also plots Figure 3A-E and Figure 4A-B. |
| Run_Figure3FG.m | Runs FGFR3 model with varying initial numbers of cells using both local and global methods to see how the wall time increases with more cells. Also plots Figure 3F-G. |
| Run_Figure4CD.m | Runs FGFR3 model using only the local method and will update a single figure with a 3D heatmap of drug concentration in the TME. Two of these are used for Figure 4C-D. |
| Run_Figure5.m | Runs IL6 model without aIL-6R to compare the two methods. Also plots Figure 5. |
| Run_Figure6.m | Runs IL6 model with aIL-6R to compare the two methods. Also plots Figure 6. |
| Run_varyInflux.m | Runs FGFR3 model while varying the vasculature assumptions about the tumor. Produces the data used to generate some of the figures in the supplement. |
| Run_varyPK.m | Runs FGFR3 model while varying the PK parameters. Produces the data used to generate some of the figures in the supplement. |
| Run_varyReactionRates.m | Runs FGFR3 model while varying the reaction rate parameters. Produces the data used to generate some of the figures in the supplement. |

### Functions

| Name | Description |
| - | - |
| agentODE_FGFR3.m | Takes in the state variables for an agent (or vector of agents) along with parameter values and outputs the rate of change according to the reaction equations in the FGFR3 model. |
| agentODE_IL6.m | Takes in the state variables for an agent (or vector of agents) along with parameter values and outputs the rate of change according to the reaction equations in the IL6 model. |
| agentODEJacobian_FGFR3.m | Takes in the state variables for an agent (or vector of agents) along with parameter values as well as the constant entries in the Jacobian of the ODE function and outputs the Jacobian. Used in `substrateSolver_FGFR3.m` to update the reaction equations rather than using an RK method or direct Euler. |
| allCombos.m | Helper function that basically takes the _n_ outputs of `ndgrid` and makes each a column of a single array. Used to determine neighbors in lattice. See, for example, https://github.com/drbergman/GlobalVsLocal/blob/29cabc8e270ec16b7d216409e1d5815280a4b49f/basePars_FGFR3.m#L11 |
| avoidWeekends.m | A function that makes sure any dosing events in the FGFR3 model occur on weekdays. |
| basePars_FGFR3.m | Creates a structure of base parameters for the FGFR3 model. |
| basePars_IL6.m | Creates a structure of base parameters for the IL6 model. |
| codensityFunction.m | Takes in locations of cells, a whole number _n_ for the nth nearest neighbors, a vector `types` identifying what type of agent is at each location, and outputs codensity calculations. |
| dosingRegimes.m | Creates the sequence of events for an FGFR3 simulation. Events include new doses of aFGFR3 and censoring. |
| eventProbabilities_IL6.m | Determines probabilities of all events in the next time step. |
| finishParameterComputation_IL6.m | Computes any input parameters that are functions of other parameters in the IL6 model. |
| fullGlobalODE_FGFR3.m | Computes the rate of change for the molecular dynamics in the FGFR3 model using the global method. |
| fullGlobalODE_IL6.m | Computes the rate of change for the molecular dynamics in the IL6 model using the global method. |
| initializeTumor_FGFR3.m | Initializes the tumor for the FGFR3 model. |
| initializeTumor_IL6.m | Initializes the tumor for the IL6 model. |
| iterativeEuler.m | Solves an ODE using direct Euler. If any of the state variables are negative, repeat the calculation with half the time step. |
| nextFileName.m | Helper function that outputs a new filename for saving data. |
| normalizeYLims.m | Helper function that normalizes _y_ limits for all axes in a figure. |
| saveSubstrateInfo.m | Saves substrate data in the IL6 model. |
| setupSolver_IL6.m | Sets up the solver in the IL6 model depending on the desired method. |
| simPatient_IL6.m | Simulates a patient in the IL6 model. |
| simTumor_FGFR3.m | Simulates a tumor until the next event in the FGFR3 model. Note: this gets called by startPatient_FGFR3.m each time a new event, i.e. new dose, occurs. |
| simTumor_IL6.m | Simulates a tumor in the IL6 model until the next event. Note: the only event in the IL6 model is censoring, so this always called once for each sample. |
| sphere_with_carveout.m | Helper function to create a sphere with the first octant carved out and color it with a heatmap. Used to create Figure 4CD. |
| startPatient_FGFR3.m | Simulates a patient in the FGFR3 model, calling simTumor_FGFR3.m between events. |
| substrateSolver_FGFR3.m | Solves the molecular dynamics in the FGFR3 model when using the local method. This includes both the PDE and the reaction ODEs. |
| substrateSolver_IL6.m | Solves the molecular dynamics in the IL6 model using either method. |
| updateTumor_FGFR3.m | Updates tumor cell agents in the FGFR3 model based on the events randomly chosen for them in the current update step. |
| updateTumor_IL6.m | Updates tumor cell agents in the IL6 model based on the events randomly chosen for them in the current update step. |
