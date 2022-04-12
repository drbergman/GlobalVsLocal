# GlobalVsLocal
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
| agentODEJacobian_FGFR3.m | Takes in the state variables for an agent (or vector of agents) along with parameter values as well as the constant entries in the Jacobian of the ODE function and outputs the Jacobian. Used in `pde_solver_FGFR3.m` to update the reaction equations rather than using an RK method or direct Euler. |
| allCombos.m | Helper function that basically takes the _n_ outputs of `ndgrid` and makes each a column of a single array. Used to determine neighbors in lattice. See, for example, basePars_FGFR3.m Line 
