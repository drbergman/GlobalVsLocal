function input = simPatient_IL6(input)

%% initialize inputs

input = finishParameterComputation_IL6(input);
if input.flags.plotFigs
    input.plot_properties = initializeFigure_IL6(input.plot_properties);
end
input = simTumor_IL6(input);

