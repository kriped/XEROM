coefficient = 0.77;
feedback_data = load("..\..\input\R4C38\EOC\Refinement1\FEEDBACK_data.mat");
%Rescale feedback coefficient. Variable n   aming must stay consistent.
K11 = feedback_data.K11*coefficient;
K12 = feedback_data.K12*coefficient;
K21 = feedback_data.K21*coefficient;
K22 = feedback_data.K22*coefficient;

save("..\..\input\R4C38\UNSTABLE\Refinement1\FEEDBACK_data.mat","K11","K12","K21","K22")

% Optionally, display a message indicating successful saving of the data
disp('Scaled feedback data has been saved successfully.');


