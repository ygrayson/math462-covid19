%{
    In this project, we plan to investigate the current situation of 
    COVID-19 using epidemic modeling. We will potentially be doing both 
    model selection and parameter estimation in order to find a good model 
    of the outbreak. One of the models that is worth trying is the SEIR 
    model that we discussed in class. We are looking at some 
    local regions in Michigan and how different regions might have a 
    different model as a good fit. With the existing data in hand, 
    we would hope to generate a good prediction into the future of the 
    outbreak. 

    We're gonna use the model we used for the ebola outbreak to help get
    started with the model.

        S is the number of susceptible people.
        E is the number of 
        I is the number of infected people.
        R is the number of recovered people.
        y is the number of case counts for infection.

%}
dataCases = load('covid_data.txt'); %we still need to get the data for this
days = dataCases(:,1);
case_count = dataCases(:,2);

%Plot the data.
figure
hold on
plot(days, case_count,'.');
xlabel('Time (days)');
ylabel('Total numbers of case counts');
title('Number of case counts as a function of time');

B = 0,9; %Beta
A = 0.9; %Alpha
r = 0.125; %Gamma
N = 500; %This is the total number of people

%Define the initial conditions.
I0 = case_count(1) / N;
E0 = 5 * I0;
S0 = 1 - I0 - E0;
R0 = 0;
y0 = N * A * E0;

%Define the vectors for the parameters and initial conditions.
params = [B,A,r,N];
x0 = [S0; E0; I0; R0; y0];
tspan = days;

%Call the function seirode to evaluate the values and then plot the
%ydot solutions vs time.
options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
fun = @(t,x) covidseirode(t,x,params);
[t,xsol] = ode45(fun, tspan, x0, options);
plot(tspan,xsol(:,5),'*');



