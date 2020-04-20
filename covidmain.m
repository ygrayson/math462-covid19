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
        E is the number of exposed people.
        I is the number of infected people.
        R is the number of recovered people.
        y is the number of case counts for infection.

%}
dataCases = load('covid_data.txt'); %we still need to get the data for this

%Plot the data.
figure
hold on
plot(dataCases(:,1),dataCases(:,2),'.');
xlabel('Time (hours)');
ylabel('Total numbers of case counts');
title('Number of case counts as a function of time');

B = .7828; %Beta
A = .3816; %Alpha 
r = .7602; %Gamma
N = 700000; %This is the total number of people 

%Define the initial conditions.
I0 = dataCases(1,2)./N;
E0 = 2.*I0;
S0 = 1 - I0 - E0;
R0 = 0;
y0 = N.*A.*E0;

%Define the vectors for the parameters and initial conditions.
params = [B,A,r,N];
x0 = [S0,E0,I0,R0,y0];
tspan = dataCases(:,1);

%Call the function seirode to evaluate the values and then plot the
%ydot solutions vs time.
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
fun = @(t,x) covidseirode(t,x,params);
[t,xsol] = ode45(fun,tspan,x0,options);
plot(tspan,xsol(:,5),'*');

%use fminsearch (nelder meed)
%seir_cost :
% inputs: alpha, beta and gamma
% outputs is -LL = {sum(y(ti))-sum(zi*ln(y(ti)))}

%This portion solves for the Poisson LL to help determine better 
%parameters for the SEIR model to better fit the data.

fun1 = @(v) covidseirPois(v);
v0 = [B,A,r];
%fminsearch helps find the minimal value for the parameters. Which in turns
%optimizes the ODE.
varpars = fminsearch(fun1,v0);
%{
B = varpars(1);
A = varpars(2);
r = varpars(3);
params = [B,A,r];

tspan = dataCases(:,1);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
fun = @(t,x) covidseirode(t,x,params);
[t,xsol] = ode45(fun,tspan,x0,options);
plot(tspan,xsol(:,5),'*');
%}

