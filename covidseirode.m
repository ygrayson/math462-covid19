%{
    This function will evaluate the SEIR model.
    
    This is the model we used for the ebola outbreak problem from lab 3:

    Sdot = -BSI                   Change in susceptible population
    Edot = BSI - AE               Change in susceptible population
    Idot = AE - rI                Change in susceptible population
    Rdot = rI                     Change in susceptible population
    ydot = NAE                    Change in susceptible population
    
%}
function dydt = covidseirode(t,x,params)

    B = params(1);
    A = params(2);
    r = params(3);
    N = params(4);
    
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    y = x(5);

    Sdot = -B.*S.*I;
    Edot = B.*S.*I - A.*E;
    Idot = A.*E - r.*I;
    Rdot = r.*I;
    ydot = N.*A.*E;
    dxdt = [Idot;Edot;Sdot;Rdot;ydot];
    
end