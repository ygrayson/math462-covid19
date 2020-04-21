%{
    This function will estimate the parameter values for the SEIR model
    using the least squares sum method.
%}
function LS = covidseirLS(v)
    dataCases = load('covid_data.txt');
    data1 = dataCases(:,2);
    B = v(1);
    A = v(2);
    r = v(3);
    N = 2000000;


    I0 = data1(1)./N;
    E0 = 2.*I0;
    S0 = 1 - I0 - E0;
    R0 = 0;
    y0 = N.*A.*E0;
    
    x0 = [S0,E0,I0,R0,y0];

    params = [B,A,r,N];
    tspan = dataCases(:,1);
    fun = @(t,x) covidseirode(t,x,params);
    [t,xsol] = ode45(fun,tspan,x0);
    zdot = xsol(:,5);
    ydot = data1;
    LS = sum((zdot - ydot).^2);
end
