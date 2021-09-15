function D = DiffusionCoefficient(D0, H, T)


D = D0*exp(-H*1000/8.314./T);



%T=Inf, D=D0; T=1000; D=D1000, H=R*ln(D0/D1000)
% D = (7e-8)*exp(-273000/8.314./T); %Quartz DTi, Cherniak et al., 2007, Chem Geol
% D = 10.^(-8.3-311000/8.314/2.303./T); %Quartz DTi, Jollands et al., 2020, Geology
% D = (2.48e-11)*exp(-199000/8.314./T); %Quartz DAl, Tailby et al., 2018, American Mineralogist

% D = (9.6e-7)*exp(-278000/8.314./T); %Zircon DLi, Trail et al., 2016, Contrib Mineral Petrol

% D = 10.^(-10.06)*exp(-229000/8.314./T); %Olivine DP, Watson et al., 2015, American Mineralogist
% D = 10.^(-5.92-12847./T); %Olivine DLi, "slow" diffusion, Dohmen et al., 2010, Geochimica et Cosmochimica Acta

% D = 0.955e-4*exp(-406000/8.314./T); %Diopside DFe-Mg, Dimanov and Sautter, 2000, Eur. J. Mineral.

% D = 0.29*exp(-455000/8.314./T); %Sanidine DBa, กอ001, Cherniak 2002, Geochimica et Cosmochimica Acta, sensitive to C boundary conditions
% D = 8.4*exp(-450000/8.314./T); %Sanidine DSr, กอ001, Cherniak 1996, Geochimica et Cosmochimica Acta

% XAn=0.83;
% D = 10^(-5.2-3.3*XAn)*exp(-264000/8.314./T); % Plagioclase DSr, McCarthy et al., 2020, Earth and Planetary Science Letters
%compile the Sr diffusion in plagioclase measurements of Cherniak and Watson (1994) and Giletti and Casserly (1994)

% P=6000; %bar
% H=179894+0.6*(P-1);
% D = (6.63e-14)*exp(-H/8.314./T); %Garnet DCa-(Fe,Mg), Vielzeuf et al., 2007, Contrib Mineral Petrol



