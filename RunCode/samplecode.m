% Note that code includes first 50 models for storage space reasons. 
% This sample code calculates the growth of two metabolic models in their original environments and in a shared environment. It then calculates the growth of a possible endosymbiosis in the shared environment with each taking turns as host.

clear;
filedir='../Agora/';
%filedir='../CarveMe/';
%filedir='../KBase/';



% Gurobi settings
params = struct();
params.OutputFlag = 0;
params.FeasibilityTol=1e-9;


% compute growth rates for host and endo in own environs

ind1=10; % integer for a metabolic model for host
eval(['load ',filedir,'models_ehmat_format_first50/ehmodel',num2str(ind1),'.mat']);
nc=size(ehmodel.rhsc,1);
ehmodel1=ehmodel; % host endo symbiont model
clear ehmodel
[resulth,hostmodel]=runehmodel(ehmodel1); % computes growth rate of host model in its env


ind2=15; % integer for a metabolic model for endo
eval(['load ',filedir,'models_ehmat_format_first50/ehmodel',num2str(ind2),'.mat']);
ehmodel2=ehmodel;
clear ehmodel
[resulte,endomodel]=runehmodel(ehmodel2); % computes growth rate of endo model in its env

% create endo model
pairmodel=create_endo_pair_model(ehmodel1,ehmodel2,nc);
resultpairhe=gurobi(pairmodel,params);
% check if feasible
if ~strcmp(resultpairhe.status,'OPTIMAL')
    resultpairhe=0;
else
	resultpairhe=abs(resultpairhe.objval);
end



% grow host in shared environment
tempmodel1=hostmodel;
hostmodel.rhslb=pairmodel.rhs(1:nc);
hostmodel.rhsub=pairmodel.rhs(nc+1:2*nc);
resulth2env = runehmodel(hostmodel); % host metabolism in shared environment

% grow endo in shared environment
endomodel.rhslb=pairmodel.rhs(1:nc);
endomodel.rhsub=pairmodel.rhs(nc+1:2*nc);
resulte2env = runehmodel(endomodel);



% swap host and endo role in endosymbiosis
clear pairmodel
pairmodel=create_endo_pair_model(ehmodel2,ehmodel1,nc);
resultpaireh=gurobi(pairmodel,params);
if ~strcmp(resultpaireh.status,'OPTIMAL')
    resultpaireh=0;
else
	resultpaireh=abs(resultpaireh.objval);
end

% return data
['The growth of the host in its environment is ',num2str(abs(resulth.objval))]
['The growth of the endo in its environment is ',num2str(abs(resulte.objval))]
['The growth of the host in the shared environment is ',num2str(abs(resulth2env.objval))]
['The growth of the endo in the shared environment is ',num2str(abs(resulte2env.objval))]
['The growth of the endosymbiosis in the shared environment is ',num2str(resultpairhe)]
['The growth of the endosymbiosis swapped host in the shared environment is ',num2str(resultpaireh)]



