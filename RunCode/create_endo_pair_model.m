function endomodel=create_endo_pair_model(ehmodel1,ehmodel2,nc)

% unmapped is k=find(ecmap(:,1)==1 & ecmap(:,2)==-1);
% find biomass index
bmih=ehmodel1.bmi;
bmie=ehmodel2.bmi;


% construct host S matrices for E and C compartments


%%% construct endo model
%nc=size(unimat,1);
nrh=size(ehmodel1.Sh,2);
nre=size(ehmodel2.Se,2);
f=zeros(nrh+nre,1);
f(bmih)=-1; 
endomodel.obj = f;

% stochiometry matrix: external compartment then org1 then org2
b1=[1:nc];
b2=b1+nc;
b3=b2+nc;
endomat=sparse(3*nc+1,nrh+nre);


% check no biomass components affected: model2.A(e2ccomp(:,1)+2*nc,bmi2)
% e compartment lower bound
endomat(b1,1:nrh)=ehmodel1.Sh;
endomat(b1,nrh+1:nre+nrh)=ehmodel2.Se2h;
% e compartment higher bound
endomat(b2,1:nrh)=ehmodel1.Sh;
endomat(b2,nrh+1:nre+nrh)=ehmodel2.Se2h;
% c compartment endo
endomat(b3,nrh+1:nre+nrh)=ehmodel2.Se;
% fix biomass rxn
%endomat(:,bmih)=endomat(:,bmih)+endomat(:,bmie+nrh);
%endomat(:,bmie+nrh)=0;
endomat(end,bmih)=1;
endomat(end,bmie+nrh)=-1; % ensures they grow at same rate

endomodel.A=sparse(endomat);


% create rhs that will be modified
endomodel.rhs=zeros(3*nc+1,1);
% move the rhs for compounds moved from environment to inside host

endomodel.rhs(b1)=ehmodel1.rhslb+ehmodel2.rhslb;
endomodel.rhs(b2)=ehmodel1.rhsub+ehmodel2.rhsub;
% endomodel.rhs(b3)=zeros(nc,1);
% endmodel.rhs(end)=0;
endomodel.sense = [repmat('>',nc,1);repmat('<',nc,1);repmat('=',nc+1,1)];
endomodel.lb=[ehmodel1.lb;ehmodel2.lb];
endomodel.ub=[ehmodel1.ub;ehmodel2.ub];

% some models have an upper limit on biomass
endomodel.lb(bmih)=0;
endomodel.ub(bmih)=1000;
endomodel.lb(bmie+nrh)=0;
endomodel.ub(bmie+nrh)=1000;
end
