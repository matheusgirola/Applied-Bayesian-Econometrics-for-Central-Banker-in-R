function States=simS(sin,pmat)
States=sin;
Rnd=rand(size(sin,1),1);
for i=2:size(sin,1)

state_past=find(States(i-1,:)==1);

    if Rnd(i,1)<pmat(state_past,state_past)

        States(i,state_past)=1; % when staying at last state

    else    % when changing to other states

        idx_other=find(States(i-1,:)==0);
        Prob2=pmat(:,state_past);

        a=[pmat(state_past,state_past) ; Prob2(idx_other)];

        cum_sum=cumsum(a);
        sorted=sort([cum_sum ; Rnd(i,1)]); % throw the prob at cumsum of other states to get
        % where it stands (where to
        % switch)

        idx=find(Rnd(i,1)==sorted)-1;      % find index

        States(i,idx_other(idx))=1;        % change state

    end  
end