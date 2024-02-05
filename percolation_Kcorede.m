function [x,comp_sizes,BTR] = percolation_Kcorede(net)
%   Code for BTR calculation
%   Please cite the following paper if you use this function:
%   "Exploring the Link between Brain Topological Resilience 
%   and Cognitive Performance in the Context of 
%   Aging and Vascular Risk Factors: A Cross-Ethnicity Population-Based
%   Study"
    loc=[];newnet=net;
    Ksmatrix=kcoreness_centrality_bu(newnet)';
    for i=1:size(net,1)
        net=newnet;
        net(loc,:)=0;net(:,loc)=0;
        kcore=kcoreness_centrality_bu(net);
        str=strengths_und(net);
        mkcore=max(kcore);
        mstr=max(str(kcore==mkcore));
        strloc=find(str==mstr);
        loc=[loc;strloc(1)];
        net(loc,:)=0;net(:,loc)=0;
        [comps,comp_s]=get_components(net);
        Ksmatrix(:,i+1)=kcoreness_centrality_bu(net)';
        comp_sizes(i)=max(comp_s)/length(net);    
        x(i)=i/size(net,1);     
    end
    BTR=(sum(comp_sizes)-0.5*(comp_sizes(1)+comp_sizes(end)))/size(net,1)*100;
end

