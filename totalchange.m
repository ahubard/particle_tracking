function file_ftoi = totalchange(folder,En)
%folder = 1;
%En = 103;
file_load =sprintf('/aline%i/rotdrum%i/o%02d/Avalanches_%i.mat',folder,folder,En,En);

load(file_load,'Number_Avalanches','Noavalanches','Avalanche_time');
deltaR = [];%Difference between positions at begging and end of the avalanche. 
NoParticles_moved = [];
Max_particle_dis = []; %Maximal displacement of a single particle. 
Initial_Angle = [];
Final_Angle = [];

for dn = 1:length(Noavalanches)
    T1 = Avalanche_time{1,dn};
    T2 = Avalanche_time{2,dn};
    fnn =sprintf('/aline%i/rotdrum%i/o%02d/Displacement_%i.mat',folder,folder,En,dn);
    load(fnn,'PX','PY','diskmove');
    if (diskmove)
        dr = sum(sqrt((PX(diskmove,T2)-PX(diskmove,T1)).^2+(PY(diskmove,T2)-PY(diskmove,T1)).^2));
        deltaR = [deltaR dr];
        dp = sum(((PX(diskmove,T2)-PX(diskmove,T1)).^2+(PY(diskmove,T2)-PY(diskmove,T1)).^2)>sqrt(2));
        NoParticles_moved = [NoParticles_moved dp];
        maxp = max(sqrt((PX(diskmove,T2)-PX(diskmove,T1)).^2+(PY(diskmove,T2)-PY(diskmove,T1)).^2));
        Max_particle_dis = [Max_particle_dis maxp];
        Initial_Angle = [Initial_Angle estimate_angle(PX(:,T1),PY(:,T1))];
        Final_Angle = [Final_Angle estimate_angle(PX(:,T2),PY(:,T2))];
        
    end
    
end
file_ftoi = sprintf('/aline%i/rotdrum%i/o%02d/Initial_to_Final_%i.mat',folder,folder,En,En);
save(file_ftoi,'deltaR','NoParticles_moved','Max_particle_dis','Initial_Angle','Final_Angle');