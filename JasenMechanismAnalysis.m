%|Numerical PVA Code|University of Illinois at Urbana-Champaign| 
%|ME 370|Created by Brian C. McGuigan
close all
clear
clc

% Section1:Input
x1 = 0;
x2 = 20i;
x3 = -70;
x4 = -72.5904+99.9664i;
x5 = -159.419+10.2127i;
x6 = -116.129-90.9111i;
x7 = -58.7787-50.7748i;
x8 = -38.6321-128.196i;
inodes = [x1; x2; x3; x4; x5; x6; x7; x8];			% Node/Joint Placement
cmat = [1 2; 1 3; 2 4; 4 3; 4 5; 3 5; 5 6; 3 7; 2 7; 6 7; 6 8; 7 8];		% Connectivity Matrix
gndlink = 2;						% Link Index of Ground Link (row of cmat denoting nodes to
									% freeze during simulation)
crlink = 1;							% Link with prescribed position over time (input crank)
motornode = 1;						% Node the motor is at
slider = {};						% Slider constraints
fxdangle = [4 6;10 12];				% Fix link/link angle during rotation. (Must share a node)
rotfix = [];						% Nodes that are fixed so that they cant rotate
slfric = [];						% Friction at slider constraints

tperiod = 30;						% total time to run simulation (seconds)
dt = 0.1;							% time step (seconds)
tnum = tperiod/dt;					%Total number of timesteps
time = 0:dt:tperiod;				%array holding the amount of time that has passed
crangvel = 30;						% Angular velocity of input crank (deg/sec) -Assumed constant.
									%+(-) angular velocity denotes CCW (CW) rotation

% Section2: Useful information about the sytem
%fxnodes = unique(setdiff(cmat([crlink,gndlink],:),motornode));
% Nodes with prescribed/fixed placement (the ones that make up the ground and crank link)
fxnodes = unique(cmat([gndlink],:));
 if ~isempty(slider)
     for j = 1:length(slider)
        fxnodes = setdiff(fxnodes,slider{j}(1));
     end
end
 
nfxnodes = setdiff(1:length(inodes),fxnodes);	% Nodes without prescribed/fixed placement
												% (the ones that can vary after we move crank link a bit)
%crpvtnd = intersect(cmat(gndlink,:),cmat(crlink,:));%Node index that input crank link rotates about
crpvtnd = motornode;							%Assigns previous variable entry
crtipnd = setdiff(cmat(crlink,:),crpvtnd);		%Node index denoting the tip of the input crank link
ndof = length(nfxnodes)*2;
%Nu(3*i,end)=1;									%Torque coeefficient on crank link
linkind = setdiff(1:length(cmat),gndlink);		%links that may have varying lengths after a timestep

% Slider information used for DFA
sliderlinks=[];
sliderendnodes=[];
if ~isempty(slider)
     for j=1:length(slider)
         linkind=setdiff(linkind,slider{j}(2));
        if length(slider{j})==3 
            linkind=setdiff(linkind,slider{j}(3));
            sliderlinks=[sliderlinks;j, slider{j}(2), slider{j}(3);];%total list of all slider links [sliderindex, link1ind linke1ind]
            cmnode=intersect(cmat(slider{j}(2),:),cmat(slider{j}(3),:));
            endnodes=setdiff([cmat(slider{1}(2),:) cmat(slider{1}(3),:)],cmnode);
            sliderendnodes=[sliderendnodes;endnodes(1),endnodes(2);];
        end
     end
end

% Section3: Allocation/Initialization
xnode=zeros(length(inodes),tnum); %allocate node positions over all time steps 
xnode(:,1)=inodes; %Initialize current timestep (put in initial node position data)

% Section4: Put in initial fixed ground link nodes positions at each time step in corresponding position array since they are fixed
for i=2:tnum
    for j=1:length(gndlink)
        xnode(cmat(gndlink(j),:),i)=inodes(cmat(gndlink(j),:));
    end
end


% Section6: Computes link lengths from initial node positions
links=zeros(length(cmat),1); %Allocate array for link lengths
for i=1:length(cmat) %loop over connectivity matrix
	currLinkNode1 = int8(cmat(i, 1));
	currLinkNode2 = int8(cmat(i, 2));
    links(i) = sqrt((real(inodes(currLinkNode2))-real(inodes(currLinkNode1)))^2 ...
	+(imag(inodes(currLinkNode2))-imag(inodes(currLinkNode1)))^2);
end

%Compute initial link angle for fixed angle constraint
angs=zeros(size(fxdangle,1),1);
angshnode=zeros(size(fxdangle,1),1);
rotangs=cell(size(rotfix,1),1);

%Calculates the initial angles for the fixed angle constraints
if ~isempty(angs)
    for i=1:size(fxdangle,1)
            
        ndind1=intersect(cmat(fxdangle(i,1),:),cmat(fxdangle(i,2),:));
        ndind2=setdiff(cmat(fxdangle(i,1),:),ndind1);
        vec1=inodes(ndind2)-inodes(ndind1);
      
        ndind2=setdiff(cmat(fxdangle(i,2),:),ndind1);
        vec2=inodes(ndind2)-inodes(ndind1);      
        angs(i)=(real(vec1)*real(vec2)+imag(vec1)*imag(vec2))/(norm(vec1)*norm(vec2));
    end
end

%Calculates the initial link vectors from each rotation fixed constraint location
if ~isempty(rotangs)
    for i=1:length(rotfix)        
        [linkrow,nodecol]=find(cmat==rotfix(i));
        rotangs{i}=zeros(length(linkrow),2);
        for k=1:length(linkrow)
            node2in=cmat(linkrow(k),setdiff(1:2,nodecol(k)));
            node1in=rotfix(i);
            vec=inodes(node2in)-inodes(node1in);
            
            rotangs{i}(k,1)=real(vec);
            rotangs{i}(k,2)=imag(vec);
        end
        
    end
end

% Section5: Assign crank tip position at all times in simulation.
crank = cmat(int8(crlink));
if cmat(int8(crlink), 1) == motornode
	crTipIdx = cmat(int8(crlink), 2);
else
	crTipIdx = cmat(int8(crlink), 1);
end
crankAngle = zeros(tnum, 1);
origAngle = atand((imag(inodes(crTipIdx))-imag(inodes(motornode)))/(real(inodes(crTipIdx)) - real(inodes(motornode))));
crankAngle(1, 1) = origAngle;
for i=2:tnum
	currAngle = origAngle + crangvel*dt*(i-1);
    crankAngle(i, 1) = currAngle;
%	if currAngle <= 0
%		currAngle = currAngle+360;
%	elseif currAngle >= 360
%		currAngle = currAngle-360;
%	end
%	currNodeReal = real(inodes(motornode))+cosd(currAngle)*links(crTipIdx);
%	currNodeImag = imag(inodes(motornode))+sind(currAngle)*links(crTipIdx);
%	xnode(crTipIdx,i) = currNodeReal + currNodeImag*i;
end


%Section7: Loop over each time step and find adjustments that preserve link lengths
str=sprintf('Running Simulation:'); %Make string
disp(str) %Display String
linksnew=zeros(length(linkind),1); %Allocate array for new link lengths (link lengths that could change in objective function after dx is solved
for i=1:tnum-1
    
    xnode(nfxnodes,i+1)=xnode(nfxnodes,i);%copy previous nfxnode quantities to current timestep
    options = optimset('Display', 'off'); %Turns off display options
    %options.OptimalityTolerance=1e-12;
    %options.FunctionTolerance=1e-12;
    %options.StepTolerance=1e-12;
    %options.Algorithm = 'levenberg-marquardt';
    func=@(dx)linklen(dx,dt,xnode(:,i+1),linkind,cmat,links(:),nfxnodes,slider,fxdangle,angs,crangvel,inodes,motornode,crtipnd,rotfix,rotangs); %Sets function handle to include additional arguments
    dxnfx=lsqnonlin(func,1e-6*ones(1,ndof),[],[],options); %Solve for nfxnode perturbation distance that fixes link length incompatibilities
    xnode(nfxnodes,i+1)=xnode(nfxnodes,i+1)+dxnfx(1:(ndof/2))'; %add real part correction
    xnode(nfxnodes,i+1)=xnode(nfxnodes,i+1)+dxnfx((ndof/2+1):end)'*1i; % add imag part correction
    %motornodevel=xnode(motornode,i+1)-xnode(motornode,i);
    
    %Loops over each link that can have varied length to compute new link lengths (They should be constant if solved correctly)
    for k=1:length(linkind) %loop over linkind
        del=xnode(cmat(linkind(k),1),i+1)-xnode(cmat(linkind(k),2),i+1);
        linksnew(k)=sqrt(real(del)^2+imag(del)^2); %compute length between nodes
    end
    linksnewtot=sum(linksnew);%Sums new link lengths (this should equal sum(links) during whole simulation)
    
    percdone=(i/tnum)*100.0; % Percentage of simulation calculated
    str=sprintf('Progress=%.2f %%: i=%i :Total Link Lengths=%.2f,',percdone,i,linksnewtot); %Make string with progress info
    disp(str)
end


%Section8: Animate Linkage
figsize=5.0; %relative figure size in inches 
xmin=min(min(real(xnode)));
xmax=max(max(real(xnode)));
ymin=min(min(imag(xnode)));
ymax=max(max(imag(xnode)));
hfig=figure(1); %initiate figure
figscale=figsize/max([xmax-xmin ymax-ymin]);
set(hfig,'units','inches','position',[2 2 (xmax-xmin)*figscale (ymax-ymin)*figscale]);
hax=axis(); %initiate axis
hold on
xlim([xmin xmax]) %Scale xaxis to fit max/min x of simulation
ylim([ymin ymax]) %Scale yaxis to fit max/min y of simulation
for k=1:tnum %loop over all timesteps
    hax = gca; %initialize current axis
    hax.ColorOrderIndex = 1; %reset color index iteration (keeps colors the same between time steps)
    for i=1:length(cmat) %loop over links
            p1in=cmat(i,1); %link i node 1
            p2in=cmat(i,2); %link i node 2
            scatter(real([xnode(p1in,1:k),xnode(p2in,1:k)]),imag([xnode(p1in,1:k),xnode(p2in,1:k)]),'.k'); %Plot position of each node as black dot 
            %If groundlink, plot dashed line, else make solid
            if sum(i==gndlink)==1
                plot(real([xnode(p1in,k),xnode(p2in,k)]),imag([xnode(p1in,k),xnode(p2in,k)]),'--') %Plot dashed line between
            else
                plot(real([xnode(p1in,k),xnode(p2in,k)]),imag([xnode(p1in,k),xnode(p2in,k)])) %Plot line between
            end
            midpos=xnode(p1in,k)+0.5*(xnode(p2in,k)-xnode(p1in,k));
            text(real(midpos),imag(midpos),[ '$\raisebox{.5pt}{\bf{\raisebox{-.9pt} {' num2str(i) '}}}$'],'Interpreter','latex');
	end
    %Puts node labels at nodes
    for j=1:length(inodes)
        text(real(xnode(j,k)),imag(xnode(j,k)),[ '$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {' num2str(j) '}}}$'], 'Interpreter', 'latex')
    end
    
    hxlabel=xlabel('X position');
    hylabel=ylabel('Y Position');
    htitle=title('Mechanism Animation');
    formatplot(hfig,hax,hxlabel,hylabel,htitle)
    pause(0.05) %how long to pause plot in seconds
    if k~=tnum %If it reaches the end, it wont clear the plot
        cla %otherwise clear axis to plot next time step
    end
end

%Section9: Compute Velocity and Acceleration of nodes
%Velocity Calculation
velocity=zeros(length(inodes),tnum-1);
for i=1:tnum-1
	if i==1
		velocity(:,i)=(xnode(:,i+1)-xnode(:,i))/dt;
	else
		velocity(:,i)=(xnode(:,i+1)-xnode(:,i-1))/(2*dt);
	end
end

%Acceleration code goes here
acceleration=zeros(length(inodes),tnum-2);
for i=1:tnum-2
    %First order differentiation
    if i==1
        acceleration(:,i)=(velocity(:,i+1)-velocity(:,i))/dt;
    end
    %Use symmetric difference differentiation (less error)
    if i~=1
       acceleration(:,i)=(velocity(:,i+1)-velocity(:,i-1))/(2*dt);
    end
end

%Additional smoothening of the position/velocity data may need to be
%implemented if there are "spikes" in the acceleration. This has to do with
%the stability of the numerical procedure used?

xPosition = zeros(tnum, 1);
yPosition = zeros(tnum, 1);
xVelocity = zeros(tnum, 1);
yVelocity = zeros(tnum, 1);
xAcceleration = zeros(tnum, 1);
yAcceleeration = zeros(tnum, 1);
xVelocity(1, 1) = real(velocity(8, 1));
yVelocity(1, 1) = imag(velocity(8, 1));
xAcceleration(1, 1) = real(acceleration(8, 1));
yAcceleration(1, 1) = imag(acceleration(8, 1));
xAcceleration(2, 1) = real(acceleration(8, 1));
yAcceleration(2, 1) = imag(acceleration(8, 1));
for i=1:tnum
    xPosition(i, 1) = real(xnode(8, i));
    yPosition(i, 1) = imag(xnode(8, i));
    if i ~= 1
        xVelocity(i, 1) = real(velocity(8, i-1));
        yVelocity(i, 1) = imag(velocity(8, i-1));
        if i ~= 2
            xAcceleration(i, 1) = real(acceleration(8, i-2));
            yAcceleration(i, 1) = imag(acceleration(8, i-2));
        end
    end
    
end
f1 = figure;
plot(crankAngle, xPosition);
title('Position v.s. crank angle (x coordinates v.s. crank angle)');
xlabel('Crank angle, degree');
ylabel('X coordinate of the foot, node 8');
f2 = figure;
plot(crankAngle, yPosition);
title('Position v.s. crank angle (y coordinates v.s. crank angle)');
xlabel('Crank angle, degree');
ylabel('Y coordinate of the foot, node 8');
f3 = figure;
plot(crankAngle, xVelocity);
title('Velocity v.s. crank angle (velocity in x direction v.s. crank angle)');
xlabel('Crank angle, degree');
ylabel('velocity of node 8 in x direction');
f4 = figure;
plot(crankAngle, yVelocity);
title('Velocity v.s. crank angle (velocity in y direction v.s. crank angle)');
xlabel('Crank angle, degree');
ylabel('velocity of node 8 in y direction');
f5 = figure;
plot(crankAngle, xAcceleration);
title('Acceleration v.s. crank angle (acceleration in x direction v.s. crank angle)');
xlabel('Crank angle, degree');
ylabel('acceleration of node 8 in x direction');
f6 = figure;
plot(crankAngle, yAcceleration);
title('Acceleration v.s. crank angle (acceleration in y direction v.s. crank angle)');
xlabel('Crank angle, degree');
ylabel('acceleration of node 8 in y direction');