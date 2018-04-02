function [sqdiff]=linklen(dx,dt,xnode,linkind,cmat,links,nfxnodes,slider,fxdangle,angs,crangvel,inodes,motornode,crtipnd,rotfix,rotangs)
    %outputs the square difference between current link lengths and
    %constrained link lengths
    %dx=[x1 x2 x3.... y1 y2 y3.....] %How to format dx| by breaking complex into x
    %and y components
    
    sqdiff=zeros(length(linkind),1); %initialize sqdiff
    xinter=xnode; %Copy to intermediate
    sizedx=length(nfxnodes)*2;
    
    %Allocation of objective function contributions from slider constraints
    sldof=0;
    for i=1:size(slider,1)
     if size(slider{i},2)==4
        sldof=sldof+1; % Middle node positon contribution to objective function 
     elseif size(slider{i},2)==3
        sldof=sldof+2; % Middle node position and length of links that make slider
     end   
    end
    sldiff=zeros(sldof,1);
    

    
    for i=1:length(nfxnodes) %loop through non-fixed nodes
        xinter(nfxnodes(i))=xinter(nfxnodes(i))+real(dx(i))+dx(sizedx/2+i)*1j; %Adds the perturbation dx to the copy
    end
    
    %Modifies dx for the slider joints according to sliding direction
    it=1;
    for i=1:size(slider,1)
        %If moving slider
        if length(slider{i})==3
            node1in=setdiff(cmat(slider{i}(2),:),slider{i}(1));
            node2in=setdiff(cmat(slider{i}(3),:),slider{i}(1));
            
            %Calculates the distance from the line (2 slider endpoints)
            %that the middle node is
            y2my1=imag(xinter(node2in)-xinter(node1in));
            x2mx1=real(xinter(node2in)-xinter(node1in));
            x2y1=real(xinter(node2in))*imag(xinter(node1in));
            y2x1=imag(xinter(node2in))*real(xinter(node1in));
            num=abs(y2my1*real(xinter(slider{i}(1)))-x2mx1*imag(xinter(slider{i}(1)))+x2y1-y2x1);
            denom=sqrt(y2my1^2+x2mx1^2);
            sldiff(it)=(num/denom); %Adds objective function contribution for distance middle slider node is from slider line
            it=it+1;
            
            %Adds objective function contribution to limit the length of
            %the two slider links to the original combined length
            sldiff(it)=(norm((xinter(node2in)-xinter(node1in)))-(links(slider{i}(2))+links(slider{i}(3))));
            it=it+1;
        end
        %If fixed slider
        if length(slider{i})==4
            y2my1=slider{i}(4);
            x2mx1=slider{i}(3);
            x2y1=(slider{i}(3)+real(inodes(slider{i}(1))))*imag(inodes(slider{i}(1),1));
            y2x1=(slider{i}(4)+imag(inodes(slider{i}(1))))*real(inodes(slider{i}(1),1));
            num=abs(y2my1*real(xinter(slider{i}(1)))-x2mx1*imag(xinter(slider{i}(1)))+x2y1-y2x1);
            denom=sqrt(y2my1^2+x2mx1^2);
            sldiff(it)=(num/denom); %Adds objective function contribution for distance middle slider node is from slider line
            it=it+1;
        end
    end
    
    %Creates objective function contributions for fixed angle constraints
    %Calculates the difference in angle at fixed angle locations between
    %adjustments
    if ~isempty(angs)
        angdiff=zeros(size(fxdangle,1),1);
        for i=1:size(fxdangle,1)

            ndind1=intersect(cmat(fxdangle(i,1),:),cmat(fxdangle(i,2),:));        
            ndind2=setdiff(cmat(fxdangle(i,1),:),ndind1);
            vec1=xinter(ndind2)-xinter(ndind1);
      
            ndind2=setdiff(cmat(fxdangle(i,2),:),ndind1);
            vec2=xinter(ndind2)-xinter(ndind1);         

            angnew=(real(vec1)*real(vec2)+imag(vec1)*imag(vec2))/(norm(vec1)*norm(vec2));            
            angdiff(i)=(angs(i)-(angnew));
        end
    end
    
    %Creates objective function contributions for rotation fixed constraints
    %Calculates the difference in angle at for each vector at a fixed rotation location between
    %adjustments   
     if ~isempty(rotangs)
        rotangdifftot=[];
        for i=1:length(rotfix)          
            [linkrow,nodecol]=find(cmat==rotfix(i));           
            for k=1:length(linkrow)
                node2in=cmat(linkrow(k),setdiff(1:2,nodecol(k)));
                node1in=rotfix(i);
                vec1=xinter(node2in)-xinter(node1in);
                vec2=rotangs{i}(k,1)+rotangs{i}(k,2)*1i;
                rotangint=1-dot(vec1,vec2)/(norm(vec1)*norm(vec2));               
                rotangdifftot=[rotangdifftot;rotangint];
            end
        end
    end
    
    for i=1:length(linkind) %loop over links that are being readjusted
        del=xinter(cmat(linkind(i),1))-xinter(cmat(linkind(i),2));
        sqdiff(i)=(abs(del)-links(linkind(i))); %computes the square difference
    end
    
    %Adds entry to objective function denoting how far sliding node is from
    %line defined by the slider link
    if ~isempty(slider)
        sqdiff=[sqdiff; sldiff];
    end
     %Adds entry to objective function for fixed angle constraints
    if ~isempty(fxdangle)
         sqdiff=[sqdiff; angdiff]; 
    end
     %Adds entry to objective function for rotation fix constraints
    if ~isempty(rotfix)
         sqdiff=[sqdiff; rotangdifftot]; 
    end
    
    %Adds objective function contribution for angular velocity of crank.
    %Calculates the angular velocity difference between adjustments to
    %match prescribed angular velocity
    newang=atan2d(imag(xinter(crtipnd)-xinter(motornode)),real(xinter(crtipnd)-xinter(motornode)));
    oldang=atan2d(imag(xnode(crtipnd)-xnode(motornode)),real(xnode(crtipnd)-xnode(motornode)));    
    vec1=xinter(crtipnd)-xinter(motornode);
    vec2=xnode(crtipnd)-xnode(motornode);
    perpdot=real(vec2)*imag(vec1)-imag(vec2)*real(vec1);
    newminold=asind(perpdot/(norm(vec1)*norm(vec2)));
    angveldiff=((newminold/dt)-crangvel);
    sqdiff=[sqdiff; angveldiff];
    %sqdiff=sqdiff.^2;
return