%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_results.m         %
% last modified:  18,2021                                  %
% author: Yixiao Zhang - yzhang34@wpi.edu                  %
% plot results and b.c. for trusses,beams, heat conduction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_results_2021(type,xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,axial,fnn,Qe,yield,g,E)
[dum1,dum2,frequence]=size(Ucomp); 
[dum1,dum2,frequence1]=size(yield);
[dum1,frequence2]=size(ien);
[dum1,frequence3]=size(nel);
[dum1,frequence4]=size(axial);
[dum1,frequence5]=size(Qe);
[dum1,frequence6]=size(E);
[dum1,dum2,frequence7]=size(xn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       TRUSS AND BEAM                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(type, 'truss') || strcmp(type,'beam'))
    
    display_factor=.08;
    eps=1e-4;
    
    % characteristic distances
    xmax=max(max(xn(1,:,1)));
    xmin=min(min(xn(1,:,1)));
    
    if (nsd >1)
        ymax=max(max(xn(2,:,1)));
        ymin=min(min(xn(2,:,1)));
        Lcar=max([xmax-xmin;ymax-ymin]);
        if (nsd >2)
            zmax=max(max(xn(3,:,1)));
            zmin=min(min(xn(3,:,1)));
            Lcar=max([xmax-xmin;ymax-ymin;zmax-zmin]);
        end
    else
        Lcar=xmax-xmin;
    end
    k=0;
    while (k~=6)
        k=menu('Results','Undeformed mesh and BCs', 'Undeformed and deformed mesh','Fracture mesh','Fracture video','Fracture Model','Exit');
        
        switch k
            case 1
                figure;
%                 figure('units','normalized','outerposition',[0 0 1 1])      
                axis equal;
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                title('Undeformed mesh and BCs');
                hold on;
                set(gcf, 'Color', [1,1,1]); %Background color white
                plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                 plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
%                   numbers(nel(1),ien{1},xn(:,:,1),nnp,nsd);
%                 numbers(nel(2),ien{2},xn,nnp,nsd);
                plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd);
                 plot_bc_predisplacement(type,Lcar,g,display_factor,nnp,xn,nsd);
                plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd);
%                
%           view(141,50);
%           view(-217.5,30);
              
                hold off;
            case 2
                figure;
                axis equal;
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                title('Undeformed and deformed mesh');
                hold on;
                set(gcf, 'Color', [1,1,1]); %Background color white
                   plot_mesh_undeformed1(nel(1),ien{1},xn(:,:,1),nnp,nsd,E{1});
%                 plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
%                 numbers(nel,ien,xn,nnp,nsd);Ucomp(:,:,2)
                [xnn]= plot_mesh_deformed(type,xn,Ucomp(:,:,frequence),Idb,display_factor,Lcar,nel(frequence3),...
                                    ien{frequence2},ndf,nsd,nen,axial{frequence4},Qe{frequence5},E{frequence6});
%                plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,axial,Qe);
               plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                view(139,76);
%                view(108.5,15.5);
%                  view(-217.5,30);
                hold off;
              case 3
                figure;
                axis equal;
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                title('Fracture mesh');
                hold on;
                set(gcf, 'Color', [1,1,1]); %Background color white
%                 plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
%                 numbers(nel,ien,xn,nnp,nsd);
%                 plot_mesh_undeformed(nel(2),ien{2},xn,nnp,nsd,E{2});
                [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence),Idb,display_factor,Lcar,nel(frequence3),...
                                   ien{frequence2},ndf,nsd,nen,axial{frequence4},Qe{frequence5},E{frequence6});
%                 plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-3),Idb,display_factor,Lcar,nel(frequence3-3),...
%                                    ien{frequence2-3},ndf,nsd,nen,axial{frequence4-4},Qe{frequence5-4},E{frequence6-3});               
                plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
            
             view(139,76);
                hold off;   
                
                case 4
                 figure('units','normalized','outerposition',[0 0 1 1])
                axis equal;
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                title('Fracture History');
                hold on;
                
                for i = 1:frequence
                    clf('reset')
                    axis equal;
%                   plot_mesh_undeformed1(nel(1),ien{1},xn,nnp,nsd,E{1})
                    [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,i),Idb,display_factor,Lcar,nel(i),...
                               ien{i},ndf,nsd,nen,axial{i},Qe{i},E{i});
                    plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
%                        view(108.5,15.5);
%                       view(-217.5,30);
                     view(139,76);
                    pause(0.1)
                    
                end
               
                hold off  
                
                
            case 5
                  k1=menu('Mode Shape Range','1','2');
                    
                  switch k1
                      case 1
                    figure;
                        subplot(2,2,1)
                        axis equal;
                        title('Fracture Model 1');
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                         [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,round(frequence*5/8)),Idb,display_factor,Lcar,nel(round(frequence3*5/8)),...
%                             ien{round(frequence2*5/8)},ndf,nsd,nen,axial{round(frequence4*5/8)},Qe{round(frequence5*5/8)},E{round(frequence6*5/8)});
                       [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-3),Idb,display_factor,Lcar,nel(frequence3-3),...
                            ien{frequence2-3},ndf,nsd,nen,axial{frequence4-3},Qe{frequence5-3},E{frequence6-3});
                          
                        
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        subplot(2,2,2)
                        axis equal;
                        title('Fracture Model 2');
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                         [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,round(frequence*6/8)),Idb,display_factor,Lcar,nel(round(frequence3*6/8)),...
%                             ien{round(frequence2*6/8)},ndf,nsd,nen,axial{round(frequence4*6/8)},Qe{round(frequence5*6/8)},E{round(frequence6*6/8)});
%                          
                        
                        [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-2),Idb,display_factor,Lcar,nel(frequence3-2),...
                            ien{frequence2-2},ndf,nsd,nen,axial{frequence4-2},Qe{frequence5-2},E{frequence6-2});
                        
                        
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        subplot(2,2,3)
                        axis equal;
                        title('Fracture Model 3')
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                         [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,round(frequence*7/8)),Idb,display_factor,Lcar,nel(round(frequence3*7/8)),...
%                             ien{round(frequence2*7/8)},ndf,nsd,nen,axial{round(frequence4*7/8)},Qe{round(frequence5*7/8)},E{round(frequence6*7/8)});
                        [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-1),Idb,display_factor,Lcar,nel(frequence3-1),...
                            ien{frequence2-1},ndf,nsd,nen,axial{frequence4-1},Qe{frequence5-1},E{frequence6-1});
                        
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        subplot(2,2,4)
                        axis equal;
                        title('Fracture Model 4')
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
                        [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence),Idb,display_factor,Lcar,nel(frequence3),...
                            ien{frequence2},ndf,nsd,nen,axial{frequence4},Qe{frequence5},E{frequence6});
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        
                      case 2
                            subplot(2,2,1)
                        axis equal;
                        title('Fracture Model 1');
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                         [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,round(frequence*5/8)),Idb,display_factor,Lcar,nel(round(frequence3*5/8)),...
%                             ien{round(frequence2*5/8)},ndf,nsd,nen,axial{round(frequence4*5/8)},Qe{round(frequence5*5/8)},E{round(frequence6*5/8)});
                       [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-12),Idb,display_factor,Lcar,nel(frequence3-12),...
                            ien{frequence2-12},ndf,nsd,nen,axial{frequence4-12},Qe{frequence5-12},E{frequence6-12});
                          
                        
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        subplot(2,2,2)
                        axis equal;
                        title('Fracture Model 2');
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                         [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,round(frequence*6/8)),Idb,display_factor,Lcar,nel(round(frequence3*6/8)),...
%                             ien{round(frequence2*6/8)},ndf,nsd,nen,axial{round(frequence4*6/8)},Qe{round(frequence5*6/8)},E{round(frequence6*6/8)});
%                          
                        
                        [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-6),Idb,display_factor,Lcar,nel(frequence3-6),...
                            ien{frequence2-6},ndf,nsd,nen,axial{frequence4-6},Qe{frequence5-6},E{frequence6-6});
                        
                        
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        subplot(2,2,3)
                        axis equal;
                        title('Fracture Model 3')
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
%                         [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,round(frequence*7/8)),Idb,display_factor,Lcar,nel(round(frequence3*7/8)),...
%                             ien{round(frequence2*7/8)},ndf,nsd,nen,axial{round(frequence4*7/8)},Qe{round(frequence5*7/8)},E{round(frequence6*7/8)});
                        [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence-3),Idb,display_factor,Lcar,nel(frequence3-3),...
                            ien{frequence2-3},ndf,nsd,nen,axial{frequence4-3},Qe{frequence5-3},E{frequence6-3});
                        
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                        subplot(2,2,4)
                        axis equal;
                        title('Fracture Model 4')
                        hold on
                        %                         plot_mesh_undeformed(nel(1),ien{1},xn,nnp,nsd,E{1});
                        [xnn]=plot_mesh_deformed(type,xn,Ucomp(:,:,frequence),Idb,display_factor,Lcar,nel(frequence3),...
                            ien{frequence2},ndf,nsd,nen,axial{frequence4},Qe{frequence5},E{frequence6});
                        plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xnn,nsd);
                        view(139,76);
                          
                  end
                        
                                %             case 6
%                 figure;
%                 subplot(2,1,1)
%                 %axis equal;
%                 title('Shear Diagram');
%                 hold on;
%                 set(gcf, 'Color', [1,1,1]); %Background color
%                 plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
% %                 numbers(nel,ien,xn,nnp,nsd)
%                 shear_diagram(xn,ien,nen,axial,nsd,nel,nnp)
%                 subplot(2,1,2)
%                 %axis equal;
%                 title('Moment Diagram')
%                 hold on
%                 set(gcf, 'Color', [1,1,1]); %Background color
%                 plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
%                 numbers(nel,ien,xn,nnp,nsd)
%                 moment_diagram(xn,ien,nen,axial,nsd,nel,nnp,fnn)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         HEAT                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(type, 'heat'))
    k=0;
    
    while (k~=4)
        k=menu('results','mesh','Temperature with numbers','Temperature w\o numbers','Fluxes','exit');
        
        switch k
            case 1
                figure;
                axis equal;
                title('Mesh');
                hold on;
                if (nen == 3)
                    triplot(ien',xn(1,:),xn(2,:));
                    numbers_tria(nel,ien,xn,nnp,nsd);
                else
                    for e=1:nel
                        for m=1:nen
                            node(m)=ien(m,e);
                        end
                        x=[xn(1,node(1)); xn(1,node(2)); xn(1,node(3)); xn(1,node(4)); xn(1,node(1))];
                        y=[xn(2,node(1)); xn(2,node(2)); xn(2,node(3)); xn(2,node(4)); xn(2,node(1))];
                        plot(x,y,'r-o');
                    end
                    numbers_quad(nel,ien,xn,nnp,nsd);
                end
                hold off;
            case 2
                figure;
                axis equal;
                %           title('Temperature Distribution');
                hold on;
                triplot(ien',xn(1,:),xn(2,:));
                plot_temp_tria(Ucomp,nel,ien,xn,nen,nnp,nsd,1);
                hold off;
            case 3
                figure;
                axis equal;
                title('Temperature Distribution');
                hold on;
                %            triplot(ien',xn(1,:),xn(2,:));
                plot_temp_tria(Ucomp,nel,ien,xn,nen,nnp,nsd,0);
                hold off;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       PLOT FUNCTIONS                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undeformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_undeformed(nel,ien,xn,nnp,nsd,E)
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    x0=[xn(1,node1);xn(1,node2)];
    if (nsd > 1)
        y0=[xn(2,node1);xn(2,node2)];
        if nsd ==3
            z0=[xn(3,node1);xn(3,node2)];
        end
    else
        y0=[0;0];
    end
    if nsd < 3
        if E(e)<= mean(E)
%      
        line(x0,y0,'Color','k','LineWidth',0.5,'LineStyle','--');
        else
      
       line(x0,y0,'Color','m','LineWidth',1,'LineStyle','-');
        end
    else
        if E(e)<= mean(E)
        line(x0,y0,z0,'Color','k','LineWidth',0.5,'LineStyle','--');
        else
        line(x0,y0,z0,'Color','m','LineWidth',1,'LineStyle','-');
        end
        
         grid on
     
    end
end
function plot_mesh_undeformed1(nel,ien,xn,nnp,nsd,E)
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    x0=[xn(1,node1);xn(1,node2)];
    if (nsd > 1)
        y0=[xn(2,node1);xn(2,node2)];
        if nsd ==3
            z0=[xn(3,node1);xn(3,node2)];
        end
    else
        y0=[0;0];
    end
    if nsd < 3
        
%      
        line(x0,y0,'Color','k','LineWidth',0.1,'LineStyle','--');
     
    else
     
        line(x0,y0,z0,'Color','k','LineWidth',0.05,'LineStyle','--');
       
        
         grid on
     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements - truss and beam %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers(nel,ien,xn,nnp,nsd)
% number the nodes in the undeformed configuration
for n=1:nnp
    if (nsd == 3)
        text(xn(1,n),xn(2,n), xn(3,n), num2str(n),'FontSize',14);
    end
    if (nsd == 2)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',14);
    end
    if (nsd == 1)
        text(xn(1,n),0, num2str(n),'FontSize',14);
    end
end
% number the elements in the undeformed configuration
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    %compute the coordinates of the middle point
    xg=0.5*(xn(:,node1)+xn(:,node2));
    s=sprintf('(%d)', e);
    if (nsd == 3)
        text(xg(1), xg(2), xg(3), s,'FontSize',14);
    end
    if (nsd == 2)
        text(xg(1),xg(2), s,'FontSize',14);
    end
    if (nsd == 1)
        text(xg(1),0, s,'FontSize',14);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements - triangular mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers_tria(nel,ien,xn,nnp,nsd)
% number the nodes
for n=1:nnp
    if (nsd > 1)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',12);
    else
        text(xn(1,n),0, num2str(n),'FontSize',12);
    end
end
% number the elements
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    node3=ien(3,e);
    %compute the coordinates of the center of mass
    xg=(1/3)*(xn(:,node1)+xn(:,node2)+xn(:,node3));
    s=sprintf('(%d)', e);
    text(xg(1),xg(2), s,'FontSize',12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements - quad mesh       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers_quad(nel,ien,xn,nnp,nsd)
% number the nodes
for n=1:nnp
    if (nsd > 1)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',12);
    else
        text(xn(1,n),0, num2str(n),'FontSize',12);
    end
end
% number the elements
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    node3=ien(3,e);
    node4=ien(4,e);
    %compute the coordinates of the center of mass
    xg=(1/4)*(xn(:,node1)+xn(:,node2)+xn(:,node3)+xn(:,node4));
    s=sprintf('(%d)', e);
    text(xg(1),xg(2), s,'FontSize',12);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% deformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnn]=plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,axial,Qe,E)
hold on
switch type
    case 'truss'
        scale=display_factor*Lcar/max(max(abs(Ucomp))); % scale factor for the displacements
        legend_comp_flag = 0;
        legend_ten_flag = 0;
        for e=1:nel
            node1=ien(1,e);
            node2=ien(2,e);
            xt(:,1)=xn(:,node1);
            xt(:,2)=xn(:,node2);
            
            if (nsd == 1)
                xt(2,1)=0;
                xt(2,2)=0;
            end
            
            for i=1:ndf
                if Idb(i,node1)~=0
                    xt(i,1)=xt(i,1)+scale*Ucomp(i,node1);
                else
                    xt(i,1)=xt(i,1)+scale*Ucomp(i,node1);
                end
                if Idb(i,node2)~=0
                    xt(i,2)=xt(i,2)+scale*Ucomp(i,node2);
                else
                    xt(i,2)=xt(i,2)+scale*Ucomp(i,node2);
                end
            end
            if nsd < 3
                xt1=xt(:,1);
                xt2=xt(:,2);           
                xnn(:,node1)=xt1;
                xnn(:,node2)=xt2;
                if (axial(1,e) > 0)
                    h(1,1) = plot(xt(1,:),xt(2,:), 'r-o','LineWidth',1.0);
                    
                    if legend_ten_flag == 0
                        legend_ten_flag = 1; 
                    end
                else
                    h(2,1) = plot(xt(1,:),xt(2,:), 'b-o','LineWidth',1.0);
                    if legend_comp_flag == 0
                        legend_comp_flag = 1;
                    end
                end
            end
            if nsd == 3
                xt1=xt(:,1);
                xt2=xt(:,2);           
                xnn(:,node1)=xt1;
                xnn(:,node2)=xt2;
                if (axial(1,e) > 0)
                    h(1,1) = plot3(xt(1,:),xt(2,:),xt(3,:), 'r','LineWidth',1.2);
                     grid on
                    if legend_ten_flag == 0
                        legend_ten_flag = 1;
                    end
                else
                    h(2,1) = plot3(xt(1,:),xt(2,:),xt(3,:), 'b','LineWidth',1.2);
                    if legend_comp_flag == 0
                        legend_comp_flag = 1;
                    end
                end
            end
        end
%         if legend_comp_flag && legend_ten_flag
%             legend(h,'Tension', 'Compression')
%         else
%             if legend_comp_flag
%                 legend(h(2,1),'Compression')
%             else
%                 if legend_ten_flag
%                     legend(h(1,1),'Tension')
%                 end
%             end
%         end
    case 'beam'
        xi=-1:0.1:1;
        nxi=size(xi,2);
        
        Umax=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                   SCALE FACTOR                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e=1:nel
            node1=ien(1,e);
            node2=ien(2,e);
            if nsd == 2
                % displacement + rotation
                Un(1,1)=Ucomp(1,node1);
                Un(1,2)=Ucomp(2,node1);
                Un(1,3)=Ucomp(3,node1);
                Un(2,1)=Ucomp(1,node2);
                Un(2,2)=Ucomp(2,node2);
                Un(2,3)=Ucomp(3,node2);
            elseif nsd == 3
                Un(1,1)=Ucomp(1,node1);
                Un(1,2)=Ucomp(2,node1);
                Un(1,3)=Ucomp(3,node1);
                Un(1,4)=Ucomp(4,node1);
                Un(1,5)=Ucomp(5,node1);
                Un(1,6)=Ucomp(6,node1);
                Un(1,7)=Ucomp(1,node2);
                Un(1,8)=Ucomp(2,node2);
                Un(1,9)=Ucomp(3,node2);
                Un(1,10)=Ucomp(4,node2);
                Un(1,11)=Ucomp(5,node2);
                Un(1,12)=Ucomp(6,node2);
            end
            
            % intial postion
            x0=zeros(nxi,nsd);
            for n=1:nxi
                x0(n,:)=((1-xi(n))/2)*xn(:,node1)'+((1+xi(n))/2)*xn(:,node2)';
            end
            
            %Le=sqrt((xn(1,node1)-xn(1,node2))^2+(xn(2,node1)-xn(2,node2))^2);
            Le = norm(xn(:,node2)-xn(:,node1));
            if nsd ==2
                % g(i,j) - vector i coordinate j
                g=zeros(nsd,nsd);
                g(1,:)=(1/Le)*(xn(:,node2)-xn(:,node1))';
                
                g(2,1)=-g(1,2);
                g(2,2)=g(1,1);
                
                % to make g_3 = z - so that the rotation are the same in both systems of coordinates
                cross_prod=g(1,1)*g(2,2)-g(1,2)*g(2,1);
                if cross_prod < 0
                    g(2,:)=-g(2,:);
                end
                
                %transformation matrices
                beta=g;
                betap=inv(beta);
                % transform displacement coordinates in the global system
                un=zeros(nen,nsd);
                
                for n=1:nen
                    for i=1:nsd
                        for j=1:nsd
                            un(n,i)=un(n,i)+Un(n,j)*betap(j,i);
                        end
                    end
                end
                
                for n=1:nen
                    un(n,3)=Un(n,3);
                end
            else
                
                %%%%%%%%%%%%%%%%
                
                % transform displacement coordinates in the global system
                un=zeros(12,1);
                
                un = Qe(:,:,e)*Un(1,:)';
                
                
            end
            %un = Qe(:,:,e)*Un;
            % transverse displacement in the local system of coordinates
            if nsd == 2
                for i=1:nxi
                    N1(i)=(1/4)*(1-xi(i))^2*(2+xi(i));
                    N2(i)=(Le/8)*(1-xi(i)^2)*(1-xi(i));
                    N3(i)=(1/4)*(1+xi(i))^2*(2-xi(i));
                    N4(i)=(Le/8)*(-1+xi(i)^2)*(1+xi(i));
                    u(i,2)=N1(i)*un(1,2)+N2(i)*un(1,3)+N3(i)*un(2,2)+N4(i)*un(2,3);
                end
                
                % horizontal displacement in the local system of coordinates
                for i=1:nxi
                    u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(2,1));
                end
            elseif nsd ==3
                for i=1:nxi
                    N1(i) = 1/2*(1-xi(i));
                    N2(i) = (1/4)*(1-xi(i))^2*(2+xi(i));
                    N3(i) = (1/4)*(1-xi(i))^2*(2+xi(i));
                    N4(i) = 1/2*(1-xi(i));
                    N5(i) = (Le/8)*(1-xi(i)^2)*(1-xi(i));
                    N6(i) = (Le/8)*(1-xi(i)^2)*(1-xi(i));
                    N7(i) = xi(i);
                    N8(i) = (1/4)*(1+xi(i))^2*(2-xi(i));
                    N9(i) = (1/4)*(1+xi(i))^2*(2-xi(i));
                    N10(i) = xi(i);
                    N11(i) = (Le/8)*(-1+xi(i)^2)*(1+xi(i));
                    N12(i) = (Le/8)*(-1+xi(i)^2)*(1+xi(i));
                    u(i,2)=N2(i)*un(2,1)+N6(i)*un(6,1)+...
                        +N8(i)*un(8,1)+N12(i)*un(12,1);
                    u(i,3)=N3(i)*un(3,1)+N5(i)*un(5,1)...
                        +N9(i)*un(9,1)+N11(i)*un(11,1);
                end
                
                % horizontal displacement in the local system of coordinates
                for i=1:nxi
                    u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(7,1));
                end
            end
            
            % transform from local to global system
            % transform from local to global system
            U=zeros(nxi,nsd);
            if nsd == 2
                for i=1:2
                    for n=1:nxi
                        for j=1:2
                            U(n,i)=U(n,i)+u(n,j)*beta(j,i);
                        end
                    end
                end
            else
                for n=1:nxi
                    U(n,:)= Qe(1:3,1:3,e)*u(n,:)';
                end
            end
            
            % scale factor
            if (max(max(abs(U)))>Umax)
                Umax=max(max(abs(U)));
            end
        end
        
        scale=display_factor*Lcar/Umax;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                       PLOT                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for e=1:nel
            node1=ien(1,e);
            node2=ien(2,e);
            if nsd == 2
                % displacement + rotation
                Un(1,1)=Ucomp(1,node1);
                Un(1,2)=Ucomp(2,node1);
                Un(1,3)=Ucomp(3,node1);
                Un(2,1)=Ucomp(1,node2);
                Un(2,2)=Ucomp(2,node2);
                Un(2,3)=Ucomp(3,node2);
            elseif nsd == 3
                Un(1,1)=Ucomp(1,node1);
                Un(1,2)=Ucomp(2,node1);
                Un(1,3)=Ucomp(3,node1);
                Un(1,4)=Ucomp(4,node1);
                Un(1,5)=Ucomp(5,node1);
                Un(1,6)=Ucomp(6,node1);
                Un(1,7)=Ucomp(1,node2);
                Un(1,8)=Ucomp(2,node2);
                Un(1,9)=Ucomp(3,node2);
                Un(1,10)=Ucomp(4,node2);
                Un(1,11)=Ucomp(5,node2);
                Un(1,12)=Ucomp(6,node2);
            end
            
            % intial postion
            x0=zeros(nxi,nsd);
            for n=1:nxi
                x0(n,:)=((1-xi(n))/2)*xn(:,node1)'+((1+xi(n))/2)*xn(:,node2)';
            end
            
            %Le=sqrt((xn(1,node1)-xn(1,node2))^2+(xn(2,node1)-xn(2,node2))^2);
            Le = norm(xn(:,node2)-xn(:,node1));
            if nsd ==2
                % g(i,j) - vector i coordinate j
                g=zeros(nsd,nsd);
                g(1,:)=(1/Le)*(xn(:,node2)-xn(:,node1))';
                
                g(2,1)=-g(1,2);
                g(2,2)=g(1,1);
                
                % to make g_3 = z - so that the rotation are the same in both systems of coordinates
                cross_prod=g(1,1)*g(2,2)-g(1,2)*g(2,1);
                if cross_prod < 0
                    g(2,:)=-g(2,:);
                end
                
                %transformation matrices
                beta=g;
                betap=inv(beta);
                % transform displacement coordinates in the global system
                un=zeros(nen,nsd);
                
                for n=1:nen
                    for i=1:nsd
                        for j=1:nsd
                            un(n,i)=un(n,i)+Un(n,j)*betap(j,i);
                        end
                    end
                end
                
                for n=1:nen
                    un(n,3)=Un(n,3);
                end
            else
                
                %%%%%%%%%%%%%%%%     3D
                % transform displacement coordinates in the global system
                un=zeros(12,1);
                
                un = Qe(:,:,e)*Un(1,:)';
                
                
            end
            %un = Qe(:,:,e)*Un;
            % transverse displacement in the local system of coordinates
            if nsd == 2
                for i=1:nxi
                    N1(i)=(1/4)*(1-xi(i))^2*(2+xi(i));
                    N2(i)=(Le/8)*(1-xi(i)^2)*(1-xi(i));
                    N3(i)=(1/4)*(1+xi(i))^2*(2-xi(i));
                    N4(i)=(Le/8)*(-1+xi(i)^2)*(1+xi(i));
                    u(i,2)=N1(i)*un(1,2)+N2(i)*un(1,3)+N3(i)*un(2,2)+N4(i)*un(2,3);
                end
                
                % horizontal displacement in the local system of coordinates
                for i=1:nxi
                    u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(2,1));
                end
            elseif nsd ==3
                for i=1:nxi
                    N1(i) = 1/2*(1-xi(i));
                    N2(i) = (1/4)*(1-xi(i))^2*(2+xi(i));
                    N3(i) = (1/4)*(1-xi(i))^2*(2+xi(i));
                    N4(i) = 1/2*(1-xi(i));
                    N5(i) = (Le/8)*(1-xi(i)^2)*(1-xi(i));
                    N6(i) = (Le/8)*(1-xi(i)^2)*(1-xi(i));
                    N7(i) = xi(i);
                    N8(i) = (1/4)*(1+xi(i))^2*(2-xi(i));
                    N9(i) = (1/4)*(1+xi(i))^2*(2-xi(i));
                    N10(i) = xi(i);
                    N11(i) = (Le/8)*(-1+xi(i)^2)*(1+xi(i));
                    N12(i) = (Le/8)*(-1+xi(i)^2)*(1+xi(i));
                    u(i,2)=N2(i)*un(2,1)+N6(i)*un(6,1)+...
                        +N8(i)*un(8,1)+N12(i)*un(12,1);
                    u(i,3)=N3(i)*un(3,1)+N5(i)*un(5,1)...
                        +N9(i)*un(9,1)+N11(i)*un(11,1);
                end
                
                % horizontal displacement in the local system of coordinates
                for i=1:nxi
                    u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(7,1));
                end
            end
            
            % transform from local to global system
            % transform from local to global system
            U=zeros(nxi,nsd);
            if nsd == 2
                for i=1:2
                    for n=1:nxi
                        for j=1:2
                            U(n,i)=U(n,i)+u(n,j)*beta(j,i);
                        end
                    end
                end
            else
                for n=1:nxi
                    U(n,:)= inv(Qe(1:3,1:3,e))*u(n,:)';
                end
            end
            
            % deformed configuration in the global system
            xt=x0+scale*U;
            xt1=xt(1,:)';
            xt2=xt(21,:)';           
            xnn(:,node1)=xt1;
            xnn(:,node2)=xt2;
%              Le = norm(xn(:,node2)-xn(:,node1));
            if nsd == 2
                if E(e)<= mean(E)
                plot(xt(:,1),xt(:,2),'b-','LineWidth',0.8);                
                else
                plot(xt(:,1),xt(:,2),'m-','LineWidth',1.2);    
                end
            else
                
                if E(e)<= mean(E)
                    plot3(xt(:,1),xt(:,2),xt(:,3), 'b-','LineWidth',1.2);
                    
                else
                    plot3(xt(:,1),xt(:,2),xt(:,3), 'm-','LineWidth',1.2);
                    grid on
                end
            end
            
        end
        
    otherwise
        disp('type not supported by plot_results');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on displacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd)
alpha=display_factor*Lcar; % scale factor  for bc symbols
switch type
    case 'truss'
        for P=1:nnp
            if (nsd >1)
                if ((Idb(1,P) == 1) && (Idb(2,P) == 1))
                    bc_symbols(xn(:,P),alpha,3,nsd);
                end
                if ((Idb(1,P) ==1 ) && (Idb(2,P) == 0))
                    bc_symbols(xn(:,P),alpha,2,nsd);
                end
                if ((Idb(1,P) == 0) && (Idb(2,P) ~= 0))
                    bc_symbols(xn(:,P),alpha,1,nsd);
                end
                 if ((Idb(1,P) == -1) && (Idb(2,P) ==-1))
                    bc_symbols(xn(:,P),alpha,7,nsd);
                end
                 if ((Idb(1,P) ==-1 ) && (Idb(2,P) == 0))
                    bc_symbols(xn(:,P),alpha,8,nsd);
                end
                      
            else
                if (Idb(1,P) ~= 0)
                    bc_symbols([xn(1,P),0],alpha,2,nsd);
                end
            end
            
            
        end
    case 'beam'
        for P=1:nnp
            if ((Idb(1,P) > 0) && (Idb(2,P) > 0) && (Idb(3,P) == 0)) %try draw positive only
                bc_symbols(xn(:,P),alpha,3,nsd);
            end
            if ((Idb(1,P) == 1) && (Idb(2,P) == 0) && (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,2,nsd);
            end
            if ((Idb(1,P) == 0) && (Idb(2,P) ~= 0) && (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,1,nsd);
            end
            if ((Idb(1,P) == 0) && (Idb(2,P) ~= 0) && (Idb(3,P) ~= 0))
                bc_symbols(xn(:,P),alpha,4,nsd);
            end
            if ((Idb(1,P) == 1) && (Idb(2,P) == 0) && (Idb(3,P) == 1))
                bc_symbols(xn(:,P),alpha,5,nsd);
            end
            if ((Idb(1,P) == 1) && (Idb(2,P) ==1) && (Idb(3,P) ==1))
                bc_symbols(xn(:,P),alpha,6,nsd);
            end
            
            if ((Idb(1,P) == -1) && (Idb(2,P)== -1) && (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,7,nsd);
            end

            if ((Idb(1,P) == -1) && (Idb(2,P)== 0) && (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,8,nsd);
            end
            
            if ((Idb(1,P) == -1) && (Idb(2,P)==-1 ) && (Idb(3,P) == -1))
                bc_symbols(xn(:,P),alpha,9,nsd);
            end
            
            if ((Idb(1,P) == -1) && (Idb(2,P)== 0 ) && (Idb(3,P) == -1))
                bc_symbols(xn(:,P),alpha,10,nsd);
            end
        end
    otherwise
        disp('type not supported by plot_results');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on Prescribed displacement %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_predisplacement(type,Lcar,g,display_factor,nnp,xn,nsd)
delta=display_factor*Lcar/max(max(abs(g))/2); % scale factor for force b.c.
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    if ( nsd == 3)
        if ( (g(1,N) ~= 0 ) || (g(2,N) ~= 0) || (g(3,N) ~= 0))
            quiver3(xn(1,N),xn(2,N),xn(3,N),g(1,N),g(2,N),g(3,N),delta,'c:','LineWidth',3);
        end
    end
    if ( nsd == 2)
        if ( (g(1,N) ~= 0 ) || (g(2,N) ~= 0))
           quiver(xn(1,N),xn(2,N),g(1,N),g(2,N),delta,'c:','LineWidth',3);
        end
    end
    if ( nsd == 1)
        if (f(1,N) ~=0)
            quiver(xn(1,N),0,f(1,N),0,delta,'r','LineWidth',2);
        end
    end
end

if (strcmp(type, 'beam'))&& nsd==3
    for N=1:nnp
        if (g(3,N) > eps)
            plot_moments([xn(1,N);xn(2,N);xn(3,N)],alpha,0,nsd);
        end
        if (g(3,N) < -eps)
            plot_moments([xn(1,N);xn(2,N);xn(3,N)],alpha,1,nsd);
        end
    end
elseif (strcmp(type, 'beam'))&& nsd==2
    for N=1:nnp
        if (g(3,N) > eps)
            plot_moments([xn(1,N);xn(2,N)],alpha,0,nsd);
        end
        if (g(3,N) < -eps)
            plot_moments([xn(1,N);xn(2,N)],alpha,1,nsd);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on force %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd)
delta=display_factor*Lcar/max(max(abs(f))); % scale factor for force b.c.
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    if ( nsd == 3)
        if ( (f(1,N) ~= 0 ) || (f(2,N) ~= 0) || (f(3,N) ~= 0))
            quiver3(xn(1,N),xn(2,N),xn(3,N),f(1,N),f(2,N),f(3,N),delta,'r','LineWidth',2);
        end
    end
    if ( nsd == 2)
        if ( (f(1,N) ~= 0 ) || (f(2,N) ~= 0))
           quiver(xn(1,N),xn(2,N),f(1,N),f(2,N),delta,'r','LineWidth',2);
        end
    end
    if ( nsd == 1)
        if (f(1,N) ~=0)
            quiver(xn(1,N),0,f(1,N),0,delta,'r','LineWidth',2);
        end
    end
end

if (strcmp(type, 'beam'))&& nsd==3
    for N=1:nnp
        if (f(3,N) > eps)
            plot_moments([xn(1,N);xn(2,N);xn(3,N)],alpha,0,nsd);
        end
        if (f(3,N) < -eps)
            plot_moments([xn(1,N);xn(2,N);xn(3,N)],alpha,1,nsd);
        end
    end
elseif (strcmp(type, 'beam'))&& nsd==2
    for N=1:nnp
        if (f(3,N) > eps)
            plot_moments([xn(1,N);xn(2,N)],alpha,0,nsd);
        end
        if (f(3,N) < -eps)
            plot_moments([xn(1,N);xn(2,N)],alpha,1,nsd);
        end
    end
end

%%%%%%%%%%%%%
% reactions %
%%%%%%%%%%%%%
function plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
beta=3*display_factor*Lcar/max(max(abs(Rcomp))); % scale factor for reactions
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    RN=zeros(ndf);
    if (nsd == 3)
        if ( (Idb(1,N) ~= 0 ) || (Idb(2,N) ~= 0) || (Idb(3,N) ~= 0))
            h = quiver3(xn(1,N),xn(2,N),xn(3,N),Rcomp(1,N),Rcomp(2,N),Rcomp(3,N),beta,'m-o','LineWidth', 1);
            %adjust_quiver_arrowhead_size(h, 0.5/beta);
        end
    end
    if (nsd == 2)
        if ((Idb(1,N) ~= 0) || (Idb(2,N) ~= 0))
            h = quiver(xn(1,N),xn(2,N),Rcomp(1,N),Rcomp(2,N),beta,'r-o', 'LineWidth', 1);
            %adjust_quiver_arrowhead_size(h, 0.5/beta);
        end
        if (strcmp(type, 'beam') && (Idb(3,N) ~= 0))
            if (Rcomp(3,N) > eps)
                plot_moments([xn(1,N);xn(2,N)],alpha,0);
            end
            if (Rcomp(3,N) < -eps)
                plot_moments([xn(1,N);xn(2,N)],alpha,1);
            end
        end
    end
    if (nsd == 1)
        if (Idb(1,N) ~= 0)
            h = quiver(xn(1,N),0,Rcomp(1,N),0,beta,'r', 'LineWidth', 1);
            adjust_quiver_arrowhead_size(h, 0.5/beta);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              BOUNDARY CONDITIONS SYMBOLS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bc_symbols(xp,alpha,symbol,nsd)
if (nsd < 3)
    switch symbol
        case 1
            % v fixed
            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3
                circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(7/8)*alpha],alpha/8);
            end
            
            
        case 2
            % u fixed
            x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
            y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3
                circle([xp(1)-(7/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
            end
            
        case 3
            % u and v fixed
            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3
                line([xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)], ...
                    [xp(2)-(3/4)*alpha;xp(2)-alpha],'Color','k','LineWidth',1.2);
            end
            
        case 4
            % v and theta fixed
            x=[xp(1)-alpha/2;xp(1)+alpha/2];
            y=[xp(2);xp(2)];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3
                circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(1/8)*alpha],alpha/8);
            end
            
        case 5
            % u and theta fixed
            x=[xp(1);xp(1)];
            y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
            
            line(x,y,'Color','k','LineWidth',1.2);
            
            for i=0:3
                circle([xp(1)-(1/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
            end
            
        case 6
            % u, v and theta fixed
            line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],'Color','k','LineWidth',1.2);
            for i=0:3
                line([xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4],[xp(2),xp(2)-alpha/4]...
                    ,'Color','k','LineWidth',1.2);
            end
        case 7
            % u, v fixed reverse up-down
            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
            xc = xp(1) ;  % Rotate about the 1/4 chord point
            yc = xp(2) ;
            ang=180;
            xr =  (x-xc)*cosd(ang) + (y-yc)*sind(ang) + xc;
            yr = -(x-xc)*sind(ang) + (y-yc)*cosd(ang) + yc;
            line(xr,yr,'Color','k','LineWidth',1.2);
                       
            for i=0:3   
                xx=[xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)];
                yy=[xp(2)-(3/4)*alpha;xp(2)-alpha];
                 xd =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                 yd = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                line(xd,yd,'Color','k','LineWidth',1.2);
            end 
            
            case 8
            % u fixed reverse L-R
            
            x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
            y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
            xc = xp(1) ;  % Rotate about the 1/4 chord point
            yc = xp(2) ;
            ang=180;
            xr =  (x-xc)*cosd(ang) + (y-yc)*sind(ang) + xc;
            yr = -(x-xc)*sind(ang) + (y-yc)*cosd(ang) + yc;
            line(xr,yr,'Color','k','LineWidth',1.2);
            
            for i=0:3
                xx=xp(1)-(7/8)*alpha;
                yy=xp(2)-(3/8)*alpha+i*alpha/4;
                xr =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                yr = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                circle([xr,yr],alpha/8);
            end
            
            case 9
              % u, v and theta fixed reverse up-down
               line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],'Color','k','LineWidth',1.2);
                 xc = xp(1) ;  % Rotate about the 1/4 chord point
                 yc = xp(2) ;
                 ang=180;
                for i=0:3
                 xx=[xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4];
                 yy=[xp(2),xp(2)-alpha/4];
                 xd =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                 yd = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                    line(xd,yd,'Color','k','LineWidth',1.2);
                end
                
        case 10
                 % u and theta fixed
            x=[xp(1);xp(1)];
            y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
            
            line(x,y,'Color','k','LineWidth',1.2);
            xc = xp(1) ;  % Rotate about the 1/4 chord point
            yc =xp(2) ;
            ang=180;
            for i=0:3
                xx=xp(1)-(1/8)*alpha;
                yy=xp(2)-(3/8)*alpha+i*alpha/4;
                xr =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                yr = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                circle([xr;yr],alpha/8);
            end   
                
          
                     
    end
end
if (nsd ==3)
        switch symbol
            case 1
                % v fixed
                x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
                y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
                z=[xp(3);xp(3);xp(3);xp(3)];
                line(x,y,z,'Color','k','LineWidth',1.2);
    
                for i=0:3
                    circle1([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(7/8)*alpha;xp(3)],alpha/8);
                end
    
    
            case 2
                % u fixed
                x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
                y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
                z=[xp(3);xp(3);xp(3);xp(3)];
                line(x,y,z,'Color','k','LineWidth',1.2);
             
    
                for i=0:3
                    circle1([xp(1)-(7/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4;xp(3)],alpha/8);
                end
    
            case 3
                % u and v fixed
                x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
                y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
                z=[xp(3);xp(3);xp(3);xp(3)];
                line(x,y,z,'Color','k','LineWidth',1.2);
    
                for i=0:3
                    line([xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)], ...
                        [xp(2)-(3/4)*alpha;xp(2)-alpha],[xp(3),xp(3)],'Color','k','LineWidth',1.2);
                end
    
            case 4
                % v and theta fixed
                x=[xp(1)-alpha/2;xp(1)+alpha/2];
                y=[xp(2);xp(2)];
                z=[xp(3);xp(3)];
                line(x,y,z,'Color','k','LineWidth',1.2);
               
    
                for i=0:3
                    circle1([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(1/8)*alpha;xp(3)],alpha/8);
                end
    
            case 5
                % u and theta fixed
                x=[xp(1);xp(1)];
                y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
                z=[xp(3);xp(3)];
                line(x,y,z,'Color','k','LineWidth',1.2);
    
                for i=0:3
                    circle1([xp(1)-(1/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4;xp(3)],alpha/8);
                end
    
            case 6
                % u, v and theta fixed
                line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],[xp(3),xp(3)],'Color','k','LineWidth',1.2);
                
                for i=0:3
                    line([xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4],[xp(2),xp(2)-alpha/4]...
                           ,[xp(3),xp(3)],'Color','k','LineWidth',1.2);
                end
             
            case 7
                % u and v fixed reverse up-down
                x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
                y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
                z=[xp(3);xp(3);xp(3);xp(3)];
                xc = xp(1) ;  % Rotate about the 1/4 chord point
                yc = xp(2) ;
                ang=180;
                xr =  (x-xc)*cosd(ang) + (y-yc)*sind(ang) + xc;
                yr = -(x-xc)*sind(ang) + (y-yc)*cosd(ang) + yc;
                line(xr,yr,z,'Color','k','LineWidth',1.2);
                
                for i=0:3
                    xx=[xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)];
                    yy=[xp(2)-(3/4)*alpha;xp(2)-alpha];
                    xd =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                    yd = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                    line(xd,yd,[xp(3),xp(3)],'Color','k','LineWidth',1.2);
                end
                
            case 8
                % u fixed reverse L-R
                x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
                y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
                z=[xp(3);xp(3);xp(3);xp(3)];
                xc = xp(1) ;  % Rotate about the 1/4 chord point
                yc = xp(2) ;
                ang=180;
                xr =  (x-xc)*cosd(ang) + (y-yc)*sind(ang) + xc;
                yr = -(x-xc)*sind(ang) + (y-yc)*cosd(ang) + yc;
                line(xr,yr,z,'Color','k','LineWidth',1.2);
                
                
                for i=0:3
                    xx=xp(1)-(7/8)*alpha;
                    yy=xp(2)-(3/8)*alpha+i*alpha/4;
                    xr =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                    yr = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                    circle1([xr;yr;xp(3)],alpha/8);
                end
                
            case 9
                % u, v and theta fixed
                line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],[xp(3),xp(3)],'Color','k','LineWidth',1.2);
                 xc = xp(1) ;  % Rotate about the 1/4 chord point
                 yc = xp(2) ;
                 ang=180;
                for i=0:3
                 xx=[xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4];
                 yy=[xp(2),xp(2)-alpha/4];
                 xd =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                 yd = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                 line(xd,yd,[xp(3),xp(3)],'Color','k','LineWidth',1.2);
                end
                
            case 10
                % u and theta fixed
                x=[xp(1);xp(1)];
                y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
                
                line(x,y,'Color','k','LineWidth',1.2);
                xc = xp(1) ;  % Rotate about the 1/4 chord point
                yc =xp(2) ;
                ang=180;
                for i=0:3
                    xx=xp(1)-(1/8)*alpha;
                    yy=xp(2)-(3/8)*alpha+i*alpha/4;
                    xr =  (xx-xc)*cosd(ang) + (yy-yc)*sind(ang) + xc;
                    yr = -(xx-xc)*sind(ang) + (yy-yc)*cosd(ang) + yc;
                     
                  circle1([xr;yr;xp(3)],alpha/8);
                
                end
                  
                
                
        end
end

function circle(x0,r)
theta=0:0.1:2*pi;
x=r*cos(theta)+x0(1);
y=r*sin(theta)+x0(2);

plot(x,y,'k','LineWidth',1.2);

function circle1(x0,r)
theta=0:0.1:2*pi;
x=r*cos(theta)+x0(1);
y=r*sin(theta)+x0(2);
z=ones(1,63)*x0(3);
plot3(x,y,z,'k','LineWidth',1.2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       MOMENTS                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_moments(xp,alpha,sign,nsd)
 if nsd==3
switch sign
    
    case 0  %positive moment
       
        r=alpha/4;
        theta=0:0.1:3*pi/2;
        x=r*cos(theta)+xp(1);
        y=r*sin(theta)+xp(2);
        z=ones(1,48)*xp(3);
        
        plot3(x,y,z,'k','LineWidth',1.2);
        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],[xp(3),xp(3)],...
            'Color','k','LineWidth',1.2);
        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],[xp(3),xp(3)],...
            'Color','k','LineWidth',1.2);
        
        
    case 1  % negative moment
        r=alpha/4;
        theta=pi:-0.1:-pi/2;
        x=r*cos(theta)+xp(1);
        y=r*sin(theta)+xp(2);
        z=ones(1,48)*xp(3);
   
        plot(x,y,'k','LineWidth',1.2);
        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],[xp(3),xp(3)],...
            'Color','k','LineWidth',1.2);
        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],[xp(3),xp(3)],...
            'Color','k','LineWidth',1.2);
    
   end  
end    
           
        if nsd==2
            switch sign
      case 0  %positive moment
       
        r=alpha/4;
        theta=0:0.1:3*pi/2;
        x=r*cos(theta)+xp(1);
        y=r*sin(theta)+xp(2);
        
        
        plot(x,y,'k','LineWidth',1.2);
        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
            'Color','k','LineWidth',1.2);
        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
            'Color','k','LineWidth',1.2);
        
        
    case 1  % negative moment
        r=alpha/4;
        theta=pi:-0.1:-pi/2;
        x=r*cos(theta)+xp(1);
        y=r*sin(theta)+xp(2);
        
   
        plot(x,y,'k','LineWidth',1.2);
        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
            'Color','k','LineWidth',1.2);
        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
            'Color','k','LineWidth',1.2);
            
    end
       

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             HEAT TEMPERATURE DISTRIBUTION                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_temp_tria(Ucomp,nel,ien,xn,nen,nnp,nsd,flag_number)
%%%%%%%%%%%%%
% tria mesh %
%%%%%%%%%%%%%
if (nen == 3)
    for e=1:nel,
        for k=1:nen,
            node(k)=ien(k,e);
        end
        X(e,:,:)=[[xn(1,node(1))   xn(1,node(3))]; [xn(1,node(2))  xn(1,node(3))]];
        Y(e,:,:)=[[xn(2,node(1))   xn(2,node(3))]; [xn(2,node(2))  xn(2,node(3))]];
        Z(e,:,:)=[[Ucomp(node(1))  Ucomp(node(3))];[Ucomp(node(2)) Ucomp(node(3))]];
        
    end
    for e=1:nel,
        for i=1:2,
            for j=1:2,
                XX(i,j)=X(e,i,j);
                YY(i,j)=Y(e,i,j);
                ZZ(i,j)=Z(e,i,j);
            end
        end
        contourf(XX,YY,ZZ,5);
        %         if (flag_number == 1)
        %            [C,h]=contourf(XX,YY,ZZ,5);
        %             clabel(C,h);
        %         end
    end
    colorbar;
    
else
    %%%%%%%%%%%%%
    % quad mesh %
    %%%%%%%%%%%%%
    for e=1:nel,
        for k=1:nen,
            node(k)=ien(k,e);
        end
        X(e,:,:)=[[xn(1,node(1))  xn(1,node(4))]; [xn(1,node(2))  xn(1,node(3))]];
        Y(e,:,:)=[[xn(2,node(1))  xn(2,node(4))]; [xn(2,node(2))  xn(2,node(3))]];
        Z(e,:,:)=[[Ucomp(node(1)) Ucomp(node(4))];[Ucomp(node(2)) Ucomp(node(3))]];
        
    end
    for e=1:nel,
        for i=1:2,
            for j=1:2,
                XX(i,j)=X(e,i,j);
                YY(i,j)=Y(e,i,j);
                ZZ(i,j)=Z(e,i,j);
            end
        end
        contourf(XX,YY,ZZ,5);
        %         if (flag_number == 1)
        %             [C,h]=contour(XX,YY,ZZ,5);
        %             clabel(C,h);
        %        end
    end
    colorbar;
    
end

if (flag_number == 1)
    for n=1:nnp
        if (nsd > 1)
            text(xn(1,n),xn(2,n), num2str(Ucomp(n)),'FontSize',12);
        else
            text(xn(1,n),0, num2str(Ucomp(n)),'FontSize',12);
        end
    end
end

function N = shape_tria(n,xi,eta)
switch n
    case 1
        N=xi;
    case 2
        N=eta;
    case 3
        N=1-xi-eta;
end




%     n=2; % number of points on the parent element
%     h=1/(n-1);
%
%     % points on the parent element
%     for i=1:n,
%         xi(1,i)=1-(h*(i-1));
%     end;
%
%     for j=1:n
%         eta(2,j)=h*(j-1);
%     end;
%
%     for k=1:(2*n-1),
%         xi(3,k)=1-(h/2)*(k-1);
%         eta(3,k)=1-xi(3,k);
%     end;
%
%     m=0;
%     alpha=((2*n-1)-rem(2*n-1,2))/2
%     for p=1:alpha,
%         for q=1:n,
%             t=(q-1)*h;
%             m=m+1;
%             xi(4,m)=(1-t)*xi(1,p)+t*xi(3,p);
%             eta(4,m)=t*eta(3,p);
%         end;
%     end;
%
%     for p=alpha+1:(2*n-1),
%         for q=1:n,
%             m=m+1;
%             t=(q-1)*h;
%             xi(4,m)=t*xi(3,p);
%             eta(4,m)=(1-t)*eta(2,p-alpha)+t*eta(3,p);
%         end;
%     end;
%
%     m=0;
%     for i=1:2*n-1,
%         for j=1:n,
%             m=m+1;
%             XI(i,j)=xi(4,m);
%             ETA(i,j)=eta(4,m);
%         end;
%     end;
%
%     % actual plot
%
%     for e=1:nel,
%         X=zeros(2*n-1,n);
%         Y=zeros(2*n-1,n);
%         Z=zeros(2*n-1,n);
%
%         for k=1:nen,
%             node(k)=ien(k,e);
%         end
%
%         for i=1:2*n-1,
%             for j=1:n,
%                 for k=1:nen,
%                     X(i,j)=X(i,j)+xn(1,node(k))*shape_tria(k,XI(i,j),ETA(i,j));
%                     Y(i,j)=Y(i,j)+xn(2,node(k))*shape_tria(k,XI(i,j),ETA(i,j));
%                     Z(i,j)=Z(i,j)+Ucomp(node(k))*shape_tria(k,XI(i,j),ETA(i,j));
%                 end;
%             end;
%         end;
%
%
%         contourf(X,Y,Z,5);
%         if (flag_number == 1)
%             [C,h]=contour(X,Y,Z,5);
%             clabel(C,h);
%         end
%         colorbar;
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shear_diagram(xn,ien,nen,axial,nsd,nel,nnp)
%------------------------------------------------------------------------
% Purpose:
% Shear  diagarma
%
% Synopsis:
%
%
% Variable Description:
% nen - number of nodes per element
%
% nsd - number of spatial dimensions
%------------------------------------------------------------------------

if nsd ==2
    for e=1:nel
        node1=ien(1,e);
        node2=ien(2,e);
        x0=[xn(1,node1);xn(1,node2)];
        shear(1,1)  = axial(2,e);
        shear(1,2)  = -axial(5,e);
        plot(x0,shear,'m-','LineWidth',0.8);  
        if e~=nel
            shear(1,1)  = -axial(5,e);
            shear(1,2)  = axial(2,e+1);
            x0=[xn(1,node2);xn(1,node2)];
            plot(x0,shear,'m-','LineWidth',0.8);
        end
    end
end

function moment_diagram(xn,ien,nen,axial,nsd,nel,nnp,fnn)
%------------------------------------------------------------------------
% Purpose:
%  Moment diagarma
%
% Synopsis:
%
%
% Variable Description:
% nen - number of nodes per element
%
% nsd - number of spatial dimensions
%------------------------------------------------------------------------

if nsd ==2
    for e=1:nel
        node1=ien(1,e);
        node2=ien(2,e);
        x0=[xn(1,node1);xn(1,node2)];
        moment(1,1) = -axial(3,e);
        moment(1,2) = axial(6,e);
        if norm(fnn(:,e)) == 0
            plot(x0,moment,'r-','LineWidth',0.8);
        else
            V1  = axial(2,e);
            V2  = -axial(5,e);
            a = (V2-V1)/2/(x0(2)-x0(1));
            b = V1-2*a*x0(1);
            c = moment(1,1)-a*x0(1)^2-b*x0(1);
            x_div = (x0(2)-x0(1))/100;
            x_x0 = x0(1):x_div:x0(2);
            Moment = a*x_x0.^2+b*x_x0+c;
            plot(x_x0,Moment,'r-','LineWidth',0.8);
        end
        if e~=nel
            moment(1,1)  = axial(6,e);
            moment(1,2)  = -axial(3,e+1);
            x0=[xn(1,node2);xn(1,node2)];
            plot(x0,moment,'r-','LineWidth',0.8);
        end
    end
end

function adjust_quiver_arrowhead_size(quivergroup_handle, scaling_factor)
% Make quiver arrowheads bigger or smaller.
%
% adjust_quiver_arrowhead_size(quivergroup_handle, scaling_factor)
%
% Example:
%   h = quiver(1:100, 1:100, randn(100, 100), randn(100, 100));
%   adjust_quiver_arrowhead_size(h, 1.5);   % Makes all arrowheads 50% bigger.
%
% Inputs:
%   quivergroup_handle      Handle returned by "quiver" command.
%   scaling_factor          Factor by which to shrink/grow arrowheads.
%
% Output: none

% Kevin J. Delaney
% December 21, 2011
% BMT Scientific Marine Services (www.scimar.com)

% Steve Kelly (kd2cca@gmail.com)
% 2014-04-09
% Worcester Polytechnic Institute

if ~exist('quivergroup_handle', 'var')
    help(mfilename);
    return
end

if isempty(quivergroup_handle) || any(~ishandle(quivergroup_handle))
    errordlg('Input "quivergroup_handle" is empty or contains invalid handles.', ...
        mfilename);
    return
end

if length(quivergroup_handle) > 1
    errordlg('Expected "quivergroup_handle" to be a single handle.', mfilename);
    return
end

if ~strcmpi(get(quivergroup_handle, 'Type'), 'hggroup')
    errrodlg('Input "quivergroup_handle" is not of type "hggroup".', mfilename);
    return
end

if ~exist('scaling_factor', 'var') || ...
        isempty(scaling_factor) || ...
        ~isnumeric(scaling_factor)
    errordlg('Input "scaling_factor" is missing, empty or non-numeric.', ...
        mfilename);
    return
end

if length(scaling_factor) > 1
    errordlg('Expected "scaling_factor" to be a scalar.', mfilename);
    return
end

if scaling_factor <= 0
    errordlg('"Scaling_factor" should be > 0.', mfilename);
    return
end

line_handles = get(quivergroup_handle, 'Children');

if isempty(line_handles) || (length(line_handles) < 3) || ...
        ~ishandle(line_handles(2)) || ~strcmpi(get(line_handles(2), 'Type'), 'line')
    errordlg('Unable to adjust arrowheads.', mfilename);
    return
end

arrowhead_line = line_handles(2);

XData = get(arrowhead_line, 'XData');
YData = get(arrowhead_line, 'YData');
ZData = get(arrowhead_line, 'ZData');

if isempty(XData) || isempty(YData)
    return
end

%   Break up XData, YData into triplets separated by NaNs.
first_nan_index = find(~isnan(XData), 1, 'first');
last_nan_index  = find(~isnan(XData), 1, 'last');

for index = first_nan_index : 4 : last_nan_index
    these_indices = index + (0:2);
    
    if these_indices(end) > length(XData)
        break
    end
    
    x_triplet = XData(these_indices);
    y_triplet = YData(these_indices);
    
    if any(isnan(x_triplet)) || any(isnan(y_triplet))
        continue
    end
    
    %   First pair.
    delta_x = diff(x_triplet(1:2));
    delta_y = diff(y_triplet(1:2));
    x_triplet(1) = x_triplet(2) - (delta_x * scaling_factor);
    y_triplet(1) = y_triplet(2) - (delta_y * scaling_factor);
    
    %   Second pair.
    delta_x = diff(x_triplet(2:3));
    delta_y = diff(y_triplet(2:3));
    x_triplet(3) = x_triplet(2) + (delta_x * scaling_factor);
    y_triplet(3) = y_triplet(2) + (delta_y * scaling_factor);
    
    XData(these_indices) = x_triplet;
    YData(these_indices) = y_triplet;
end

set(arrowhead_line, 'XData', XData, 'YData', YData);
