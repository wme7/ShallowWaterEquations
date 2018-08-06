function [qE,xE] = ExactRiemannSWE_Toro2001(x,TimeOut,ul,ur,dl,dr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Exact Riemann solver of the Saint Venant (shallow water) 
%
%                     coded by Manuel A. Diaz [2018.08.14]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: 
% [1] Toro, Eleuterio F., and Eleuterio Toro. Shock-capturing methods for
%     free-surface shallow flows. New York: John Wiley, 2001. pp:128. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global g; %gravity
global cl cr;
global nIter TOL nCells ChalLen;
global gate x0; %position of the gate

%----------------parse initial data----------------
nIter=50;
TOL=1E-6;
nCells=500;
x0 = x(1); % starting point of domain
ChalLen=x(end)-x(1); %channel lenght
gate=0.5*(x(end)+x(1)); %gate position

%---------compute celerity on left and right states--------
cl=(g*dl)^0.5;
cr=(g*dr)^0.5;

%---------check the depth positivity condition---------
dCrit=(ur-ul)-2*(cl+cr);
if (dl<=0||dr<=0||dCrit>=0)	%dry bed cases
    [xE,D,U] = drybed(TimeOut,ul,ur,dl,dr);
else
    [xE,D,U] = wetbed(TimeOut,ul,ur,dl,dr);
end

%----------------------------plot/debug------------------------------------
% figure(2);
% subplot(2,1,1); plot(xpos,D,'b-'); title('Depth');
% subplot(2,1,2); plot(xpos,U,'b-'); title('Velocity');

%------------------------------output--------------------------------------
qE = [D;U];

end 

function [xpos,D,U] = drybed(TimeOut,ul,ur,dl,dr)
    %to compute the exact solution in the case in which a portion of dry
    %bed is present
    global nCells;
    global gate;	%position of the gate
    global ChalLen;
    global x0;

    D=zeros(1,nCells);
    U=zeros(1,nCells);
    xpos=zeros(1,nCells);

    for i=1:1:nCells
        xpos(i)=x0 + i*ChalLen/nCells;
        xcoord=xpos(i)-gate;
        s=xcoord/TimeOut;

        if dl<=0
            %left state is dry
            [D(i),U(i)]=SamLef(s,ul,ur,dl,dr);
        else
            if dr<=0
                %right state is dry
                [D(i),U(i)]=SamRig(s,ul,ur,dl,dr);
            else
                %middle state is dry
                [D(i),U(i)]=SamMid(s,ul,ur,dl,dr);
            end
        end
    end
end

function [xpos,D,U] = wetbed(TimeOut,ul,ur,dl,dr)
    %solve the Riemann problem exactly for the wet-bed case
    global nIter TOL nCells ChalLen;
    global cl cr;
    global g;
    global gate;         %position of the gate
    global x0;

    depth0=starte(ul,ur,dl,dr);	%get start value of depth in star region
    ds=depth0;

    D=zeros(1,nCells);
    U=zeros(1,nCells);
    xpos=zeros(1,nCells);

    %-------start iteration-------
    for i=1:1:nIter
        [FL,FLD]=GEOFUN(ds,dl,cl);
        [FR,FRD]=GEOFUN(ds,dr,cr);
        ds=ds-(FL+FR+ur-ul)/(FLD+FRD);
        cha=abs(ds-depth0)/(0.5*(ds+depth0));         %convergence criteria
        if cha<=TOL
            break;
        end

        if ds<0
            ds=TOL;
        end

        depth0=ds;
    end

    if i>=nIter
        disp('number of iterations exceeded, STOP');
        pause;
    end

    %------compute velocity US in star region------
    us=0.5*(ul+ur)+0.5*(FR-FL);
    cs=(g*ds)^0.5;
    %----------evaluate exact solution at time TimeOut-----------
    for i=1:nCells
        xpos(i)= x0 + i*ChalLen/nCells;
        xcoord=xpos(i)-gate;
        s=xcoord/TimeOut;
        [D(i),U(i)]=SamWet(s,ds,us,cs,ul,ur,dl,dr);
    end
end

function d_guess=starte(ul,ur,dl,dr)
    %to provide starting value for Newton iteration. The Two-Rarefaction
    %Riemann solver (TRRS) and Two-Shock Riemann Solver (TSRS) are used
    %adaptively

    global cl cr;
    global g;

    dmin=min(dl,dr);
    %--------------------use TRRS solution as staring value--------------------
    d_guess=1/g*(0.5*(cl+cr)-0.25*(ur-ul))^2;
    if d_guess<=dmin
        %TRRS is suitable
    else
        %use TSRS solution as starting value with DS as computed from TRRS as
        %estimate
        gel=(0.5*g*(d_guess+dl)/(d_guess*dl))^0.5;
        ger=(0.5*g*(d_guess+dr)/(d_guess*dr))^0.5;
        d_guess=(gel*dl+ger*dr-(ur-ul))/(gel+ger);
    end
end

function [D,U]=SamLef(s,ul,ur,dl,dr)
    %to sample the solution for the case in which the left state is dry

    global cr;
    global g;

    shr=ur+cr;       %speed of right rarefaction head
    if s>=shr
        D=dr;
        U=ur;
    else
        str=ur-2*cr;         %speed of right rarefaction tail
        if s>=str
            %sampling point lies inside the rarefaction
            U=(ur-2*cr+2*s)/3;
            C=(-ur+2*cr+s)/3;
            D=C*C/g;
        else
            %sampling point lies in dry-bed state
            D=dl;
            U=ul;
        end
    end
end

function [D,U]=SamMid(s,ul,ur,dl,dr)
    %to sample the solution for the case in which the middle state is dry

    global cl cr;
    global g;

    %compute wave speeds
    shl=ul-cl;
    ssl=ul+2*cl;
    ssr=ur-2*cr;
    shr=ur+cr;

    if s<=shl
        %sampling point lies to the left of the left rarefaction
        D=dl;
        U=ul;
    end

    if (s>shl)&&(s<=ssl)
        %sampling point lies inside the left rarefaction
        U=(ul+2*cl+2*s)/3;
        C=(ul+2*cl-s)/3;
        D=C*C/g;
    end

    if (s>ssl)&&(s<=ssr)
        %sampling point lies inside the middle dry bed regions
        D=0;
        U=0;
    end

    if (s>ssr)&&(s<=shr)
        %sampling point lies inside the right rarefaction
        U=(ur-2*cr+2*s)/3;
        C=(-ur+2*cr+s)/3;
        D=C*C/g;
    end

    if s>shr
        %sampling point lies to the right of the right rarefaction
        D=dr;
        U=ur;
    end
end

function [D,U]=SamRig(s,ul,ur,dl,dr)
    %to sample the solution for the case in which the right state is dry

    global cl;
    global g;

    shl=ul-cl;    %speed of left rarefaction head

    if s<=shl 
        %sampling point lies to the left of the rarefaction
        D=dl;
        U=ul;
    else
        stl=ul+2*cl;
        if s<=stl
            %sampling point lies inside the rarefaction
            U=(ul+2*cl+2*s)/3;
            C=(ul+2*cl-s)/3;
            D=C*C/g;
        else
            %sample point lies in right dry-bed state
            D=dr;
            U=ur;
        end
    end
end

function [D,U]=SamWet(s,ds,us,cs,ul,ur,dl,dr)
    %to sample solution through wave structure at TimeOut for wet-bed case

    global cl cr;
    global g;

    if s<=us
        %---------------------sample left wave------------------------
        if ds>=dl
            %-----------left shock----------
            ql=((ds+dl)*ds/(2*dl*dl))^0.5;
            sl=ul-cl*ql;           %shock position
            if s<=sl
                D=dl;
                U=ul;
            else
                D=ds;
                U=us;
            end
        else
            %-----------left rarefaction----------
            shl=ul-cl;            %position of the rarefaction head
            if s<=shl
                %sample point lies to the left of the rarefaction
                D=dl;
                U=ul;
            else
                stl=us-cs;      %position of the rarafaction tail
                if s<=stl
                    %sample point lies inside the rarafaction
                    U=(ul+2*cl+2*s)/3;
                    C=(ul+2*cl-s)/3;
                    D=C*C/g;
                else
                    %sample point lies in the star region
                    D=ds;
                    U=us;
                end
            end
        end
    else
        %---------------------sample right wave------------------------
        if ds>dr
            %-----------right shock----------
            qr=((ds+dr)*ds/(2*dr*dr))^0.5;
            sr=ur+cr*qr;
            if (s>sr)
                D=dr;
                U=ur;
            else
                D=ds;
                U=us;
            end
        else
            %-----------right rarefaction----------
            shr=ur+cr;
            if s>=shr
                %sample point lies to the right of the rarefaction
                D=dr;
                U=ur;
            else
                str=us+cs;
                if s>=str
                    %sample point lies inside the rarefaction
                    U=(ur-2*cr+2*s)/3;
                    C=(-ur+2*cr+s)/3;
                    D=C*C/g;
                else
                    %sample point lies in the star region
                    D=ds;
                    U=us;
                end
            end
        end
    end
end

function [F,FD]=GEOFUN(D,dk,ck)
    %to evaluate function FL, FR and their derivatives in iterative Riemann
    %solver, for wet-bed case

    global g;

    if D<=dk
        %-------rarefaction wave-------
        C=(g*D)^0.5;
        F=2*(C-ck);
        FD=g/C;
    else
        %--------shock wave-------
        ges=(0.5*g*(D+dk)/(D*dk))^0.5;
        F=(D-dk)*ges;
        FD=ges-0.25*g*(D-dk)/(ges*D*D);
    end
end
