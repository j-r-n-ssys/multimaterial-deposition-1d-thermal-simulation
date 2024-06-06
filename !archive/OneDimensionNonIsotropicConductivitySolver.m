clc; clear

Model.Name = 'F375 Sinterable Model Material';
Model.Structure = 'Semi-Crystalline, Filled';
Model.BaseResin = 'POM';
Model.Density = .72*8030+(1-.72)*1540;
Model.ThermalConductivity = 10.6;
Model.SpecficHeatCapacity = 943;
Model.VolumetricHeatCapacity = Model.Density * Model.SpecficHeatCapacity;
Model.ThermalEffusivity = sqrt(Model.ThermalConductivity * Model.VolumetricHeatCapacity);
Model.Temperature.Glass = -30;
Model.Temperature.Melt = 170;
Model.Temperature.Extrusion = 235;
Model.ThermalDiffusivity = Model.ThermalConductivity / Model.VolumetricHeatCapacity;

Support.Name = 'QSR Support Material';
Support.Structure = 'Amorphous, Unfilled';
Support.BaseResin = 'PMMA (Approximation';
Support.Density = 1180;
Support.ThermalConductivity = 0.1;
Support.SpecficHeatCapacity = 1950;
Support.VolumetricHeatCapacity = Support.Density * Support.SpecficHeatCapacity;
Support.ThermalEffusivity = sqrt(Support.ThermalConductivity * Support.VolumetricHeatCapacity);
Support.Temperature.Glass = 165;
Support.Temperature.Melt = NaN;
Support.Temperature.Extrusion = 295;
Support.ThermalDiffusivity = Support.ThermalConductivity / Support.VolumetricHeatCapacity;

Process.Slice.Height = convlength(0.010,'in','m');
Process.Slice.Width = convlength(0.016,'in','m');
Process.Temperature.Envirionment = 115;
Process.ConvectionCoefficient = 200;

Simulation.LayerCount = 10;
Simulation.NodesPerLayer = 10;
Simulation.TotalNodes = Simulation.LayerCount * Simulation.NodesPerLayer;
Simulation.TimeStep = 10/1E5;
Simulation.EndTime = 1E1;
Simulation.TotalTimeSteps = floor(Simulation.EndTime/Simulation.TimeStep);

Coefficient.Conduction = 0;
Coefficient.Conduction = 0;

Results = LocalSolver(Model,Model,Process,Simulation);
StudyResults.MM.Time = Results.Time;
StudyResults.MM.Temperature = Results.Temperature.Interface;
StudyResults.MM.AdhesionRatio = Results.AdhesionRatio;

Results = LocalSolver(Model,Support,Process,Simulation);
StudyResults.MS.Time = Results.Time;
StudyResults.MS.Temperature = Results.Temperature.Interface;
StudyResults.MS.AdhesionRatio = Results.AdhesionRatio;

Results = LocalSolver(Support,Model,Process,Simulation);
StudyResults.SM.Time = Results.Time;
StudyResults.SM.Temperature = Results.Temperature.Interface;
StudyResults.SM.AdhesionRatio = Results.AdhesionRatio;

Results = LocalSolver(Support,Support,Process,Simulation);
StudyResults.SS.Time = Results.Time;
StudyResults.SS.Temperature = Results.Temperature.Interface;
StudyResults.SS.AdhesionRatio = Results.AdhesionRatio;


semilogx(StudyResults.MM.Time,StudyResults.MM.Temperature,'LineWidth',2)
hold on
semilogx(StudyResults.MS.Time,StudyResults.MS.Temperature,'LineWidth',2)
semilogx(StudyResults.SM.Time,StudyResults.SM.Temperature,'LineWidth',2)
semilogx(StudyResults.SS.Time,StudyResults.SS.Temperature,'LineWidth',2)
hold off
legend(sprintf('Model-Model (%3.1f%%)',StudyResults.MM.AdhesionRatio*100),...
   sprintf('Model-Support (%3.1f%%)',StudyResults.MS.AdhesionRatio*100),...
   sprintf('Support-Model (%3.1f%%)',StudyResults.SM.AdhesionRatio*100),...
   sprintf('Support-Support (%3.1f%%)',StudyResults.SS.AdhesionRatio*100));
axis([1e-4 1e1 100 250])
xlabel('Time [s]')
ylabel('Temperature [deg-C]')
set(gcf,'color','white')

function [Results] = LocalSolver(Extrudate,Base,Process,Simulation)
    
    %calcualte distance between nodes in Z
    Simulation.NodeSpacing = Process.Slice.Height / Simulation.NodesPerLayer;
   
    %pre-allocate explicity solution matrix
    ExplicitConductionMatrix = zeros(Simulation.TotalNodes);
    
    %pre-allocate explicit convection solution matrix
    ExplicitConvectionMatrix = zeros(1,Simulation.TotalNodes);
                     
    %pre-allocate temperature results matrix
    Results.Temperature.All = zeros(Simulation.TotalTimeSteps,Simulation.TotalNodes);
    
    for index=1:Simulation.TotalNodes

        if(index==1) % Extrudate (first node)            
                        
            ExplicitConductionMatrix(1:3,1)=Extrudate.ThermalDiffusivity*[1 -2 1]/Simulation.NodeSpacing^2;
            ExplicitConvectionMatrix(index) = (Process.ConvectionCoefficient)/(Extrudate.VolumetricHeatCapacity*Process.Slice.Height);
            Results.Temperature.All(1,index) = Extrudate.Temperature.Extrusion;
            
        elseif(index<Simulation.NodesPerLayer+1 && not(index==1)) % Extrudate    
            
            ExplicitConductionMatrix(index-1:index+1,index)=Extrudate.ThermalDiffusivity*[1 -2 1]/Simulation.NodeSpacing^2;            
            ExplicitConvectionMatrix(index) = (Process.ConvectionCoefficient)/(Extrudate.VolumetricHeatCapacity*Process.Slice.Height);
            Results.Temperature.All(1,index) = Extrudate.Temperature.Extrusion;
            
        elseif(index==Simulation.NodesPerLayer+1) % Interface
            
            ExplicitConductionMatrix(index-2:index+2,index)=InterfaceConductivity(Simulation.NodeSpacing,Extrudate.ThermalConductivity,Base.ThermalConductivity,Extrudate.ThermalDiffusivity,Base.ThermalDiffusivity);
            %Results.Temperature.All(1,index) = Process.Temperature.Envirionment+(Extrudate.Temperature.Extrusion-Process.Temperature.Envirionment)*(Base.ThermalEffusivity/(Base.ThermalEffusivity+Extrudate.ThermalEffusivity));
            ExplicitConvectionMatrix(index) = (Process.ConvectionCoefficient)/(mean([Extrudate.VolumetricHeatCapacity,Base.VolumetricHeatCapacity])*Process.Slice.Height);
            Results.Temperature.All(1,index) = (Process.Temperature.Envirionment+Extrudate.Temperature.Extrusion)/2;

        elseif(index>Simulation.NodesPerLayer+1 && not(index==Simulation.TotalNodes)) % Base

            ExplicitConductionMatrix(index-1:index+1,index)=Base.ThermalDiffusivity*[1 -2 1]/Simulation.NodeSpacing^2;
            ExplicitConvectionMatrix(index) = (Process.ConvectionCoefficient)/(Base.VolumetricHeatCapacity*Process.Slice.Height);
            Results.Temperature.All(1,index) = Process.Temperature.Envirionment;
            
        elseif(index==Simulation.TotalNodes) % Base (last node)
            
            ExplicitConductionMatrix(end-2:end,index)=Base.ThermalDiffusivity*[1 -2 1]/Simulation.NodeSpacing^2;            
            ExplicitConvectionMatrix(index) = (Process.ConvectionCoefficient)/(Base.VolumetricHeatCapacity*Process.Slice.Height);        
            Results.Temperature.All(1,index) = Process.Temperature.Envirionment;
            
        end

    end
        
    Results.Time = zeros(1,Simulation.TotalNodes);

    for index = 2:Simulation.TotalTimeSteps
        
        Results.Time(index) = (index-1)*Simulation.TimeStep;
        
        Results.Temperature.All(index,:) = Results.Temperature.All(index-1,:) + ...
                                        Simulation.TimeStep*(ExplicitConductionMatrix'*Results.Temperature.All(index-1,:)')' - ...
                                        Simulation.TimeStep*ExplicitConvectionMatrix.*(Results.Temperature.All(index-1,:)-Process.Temperature.Envirionment);
                                    
    
    end
    
    Results.Temperature.Interface = Results.Temperature.All(:,Simulation.NodesPerLayer+1);

    switch(Extrudate.Name)
        case('F375 Sinterable Model Material')
            aT=0.013*10.^(2541*(1./(Results.Temperature.Interface+273.14)-1/(195+273.14)));
        case('QSR Support Material')
            aT=0.217*10.^(-5.78*(Results.Temperature.Interface-200)./(182+Results.Temperature.Interface-200));
        otherwise
            error('Extrudate name not recognized')
    end
    

    Z=1./(2*pi*aT/1); 

    Z(Results.Temperature.Interface<max(Extrudate.Temperature.Glass,Extrudate.Temperature.Melt)-10)=0;

    Results.AdhesionRatio = abs(trapz(Z,Results.Time))^0.25;
    
    Results.AdhesionRatio = min(Results.AdhesionRatio,1);

    clear Z aT

end

function Array = InterfaceConductivity(dZ,K1,K2,D1,D2)
            
    Array =(1/(6*(K1+K2)*dZ^2))*...
        horzcat((2*K1+3*K2)*D1-D2*K1,...
        4*D2*K1-2*(K1+3*K2)*D1,...
        0,...
        4*D1*K2-2*(3*K1+K2)*D2,...
        (3*K1+2*K2)*D2-D1*K2);
    
end
