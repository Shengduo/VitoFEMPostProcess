function generateTrainingVN_t_function(totalprefix, target_x, z_locations, theta_flag)
    Time = linspace(0, 100, 1000) * 1.e-6;
    Vsave = zeros(length(totalprefix) * length(target_x), length(Time));
    Nsave = Vsave;
    
    start_idx = 0;
    
    for i = 1:1:length(totalprefix)
        disp(totalprefix)
        videoprefix = totalprefix(i); 
        faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');

        % Read time
        time = h5read(faultFileName, '/time');
        time = reshape(time, [1, size(time, 3)]);
        nOfTimeSteps = size(time, 2);

        % Read node geometry
        nodalXYZ = h5read(faultFileName, '/geometry/vertices');
        nOfNodes = size(nodalXYZ, 2);
        
        % Read nodal slip, slip rate and traction
        slip = h5read(faultFileName, '/vertex_fields/slip');
        SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
        traction = h5read(faultFileName, '/vertex_fields/traction');
        if theta_flag == 1
            theta = h5read(faultFileName, '/vertex_fields/state_variable');
        end
        
        connection = h5read(faultFileName, '/topology/cells');
        connection = connection + 1;

        % Read nodal slip
        % Input the first wire position
        WirePos1 = [-0.025657; -0.014222; 0];
        VSstart = [0.006354, 0.003522, 0]';
        
        % Calculate NodalXYZ2D 
        NodalXYZ2D = zeros(2, size(nodalXYZ, 2));
        for ii = 1:1:size(nodalXYZ, 2)
            NodalXYZ2D(1, ii) = sign(nodalXYZ(1, ii) - VSstart(1)) * norm(nodalXYZ(1:2, ii) - VSstart(1:2), 2);
            NodalXYZ2D(2, ii) = nodalXYZ(3, ii);
        end

        % Fault Range
        xrange = [-100, 150] - 1e3 * norm(VSstart - WirePos1);
        
        % Magnitude of slip
        slipMag = zeros(nOfNodes, nOfTimeSteps);
        
        % Magnitude of slip rate
        slipRateMag = zeros(nOfNodes, nOfTimeSteps);

        for t = 1:1:nOfTimeSteps
            for ii = 1:1:nOfNodes
                slipRateMag(ii, t) = norm(SlipRate(:, ii, t), 2);
                slipMag(ii, t) = norm(slip(:, ii, t), 2);
            end
        end
        
        % Loop through z locations
        for zz = 1:1:length(z_locations)
            z_location = z_locations(zz);
            
            % Get the mesh stored as triangulation, get the location of places
            % of interest
            TR = triangulation(connection', NodalXYZ2D');
            QueryXYZ = [target_x / 1e3; z_location / 1e3 * ones(1, size(target_x, 2))]';
            elementID = pointLocation(TR, QueryXYZ);
            nOfQueryPts = size(target_x, 2);

            %% Get sliprate, shear stress aned normal stress history at target points
            QuerySlip = zeros(nOfQueryPts, nOfTimeSteps);
            QueryV = zeros(nOfQueryPts, nOfTimeSteps);
            QueryShearStress = ones(nOfQueryPts, nOfTimeSteps);
            QueryNormalStress = ones(nOfQueryPts, nOfTimeSteps);
            if theta_flag == 1
                QueryTheta = ones(nOfQueryPts, nOfTimeSteps);
            end

            % Get the interpolation
            for t = 1:1:nOfTimeSteps
                for pt = 1:1:nOfQueryPts
                    % Slip
                    slip = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(slipMag(connection(:, elementID(pt)), t)), 'natural');
                    QuerySlip(pt, t) = slip(QueryXYZ(pt, 1), QueryXYZ(pt, 2));

                    % Slip rate
                    vel_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(slipRateMag(connection(:, elementID(pt)), t)), 'natural');
                    QueryV(pt, t) = vel_x(QueryXYZ(pt, 1), QueryXYZ(pt, 2));

                    % Shear stress
                    shear = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(traction(1, connection(:, elementID(pt)), t))', 'natural');
                    QueryShearStress(pt, t) = shear(QueryXYZ(pt, 1), QueryXYZ(pt, 2));

                    % Normal stress
                    normal = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(traction(3, connection(:, elementID(pt)), t))', 'natural');
                    QueryNormalStress(pt, t) = normal(QueryXYZ(pt, 1), QueryXYZ(pt, 2));

                    if theta_flag == 1
                        % State variable
                        stateVar = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(theta(1, connection(:, elementID(pt)), t))', 'natural');
                        QueryTheta(pt, t) = stateVar(QueryXYZ(pt, 1), QueryXYZ(pt, 2));
                    end
                end
            end

            % Test and Train data set
            Vsave(start_idx + 1 : start_idx + nOfQueryPts, :) = (interp1(time, QueryV', Time))';
            Nsave(start_idx + 1 : start_idx + nOfQueryPts, :) = -(interp1(time, QueryNormalStress', Time))' / 1e6;
            start_idx = start_idx + nOfQueryPts;
        end
    end
    
    % Save Nsave Vsave
    save('~/InverseProblems/RateStateFrictionTraining/data/realDatas1018.mat', 'Nsave', 'Vsave', 'Time');
end
