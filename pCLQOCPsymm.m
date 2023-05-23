classdef pCLQOCPsymm < pCLQOCP
%PCLQOCPSYMM represents symmetric constrained linear-quadratic optimal
%   control problems (OCPs).
%
% Syntax:
%   myQP = PCLQOCPSYMM(A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N,allTheta,allOmega)
%   myQP = PCLQOCPSYMM(myLTISystem,N,allTheta,allOmega)
%
%
% Description:
%
%   myQP = pCLQOCP(A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N,allTheta,allOmega) creates
%       an object for a constrained linear-quatratic OCP of the form
%
%       min(U,X) x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
%       s.t.     x(k+1) = A*x(k) + B*u(k), k=0,...,N-1
%                Ou*u(k) <= cu, k=0,...,N-1
%                Ox*x(k) <= cx, k=0,...,N-1
%                Ot*x(N) <= ct.
%
%       For more information on entering the terminal set Ot*x(N) <= ct and
%       the terminal penalty P, see documentation of pCLQOCP.
%
%       Symmetries of the object are defined as symmetric pairs (Theta,Omega)
%       (see Danielson2015*). All matrices Theta and Omega are defined in 
%       the inputs allTheta and allOmega, respectively, such that their 
%       ith sheets define the matrices Theta and Omega of the ith symmetric pair.
%
%   myQP = pCLQOCP(myLTISystem,N,allTheta,allOmega) creates an object for 
%       an object of class LTISystem from the multi-parametric toolbox** and
%       a given horizon N,
%
%       min(U,X) x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
%       s.t.     x(k+1) = A*x(k) + B*u(k), k=0,...,N-1
%                umin <= u(k) <= umax, k=0,...,N-1
%                xmin <= x(k) <= xmax, k=0,...,N-1
%                Ot*x(N) <= ct.
%
%       For more information on entering the terminal set Ot*x(N) <= ct and
%       the terminal penalty P, see documentation of pCLQOCP.
%
%       Symmetries of the object are defined as symmetric pairs (Theta,Omega)
%       (see Danielson2015*). All matrices Theta and Omega are defined in 
%       the inputs allTheta and allOmega, respectively, such that their 
%       ith sheets define the matrices Theta and Omega of the ith symmetric pair.
%
%
% Examples:
%
%   create object similar to Example 2 from Danielson2015*
%      A        = [1 1;-1 1];
%      B        = eye(2);
%      Q        = eye(2);
%      R        = eye(2);
%      N        = 3; 
%      Ox       = [0 1;0 -1;1 0; -1 0; 1/sqrt(2)*[1 1; -1 1;1 -1;-1 -1]];
%      cx       = 10*ones(8,1);
%      Ou       = [0 1;0 -1;1 0; -1 0; 1/sqrt(2)*[1 1; -1 1;1 -1;-1 -1]];
%      cu       = ones(8,1);
%      allTheta = zeros(2,2,8); allTheta(:,:,1) = eye(2); allTheta(:,:,2) = -eye(2); allTheta(:,:,3) = [0 -1; 1 0]; allTheta(:,:,4) = [0 1; -1 0]; allTheta(:,:,5) = 1/sqrt(2)*[1 -1; 1 1]; allTheta(:,:,6) = 1/sqrt(2)*[-1 -1; 1 -1]; allTheta(:,:,7) = 1/sqrt(2)*[-1 1; -1 -1]; allTheta(:,:,8) = 1/sqrt(2)*[1 1; -1 1];
%      allOmega = zeros(2,2,8); allOmega(:,:,1) = eye(2); allOmega(:,:,2) = -eye(2); allOmega(:,:,3) = [0 -1; 1 0]; allOmega(:,:,4) = [0 1; -1 0]; allOmega(:,:,5) = 1/sqrt(2)*[1 -1; 1 1]; allOmega(:,:,6) = 1/sqrt(2)*[-1 -1; 1 -1]; allOmega(:,:,7) = 1/sqrt(2)*[-1 1; -1 -1]; allOmega(:,:,8) = 1/sqrt(2)*[1 1; -1 1];
%      myQP     = pCLQOCPsymm(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N,allTheta,allOmega);
%
%   create object from Example 1 from Danielson2012***
%      myLTISystem           = LTISystem('A',[2 1;-1 2],'B',eye(2));
%      myLTISystem.x.min     = [-1; -1];
%      myLTISystem.x.max     = [1; 1];
%      myLTISystem.u.min     = [-1; -1];
%      myLTISystem.u.max     = [1; 1];
%      myLTISystem.x.penalty = QuadFunction(eye(2));
%      myLTISystem.u.penalty = QuadFunction(5000*eye(2));
%      N = 3;
%      allTheta = zeros(2,2,4); allTheta(:,:,1) = eye(2); allTheta(:,:,2) = -eye(2); allTheta(:,:,3) = [0 -1; 1 0]; allTheta(:,:,4) = [0 1; -1 0];
%      allOmega = zeros(2,2,4); allOmega(:,:,1) = eye(2); allOmega(:,:,2) = -eye(2); allOmega(:,:,3) = [0 -1; 1 0]; allOmega(:,:,4) = [0 1; -1 0];
%      myQP     = pCLQOCPsymm(myLTISystem,N,allTheta,allOmega);
%
%
% *C. R. Danielson and F. Borrelli, "Symmetric Linear Model Predictive 
%  Control, " IEEE Transactions on Automatic Control," vol. 60(5), pp. 
%  1244-1259, 2015.
%
% **M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari, "Multi-
%  Parametric Toolbox 3.0," Proc. of the European Control
%  Conference, pp. 502-510, 2013.
%
% ***C. R. Danielson and F. Borrelli, "Symmetric explicit model predictive
%  control," IFAC Proceedings Volumes, vol. 45(17), pp. 132-137, 2012.
%
% 
% See also LTISystem and pCLQOCP and pQP

%  Authors
%    
%   2022    Ruth Mitze and Martin Mönnigmann:
%           Ruhr-Universität Bochum
%           Systems Theory and Automatic Control
%   mailto: ruth.mitze@rub.de 
%   mailto: martin.moennigmann@rub.de

%  License
%    
%    This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation; either version 3 of the License, or (at
%  your option) any later version.
%    This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
%  General Public License for more details.
%    You should have received a copy of the GNU Lesser General Public
%  License along with this library; if not, write to the  Free Software
%  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%  MA 02110-1335 USA
    
    properties (SetAccess = private)
        allTheta           % each sheet allTheta(:,:,i) represents matrix Theta of ith symmetric pair (Theta,Omega)
        allOmega           % each sheet allOmega(:,:,i) represents matrix Omega of ith symmetric pair (Theta,Omega)
        
        solutionReduced    % set of only one active set per orbit of solution
        nSymm              % number of symmetries
        allTransformations % each sheet allTranformations(:,:,i) represents the constraint transformation matrix corresponding to the ith symmetric pair (Theta,Omega)
    end
    
    methods        
        function obj = pCLQOCPsymm(varargin)
            %PCLQOCPsymm constructs an instance of class pCLQOCPsymm
            %
            % 
            % See also pCLQOCPsymm
            
            % case invalid number of inputs
            if nargin~=4 && nargin~=14
                error('Inputs musst follow either (myLTISystem,N,allTheta,allOmega) or (A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N,allTheta,allOmega)');
            end
                    
            obj = obj@pCLQOCP(varargin{1:end-2});
            
            % test if input allTheta is correct
            if ~isnumeric(varargin{end-1}) || ~isreal(varargin{end-1}) || size(varargin{end-1},1)~=obj.n || size(varargin{end-1},2)~=obj.n
                error('Input allTheta must be a threedimensional numerical matrix with the same number of lines and columns as A')
            end
            obj.allTheta = varargin{end-1};
            obj.nSymm = size(obj.allTheta,3);
            
            % test if input allOmega is correct
            if ~isnumeric(varargin{end}) || ~isreal(varargin{end}) || size(varargin{end},3)~=obj.nSymm || size(varargin{end},1)~=obj.m || size(varargin{end},2)~=obj.m
                error('Input allOmega must be a threedimensional numerical matrix with the same number of lines and columns as the number of columns in B, and the same number of sheets as allTheta.')
            end
            obj.allOmega = varargin{end};
        end
        
        
        function obj = solveDP(obj)
            %SOLVEDP determines the property solution for problems of class
            %   pCLOCPsymm using dynamic programming (first introduced in
            %   Mitze2020*) with the extension for symmetric problems
            %   proposed in Mitze2022**. The algorithm also implements the
            %   cardinality limit for processed active sets introduced in 
            %   Mitze2022***.
            %
            %   A necessary condition for using this algorithm is that the
            %   problem considers the LQR penalty as terminal penalty and
            %   the LQR invariant set as terminal set.
            % 
            % Syntax:
            %   myQP = myQP.SOLVEDP
            %
            %
            % Examples:
            %
            %   create object from Example 1 from Danielson2012****
            %      A        = [2 1;-1 2];
            %      B        = eye(2);
            %      Q        = eye(2);
            %      R        = 5000*eye(2);
            %      N        = 3; 
            %      Ox       = [-eye(2); eye(2)];
            %      cx       = ones(4,1);
            %      Ou       = [-eye(2); eye(2)];
            %      cu       = ones(4,1);
            %      allTheta = zeros(2,2,4); allTheta(:,:,1) = eye(2); allTheta(:,:,2) = -eye(2); allTheta(:,:,3) = [0 -1; 1 0]; allTheta(:,:,4) = [0 1; -1 0];
            %      allOmega = zeros(2,2,4); allOmega(:,:,1) = eye(2); allOmega(:,:,2) = -eye(2); allOmega(:,:,3) = [0 -1; 1 0]; allOmega(:,:,4) = [0 1; -1 0];
            %      myQP     = pCLQOCPsymm(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N,allTheta,allOmega);
            %
            %   compute solution
            %     myQP = myQP.solveDP;
            %     myQP.solution
            %
            %
            % *R. Mitze and M. Mönnigmann, "A dynamic programming
            %  approach to solving constrained linear-quadratic optimal
            %  control problems", Automatica, vol. 120, p. 109132, 2020.
            %
            % **R. Mitze, M. Kvasnica, and M. Mönnigmann, "Exploiting
            %  symmetries in active set enumeration for constrained
            %  linear-quadratic optimal control", Automatica, vol. 151,
            %  p. 110900, 2023.
            %
            % ***R. Mitze and M. Mönnigmann, "Improved active set dynamic
            %  programming for solving linear-quadratic optimal control
            %  problems", 2022 IEEE 61st Conference on Decision and Control
            %  (CDC), pp. 1764-1769, 2022.
            %
            % ****C. R. Danielson and F. Borrelli, "Symmetric explicit 
            %  model predictive control," IFAC Proceedings Volumes, vol. 
            %  45(17), pp. 132-137, 2012.
            
            % test if terminal set and terminal penalty are defined
            if obj.qT==0 || ~any(any(obj.P))
                error('The terminal penalty of the system must be the LQR penalty and the terminal set must be the LQR invariant set.')
            else
                fprintf('\nNote: the terminal penalty of the system must be the LQR penalty and the terminal set must be the LQR invariant set\n')
            end
            
            % save target horizon
            N_max = obj.N;
            
            % save maximal cardinality of processed active sets
            maxCardinality = size(obj.G,2);

            % determine reduced solution for initial horizon
            obj.N = 1;
            fprintf('\ndetermine primary solution for initial horizon N=1\n')
            [obj,solution_N,degenerate_N] = obj.solveCombinatorialConsiderDegeneracy(maxCardinality);

            % increment horizon
            while obj.N+1<=N_max
                obj.N          = obj.N+1;
                solution_Nm1   = solution_N;
                degenerate_Nm1 = degenerate_N;
                
                % determine solution for incremented horizon
                fprintf(['determine solution for increased horizon N=' num2str(obj.N) '\n'])
                [obj,solution_N,degenerate_N] = obj.solveWithPreviousSolution(solution_Nm1,degenerate_Nm1,maxCardinality);

                % test if solution of infinite horizon is reached
                if isempty(find(solution_N(:,(obj.N-1)*obj.qUX+1:end),1))
                    fprintf(['\nfeasible set for N_max=' num2str(N_max) ' was reached for N=' num2str(obj.N-1) '\n'])
                    
                    % redo QP matrices for horizon N_max and add inactive constraints to solution
                    solution_N = [solution_N false(size(solution_N,1),(N_max-obj.N)*obj.qUX)];
                    obj.N = N_max;
                    break
                end
            end            
            solution   = solution_N;
            degenerate = degenerate_N;

            delete = false(size(solution,1),1);
            % for every active set in the solution
            for activeSet_i=1:size(solution,1)     
                
                % test if obj.G(active set,:) has full row rank
                if sum(solution(activeSet_i,:))~=rank(obj.G(solution(activeSet_i,:),:))
                    
                    % if not satisfied, remove from solution
                    delete(activeSet_i,1) = true;
                else
                    % if satisfied and if degenerate (optimality LP resulted t=0)
                    if degenerate(activeSet_i,:)
                        % compute polytope
                        [T,d] = computePolytope(obj,solution(activeSet_i,:));
                        if ~isFullDim(Polyhedron(T,d))
                            % remove from solution if polytope is lower dimensional
                            delete(activeSet_i,1) = true;
                        end
                    end
                end
            end
            solution(delete,:) = [];
            
            % store current solution as reduced solution
            obj.solutionReduced = solution;
            
            % complete solution by adding all active sets of the corresponding orbits
            complete_solution = false(size(solution,1)*obj.nSymm,obj.q);
            obj = determineTransformationMatrix(obj);
            counter_temp = 0;
            for i=1:size(solution,1)
                [~,orbit] = determineOrbit(obj,solution(i,:));
                counter_temp = counter_temp+size(orbit,1);
                complete_solution(counter_temp-size(orbit,1)+1:counter_temp,:) = orbit;
            end
            complete_solution(counter_temp+1:end,:)=[];            
            obj.solution = complete_solution;
        end
          
        function obj = solveCombinatorial(obj)
            %SOLVECOMBINATORIAL determines the property solution with
            %   combinatorial mpQP (algorithm from Gupta2011*). The
            %   algorithm is specially tailored for symmetric problems with
            %   the improvement proposed in Mitze2023**.
            % 
            % Syntax:
            %   myQP = myQP.SOLVECOMBINATORIAL
            %
            %
            % Examples:
            %
            %   create object from Example 1 from Danielson2012***
            %      A        = [2 1;-1 2];
            %      B        = eye(2);
            %      Q        = eye(2);
            %      R        = 5000*eye(2);
            %      N        = 3; 
            %      Ox       = [-eye(2); eye(2)];
            %      cx       = ones(4,1);
            %      Ou       = [-eye(2); eye(2)];
            %      cu       = ones(4,1);
            %      allTheta = zeros(2,2,4); allTheta(:,:,1) = eye(2); allTheta(:,:,2) = -eye(2); allTheta(:,:,3) = [0 -1; 1 0]; allTheta(:,:,4) = [0 1; -1 0];
            %      allOmega = zeros(2,2,4); allOmega(:,:,1) = eye(2); allOmega(:,:,2) = -eye(2); allOmega(:,:,3) = [0 -1; 1 0]; allOmega(:,:,4) = [0 1; -1 0];
            %      myQP     = pCLQOCPsymm(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N,allTheta,allOmega);
            %
            %   compute solution
            %     myQP = myQP.solveCombinatorial;
            %     myQP.solution
            %
            %
            % *A. Gupta, S. Bhartiya, and P. Nataraj, "A novel approach to
            %  multiparametric quadratic programming," Automatica, vol. 47,
            %  pp. 2112-2117, 2011.
            %
            % **R. Mitze, and M. Mönnigmann, "Exploiting symmetries in 
            %  tree-based combinatorial calculation of explicit linear 
            %  MPC solutions symmetries in active set enumeration for 
            %  constrained linear-quadratic optimal control", accepted to 
            %  the 24th International Conference on Process Control (PC),
            %  2023.
            %
            % ***C. R. Danielson and F. Borrelli, "Symmetric explicit model
            %  predictive control," IFAC Proceedings Volumes, vor. 45(17), 
            %  pp. 132-137, 2012.
            
            fprintf('\ndetermine solution with combinatorial mpQP for symmetrical problems\n')

            % save maximal cardinality of processed active sets
            maxCardinality = min(size(obj.G,2),obj.q);

            % initialize
            chunkSize = 1000; % number of preallocated lines (overestimated guess)
            solution             = false(chunkSize,obj.q);
            infeasibleActiveSets = false(chunkSize,obj.q);
            nonprimaryActiveSets = false(chunkSize,obj.q);
            counter_solution             = 0;
            counter_infeasibleActiveSets = 0;
            counter_nonprimaryActiveSets = 0;

            % determine constraint transformation matrices
            obj = determineTransformationMatrix(obj);
            
            % for every active set in combinatorial tree
            candidateActiveSet = false(1,obj.q);
            while true
                
                % if obj.G(active set,:) has full row rank
                if sum(candidateActiveSet)==rank(obj.G(candidateActiveSet,:))
                        
                    % if no subset is a detected infeasible active set
                    [obj,exitflag] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),candidateActiveSet);
                    if ~exitflag
                        
                        % if active set is not a childnode of detected non-primary active set
                        nonprimary = false;
                        if any(candidateActiveSet)  
                            [obj,nonprimary] = isNonprimary(obj,nonprimaryActiveSets(1:counter_nonprimaryActiveSets,:),candidateActiveSet);
                        end
                        if ~nonprimary
                            
                            % solve optimality LP
                            [obj,t,exitflagOpt] = isOptimal(obj,candidateActiveSet);
                            if exitflagOpt==1
                                
                                % store if optimal and t>0
                                if t>eps
                                    counter_solution = counter_solution+1;
                                    if ~mod(counter_solution,chunkSize)
                                        solution = [solution; false(chunkSize,obj.q)]; % allocate more lines
                                    end
                                    solution(counter_solution,:) = candidateActiveSet;
                                
                                % compute polytope if optimal and t=0
                                else
                                    [T,d] = computePolytope(obj,candidateActiveSet);
                                    if isFullDim(Polyhedron(T,d))
                                        % store if polytope is full-dimensional
                                        counter_solution = counter_solution+1;
                                        if ~mod(counter_solution,chunkSize)
                                            solution = [solution; false(chunkSize,obj.q)]; % allocate more lines
                                        end
                                        solution(counter_solution,:) = candidateActiveSet;
                                    end
                                end
                                
                                % for every element in the orbit of the processed active set
                                [obj,orbit] = determineOrbit(obj,candidateActiveSet);
                                for numberSymmetries = 1:size(orbit,1)                        

                                    % if it is not a childnode of a detected non-primary active set
                                    if ~all(orbit(numberSymmetries,:)==candidateActiveSet)
                                        nonprimary = false;
                                        if any(orbit(numberSymmetries,:))  
                                            [obj,nonprimary] = isNonprimary(obj,nonprimaryActiveSets(1:counter_nonprimaryActiveSets,:),orbit(numberSymmetries,:));
                                        end
                                        if ~nonprimary
                                            % add to set of non-primary active sets
                                            counter_nonprimaryActiveSets = counter_nonprimaryActiveSets+1;
                                            if ~mod(counter_nonprimaryActiveSets,chunkSize)
                                                nonprimaryActiveSets = [nonprimaryActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                            end
                                            nonprimaryActiveSets(counter_nonprimaryActiveSets,:) = orbit(numberSymmetries,:);
                                        end
                                    end
                                end

                            % solve feasibility LP if not optimal
                            elseif exitflagOpt==-2
                                [obj,exitflagFeas] = isFeasible(obj,candidateActiveSet);

                                if exitflagFeas~=-2
                                    
                                    % for every element in the orbit of the processed active set
                                    [obj,orbit] = determineOrbit(obj,candidateActiveSet);
                                    for numberSymmetries = 1:size(orbit,1)                        

                                        % if it is not a childnode of a detected non-primary active set
                                        if ~all(orbit(numberSymmetries,:)==candidateActiveSet)
                                            nonprimary = false;
                                            if any(orbit(numberSymmetries,:))  
                                                [obj,nonprimary] = isNonprimary(obj,nonprimaryActiveSets(1:counter_nonprimaryActiveSets,:),orbit(numberSymmetries,:));
                                            end
                                            if ~nonprimary
                                                % add to set of non-primary active sets
                                                counter_nonprimaryActiveSets = counter_nonprimaryActiveSets+1;
                                                if ~mod(counter_nonprimaryActiveSets,chunkSize)
                                                    nonprimaryActiveSets = [nonprimaryActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                                end
                                                nonprimaryActiveSets(counter_nonprimaryActiveSets,:) = orbit(numberSymmetries,:);
                                            end
                                        end
                                    end
                                
                                % add to set of infeasible active sets if infeasible
                                else
                                    % add active set itself to pruned
                                    counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                                    if ~mod(counter_infeasibleActiveSets,chunkSize)
                                        infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                    end
                                    infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidateActiveSet;
                                    
                                    % add whole orbit to pruned if not already a superset
                                    % for every element in the orbit of the processed active set
                                    [obj,orbit] = determineOrbit(obj,candidateActiveSet);
                                    for numberSymmetries = 1:size(orbit,1)
                                        % if it is not a childnode of a detected non-primary active set
                                        if ~all(orbit(numberSymmetries,:)==candidateActiveSet)
                                            exitflag = false;
                                            if any(orbit(numberSymmetries,:))  
                                                [obj,exitflag] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),orbit(numberSymmetries,:));
                                            end
                                            if ~exitflag
                                                % add to set of infeasible active sets
                                                counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                                                if ~mod(counter_infeasibleActiveSets,chunkSize)
                                                    infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                                end
                                                infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidateActiveSet;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                % compute next active set in combinatorial tree
                [obj,candidateActiveSet,finished] = nextActiveSet(obj,candidateActiveSet,maxCardinality);
                if finished
                    break
                end
            end
            
            % store current solution as reduced solution
            obj.solutionReduced = solution(1:counter_solution,:);
            
            % complete solution by adding active sets of orbits
            complete_solution = false(counter_solution*obj.nSymm,obj.q);
            counter_temp = 0;
            for i=1:counter_solution
                [~,orbit] = determineOrbit(obj,solution(i,:));
                counter_temp = counter_temp+size(orbit,1);
                complete_solution(counter_temp-size(orbit,1)+1:counter_temp,:) = orbit;
            end
            complete_solution(counter_temp+1:end,:)=[];            
            obj.solution = complete_solution;
        end
        
        
        function testSymmetries(obj)
            %TESTSYMMETRIES tests if the symmetric pairs entered as inputs
            % allTheta and allOmega for object of class pCLQOCPsymm are
            % symmetric pairs of the problem (see Danielson2015*) and if the
            % symmetric pairs satisfy group properties.
            %
            %   Syntax:
            %       myQP.testSymmetries
            %
            %
            % Examples:
            %
            %   create object from Example 1 from Danielson2012**
            %      A        = [2 1;-1 2];
            %      B        = eye(2);
            %      Q        = eye(2);
            %      R        = 5000*eye(2);
            %      N        = 3; 
            %      Ox       = [-eye(2); eye(2)];
            %      cx       = ones(4,1);
            %      Ou       = [-eye(2); eye(2)];
            %      cu       = ones(4,1);
            %      allTheta = zeros(2,2,4); allTheta(:,:,1) = eye(2); allTheta(:,:,2) = -eye(2); allTheta(:,:,3) = [0 -1; 1 0]; allTheta(:,:,4) = [0 1; -1 0];
            %      allOmega = zeros(2,2,4); allOmega(:,:,1) = eye(2); allOmega(:,:,2) = -eye(2); allOmega(:,:,3) = [0 -1; 1 0]; allOmega(:,:,4) = [0 1; -1 0];
            %      myQP     = pCLQOCPsymm(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N,allTheta,allOmega);
            %
            %   test entered symmetric pairs
            %      myQP.testSymmetries
            %
            %
            % *C. R. Danielson and F. Borrelli, "Symmetric Linear Model Predictive 
            %  Control, " IEEE Transactions on Automatic Control," vol. 60(5), pp. 
            %  1244-1259, 2015
            %
            % **C. R. Danielson and F. Borrelli, "Symmetric explicit model predictive
            %  control," IFAC Proceedings Volumes, vor. 45(17), pp. 132-137, 2012
            %
            % ***C. R. Danielson, "Symmetric Constrained Optimal Control: Theory,
            %  Algorithms, and Applications," PhD thesis, University of California,
            %  Berkeley, 2014
            
            % initialize
            flag_global = true;
            
            % Step 1:test if set of symmetric pairs are symmetries of the problem
            % for every symmetric pair
            for i=1:obj.nSymm
                
                % test if symmetry satisfies the symmetry of the dynamics
                if ~all(all(obj.A*obj.allTheta(:,:,i)-obj.allTheta(:,:,i)*obj.A<=1e-6)) || ~all(all(obj.B*obj.allOmega(:,:,i)-obj.allTheta(:,:,i)*obj.B<=1e-6)) % !! tolerance manually chosen !!
                    flag_global = false;
                    fprintf(['\nViolation of symmetry of dynamics for symmetric pair number ' num2str(i) ' detected\n'])
                end
                
                % test if symmetry satisfies the symmetry of the cost (if terminal penalty is LQR penalty, this is always satisfied, see Prop. 14 in Danielson2014***)
                if ~all(all(obj.allTheta(:,:,i)'*obj.Q*obj.allTheta(:,:,i)-obj.Q<=1e-6)) || ~all(all(obj.allOmega(:,:,i)'*obj.R*obj.allOmega(:,:,i)-obj.R<=1e-6)) % !! tolerance manually chosen !!
                    flag_global = false;
                    fprintf(['\nViolation of symmetry of cost for symmetric pair number ' num2str(i) ' detected\n'])
                end
                
                % test if symmetry satisfies the symmetry of the constraints (if terminal set is LQR set, this is always satisfied, see Prop. 10 in Danielson2014***)
                for j=1:length(obj.cx)
                    flag = false;
                    for k=1:length(obj.cx)
                        if all(obj.Ox(j,:)/obj.cx(j)-obj.Ox(k,:)*obj.allTheta(:,:,i)/obj.cx(k)<=1e-6) % !! tolerance manually chosen !!
                            flag = true;
                            break
                        end
                    end
                    if ~flag
                        flag_global = false;
                        fprintf(['\nViolation of symmetry of constraints for symmetric pair number ' num2str(i) ' detected\n'])
                        break
                    end
                end
                if flag
                    for j=1:length(obj.cu)
                        flag = false;
                        for k=1:length(obj.cu)
                            if all(obj.Ou(j,:)/obj.cu(j)-obj.Ou(k,:)*obj.allOmega(:,:,i)/obj.cu(k)<=1e-6) % !! tolerance manually chosen !!
                                flag = true;
                                break
                            end
                        end
                        if ~flag
                            flag_global = false;
                            fprintf(['\nViolation of symmetry of constraints for symmetric pair number ' num2str(i) ' detected\n'])
                            break
                        end
                    end
                end
            end
            
            % Step 2:test if set of symmetric pairs is a group under matrix multiplication
            % test for the group property "identity"
            flag = false;
            for i=1:obj.nSymm
                if all(all(obj.allTheta(:,:,i)-eye(obj.n)<=1e-6)) && all(all(obj.allOmega(:,:,i)-eye(obj.n)<=1e-6)) % !! tolerance manually chosen !!
                    flag = true;
                    break
                end
            end
            if ~flag
                flag_global = false;
                fprintf('\nViolation of group properties detected: missing identity element Theta = eye(obj.n), Omega = eye(obj.m)\n')
            end
            
            % test for the group property "inverse"
            for i=1:obj.nSymm
                flag = false;
                ThetaInv = inv(obj.allTheta(:,:,i));
                OmegaInv = inv(obj.allOmega(:,:,i));
                for j=1:obj.nSymm
                    if all(all(ThetaInv-obj.allTheta(:,:,j)<=1e-6)) && all(all(OmegaInv-obj.allOmega(:,:,j)<=1e-6)) % !! tolerance manually chosen !!
                        flag = true;
                        break
                    end
                end
                if ~flag
                    flag_global = false;
                    fprintf(['\nViolation of group properties detected: no inverse for symmetric pair number ' num2str(i) '\n'])
                end
            end
            
            % test for the group property "closure"
            for i=1:obj.nSymm
                for j=1:obj.nSymm
                    flag = false;
                    ThetaProd = obj.allTheta(:,:,i)*obj.allTheta(:,:,j);
                    OmegaProd = obj.allOmega(:,:,i)*obj.allOmega(:,:,j);
                    for k=1:obj.nSymm
                        if all(all(ThetaProd-obj.allTheta(:,:,k)<=1e-6)) && all(all(OmegaProd-obj.allOmega(:,:,k)<=1e-6)) % !! tolerance manually chosen !!
                            flag = true;
                            break
                        end
                    end
                    if ~flag
                        flag_global = false;
                        fprintf(['\nViolation of group properties detected: no closure for combination of symmetric pairs with numbers' num2str(i) ' and ' num2str(j) '\n'])
                    end
                end
            end
            
            % case no violations detected
            if flag_global
                fprintf('\nNo violation of symmetries or group properties detected\n')
            end
        end
    end
    
    methods (Access = protected)        
        function [obj,solution,degenerate] = solveCombinatorialConsiderDegeneracy(obj,maxCardinality)
            %SOLVECOMBINATORIALCONSIDERDEGENERACY determines all primary optimal
            %   active sets with combinatorial mpQP (Alg. 1 from Mitze2022*,
            %   improved by strategy similar to strategy in Mitze2023***)
            %   and stores them in solution. Active sets such that
            %   obj.G(active set,:) is not of full row rank are marked (see
            %   also Mitze2020**).
            % 
            %
            % Input Arguments:
            %   maxCardinality: maximal cardinality of active sets.
            %
            % Output Arguments:
            %   solution: matrix containing all optimal active sets 
            %   (expressed as logical line vectors of length obj.q,
            %   inactive constraint = false, active constraint = true)
            %
            %   degenerate: logical vector. Every entry corresponds to a
            %   line in solution. An entry is true if the active set in
            %   that line has the solution t=0 to the optimality LP, and 
            %   false otherwise
            %
            %
            % *R. Mitze, M. Kvasnica, and M. Mönnigmann, "Exploiting
            %  symmetries in active set enumeration for constrained
            %  linear-quadratic optimal control", Automatica, vol. 151,
            %  p. 110900, 2023.
            %
            % **R. Mitze and M. Mönnigmann, "A dynamic programming
            %  approach to solving constrained linear-quadratic optimal
            %  control problems", Automatica, vol. 120, p. 109132, 2020.
            %
            % ***R. Mitze, and M. Mönnigmann, "Exploiting symmetries in 
            %  tree-based combinatorial calculation of explicit linear 
            %  MPC solutions symmetries in active set enumeration for 
            %  constrained linear-quadratic optimal control", accepted to 
            %  the 24th International Conference on Process Control (PC),
            %  2023.

            % save maximal cardinality for processed active sets
            maxCardinality = min(maxCardinality,obj.q);
            
            % initialize
            chunkSize = 1000; % number of preallocated lines (overestimated guess)
            solution             = false(chunkSize,obj.q);
            infeasibleActiveSets = false(chunkSize,obj.q);
            degenerate           = false(chunkSize,1);
            nonprimaryActiveSets = false(chunkSize,obj.q);
            counter_solution             = 0;
            counter_infeasibleActiveSets = 0;
            counter_nonprimaryActiveSets = 0;

            % determine constraint transformation matrices
            obj = determineTransformationMatrix(obj);
            
            % for every active set in the combinatorial tree
            candidateActiveSet = false(1,obj.q);
            while true

                % if no subset is a detected infeasible active set
                [obj,exitflag] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),candidateActiveSet);
                if ~exitflag
                
                    % if active set is not a childnode of detected non-primary active set
                    nonprimary = false;
                    if any(candidateActiveSet)  
                        [obj,nonprimary] = isNonprimary(obj,nonprimaryActiveSets(1:counter_nonprimaryActiveSets,:),candidateActiveSet);
                    end
                    if ~nonprimary

                        % solve optimality LP
                        [obj,t,exitflagOpt] = isOptimal(obj,candidateActiveSet);

                        % store in solution if optimal
                        if exitflagOpt==1
                            counter_solution = counter_solution+1;
                            if ~mod(counter_solution,chunkSize)
                                solution   = [solution; false(chunkSize,obj.q)]; % allocate more lines
                                degenerate = [degenerate; false(chunkSize,1)];   % allocate more lines
                            end
                            solution(counter_solution,:) = candidateActiveSet; 
                            
                            if t<=eps
                                % set degenerate flag to true if optimal and t=0
                                degenerate(counter_solution,:) = true;
                            end
                                
                            % for every element in the orbit of the processed active set
                            [obj,orbit] = determineOrbit(obj,candidateActiveSet);
                            for numberSymmetries = 1:size(orbit,1)                        

                                % if it is not a childnode of a detected non-primary active set
                                if ~all(orbit(numberSymmetries,:)==candidateActiveSet)
                                    nonprimary = false;
                                    if any(orbit(numberSymmetries,:))  
                                        [obj,nonprimary] = isNonprimary(obj,nonprimaryActiveSets(1:counter_nonprimaryActiveSets,:),orbit(numberSymmetries,:));
                                    end
                                    if ~nonprimary
                                        % add to set of non-primary active sets
                                        counter_nonprimaryActiveSets = counter_nonprimaryActiveSets+1;
                                        if ~mod(counter_nonprimaryActiveSets,chunkSize)
                                            nonprimaryActiveSets = [nonprimaryActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                        end
                                        nonprimaryActiveSets(counter_nonprimaryActiveSets,:) = orbit(numberSymmetries,:);
                                    end
                                end
                            end

                        % solve feasibility LP if not optimal
                        elseif exitflagOpt==-2
                            [obj,exitflagFeas] = isFeasible(obj,candidateActiveSet);

                            if exitflagFeas~=-2
                                
                                % for every element in the orbit of the processed active set
                                [obj,orbit] = determineOrbit(obj,candidateActiveSet);
                                for numberSymmetries = 1:size(orbit,1)                        

                                    % if it is not a childnode of a detected non-primary active set
                                    if ~all(orbit(numberSymmetries,:)==candidateActiveSet)
                                        nonprimary = false;
                                        if any(orbit(numberSymmetries,:))  
                                            [obj,nonprimary] = isNonprimary(obj,nonprimaryActiveSets(1:counter_nonprimaryActiveSets,:),orbit(numberSymmetries,:));
                                        end
                                        if ~nonprimary
                                            % add to set of non-primary active sets
                                            counter_nonprimaryActiveSets = counter_nonprimaryActiveSets+1;
                                            if ~mod(counter_nonprimaryActiveSets,chunkSize)
                                                nonprimaryActiveSets = [nonprimaryActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                            end
                                            nonprimaryActiveSets(counter_nonprimaryActiveSets,:) = orbit(numberSymmetries,:);
                                        end
                                    end
                                end
                                
                            % add to set of infeasible active sets if infeasible
                            else
                                % add active set itself to pruned
                                counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                                if ~mod(counter_infeasibleActiveSets,chunkSize)
                                    infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                end
                                infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidateActiveSet;


                                % add whole orbit to pruned if not already a superset
                                % for every element in the orbit of the processed active set
                                [obj,orbit] = determineOrbit(obj,candidateActiveSet);
                                for numberSymmetries = 1:size(orbit,1)
                                    % if it is not a childnode of a detected non-primary active set
                                    if ~all(orbit(numberSymmetries,:)==candidateActiveSet)
                                        exitflag = false;
                                        if any(orbit(numberSymmetries,:))  
                                            [obj,exitflag] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),orbit(numberSymmetries,:));
                                        end
                                        if ~exitflag
                                            % add to set of infeasible active sets
                                            counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                                            if ~mod(counter_infeasibleActiveSets,chunkSize)
                                                infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                            end
                                            infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidateActiveSet;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                % compute next active set in combinatorial tree
                [obj,candidateActiveSet,finished] = nextActiveSet(obj,candidateActiveSet,maxCardinality);
                if finished
                    break
                end
            end
            solution   = solution(1:counter_solution,:);   % delete unnecessary lines
            degenerate = degenerate(1:counter_solution,:); % delete unnecessary lines
        end
        
        function obj = determineTransformationMatrix(obj)
            %DETERMINETRANSFORMATIONMATRIX determines the constraint 
            %   transformation matrices (property obj.allTransformations)
            %   for each symmetric pair (Theta,Omega). Each sheet of 
            %   obj.allTransformations corresponds to the symmetry given by
            %   just that sheets in obj.allTheta and obj.allOmega (see 
            %   Mitze2022* for more details).
            %
            %
            %
            % *R. Mitze, M. Kvasnica, and M. Mönnigmann, "Exploiting
            %  symmetries in active set enumeration for constrained
            %  linear-quadratic optimal control", Automatica, vol. 151,
            %  p. 110900, 2023.
            
            % verstaendlicher constraint transformation matrix erklären.
            
            % initialize
            obj.allTransformations = false(obj.q,obj.q,obj.nSymm);
            
            % for each symmetric pair (Theta,Omega) and for each constraint i
            for k=1:obj.nSymm
                Theta = obj.allTheta(:,:,k);
                Omega = obj.allOmega(:,:,k);
                for i=1:obj.q
                    
                    % find symmetric constraint j
                    for j=1:obj.q
                        if norm(obj.G(i,:)/obj.w(i)-obj.G(j,:)*kron(eye(obj.N),Omega)/obj.w(j))<1e-6 && norm(obj.E(i,:)/obj.w(i)-obj.E(j,:)*Theta/obj.w(j))<1e-6 % !! tolerance manually chosen !!
                            
                            % store relation between symmetry k and the constraints i and j in transformation matrix
                            obj.allTransformations(i,j,k)=true;
                            break
                        end
                    end
                    
                    % case no symmetric constraint was found
                    if j==obj.q && ~obj.allTransformations(i,j,k)
                        error('no symmetric constraint was found. Run obj.testSymmetries to analyze your entered symmetric pairs in obj.allTheta and obj.allOmega')
                    end
                end
            end
        end
        
        function [obj,exitflag] = isNonprimary(obj,nonprimaryActiveSets,activeSet)
            %ISNONPRIMARY tests if an activeSet is a descendant of a 
            %   known nonprimary active set in the combinatorial tree.
            %
            %
            % Input Arguments:
            %   activeSet: active set in form of a logical line vector of 
            %   length obj.q (inactive constraint corresponds false, active
            %   constraint corresponds true)
            %
            %   nonprimaryActiveSets: each line contains a known nonprimary
            %   active set (inactive constraint corresponds false, active
            %   constraint corresponds true)
            %
            % Output Arguments:
            %   exitflag: is true if activeSet is a descendant of a known
            %   nonprimary active set in the combinatorial tree, and false
            %   otherwise
            
            % initialize
            exitflag = false;
            
            % test if active set is childnode of a set in nonprimaryActiveSets
            [~,lastIndices] = max(fliplr(nonprimaryActiveSets),[],2);
            lastIndices = obj.q+1-lastIndices;
            for i = unique(lastIndices)'
                nonprimaryActiveSets_i = nonprimaryActiveSets(lastIndices==i,:);
                if any(all(nonprimaryActiveSets_i(:,1:i)==repmat(activeSet(1:i),size(nonprimaryActiveSets_i,1),1),2))
                    exitflag = true;
                    break;
                end
            end
        end
        
        function [obj,orbit] = determineOrbit(obj,activeSet)
            %DETERMINEORBIT determines all active sets in the orbit of the
            %   active set that is entered.
            %
            %
            % Input Arguments:
            %   activeSet: active set must be expressed as a logical line
            %   vector of length obj.q (inactive constraint corresponds
            %   false, active constraint corresponds true)
            %
            % Output Arguments:
            %   orbit: lines of contain the active sets in the orbit
            %   (inactive constraint corresponds false, active constraint
            %   corresponds true)
            
            % initialize
            orbit = false(obj.nSymm,obj.q);
            
            % determine symmetric active set under all symmetries and add to orbit
            for i = 1:obj.nSymm
                orbit(i,:)= (any(obj.allTransformations(:,:,i)&activeSet,2))';
            end
            
            % delete redundant active sets in the orbit
            orbit = unique(orbit,'rows');
        end
        
        function plotPolytopes(obj)
            %PLOTPOLYTOPES plots polytopes defined by the active sets in
            %   obj.solution in current figure using Polyhedron from the
            %   multi-parametric toolbox*. Polytopes defined by active sets
            %   of reduced solution (one polytope per symmetry) are shown
            %   in green.
            %
            %
            % *M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
            %  Parametric Toolbox 3.0," Proc. of the European Control
            %  Conference, 2013, pp. 502-510.
            %
            % 
            % See also Polyhedron
            
            fprintf('\nPolytopes definded by active sets of reduced solution are shown in green.\n')
            
            % for all active sets in solution, plot polytopes defined by active sets in white color
            for activeSet_i=1:size(obj.solution,1)
                [T,d] = computePolytope(obj,obj.solution(activeSet_i,:));
                Poly  = Polyhedron(T,d);
                plot(Poly,'Color',[1 1 1]) % white
            end
            
            % for all active sets in reduced solution, plot polytopes defined by active sets in green color
            for activeSet_i=1:size(obj.solutionReduced,1)
                [T,d] = computePolytope(obj,obj.solutionReduced(activeSet_i,:));
                Poly  = Polyhedron(T,d);
                plot(Poly,'Color','green') % green
            end
        end
    end
        
end