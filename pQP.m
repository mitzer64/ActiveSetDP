classdef pQP
%PQP represents parametric quadratic programs (pQPs).
%
% Syntax:
%   myQP = PQP(H,G,S,w)
%   myQP = PQP(H,F,Y,G,E,w)
%
%
% Description:
%
%   myQP = pQP(H,G,S,w) creates an object for a QP of the form
%
%       min(z)   1/2*z'*H*z
%       s.t.     G*z <= w + S*x0
%       
%       with z = U+inv(H)*F'*x0.
%
%   myQP = pQP(H,F,Y,G,E,w) creates an object for a QP of the form
%
%       min(U)   1/2*U'*H*U + x0'*F*U + 1/2*x0'*Y*x0
%       s.t.     G*U <= w + E*x0.
%
%
% Examples:
%
%   create object for problem given in Tondel2003*
%     H = eye(3);
%     G = [1 0 -1; -1 0 -1; 0 1 -1; 0 -1 -1; zeros(4,3)];
%     S = [1 0; -1 0; 0 -1; 0 1; eye(2); -eye(2)];
%     w = [-ones(4,1); 5*ones(4,1)];
%     myQP = pQP(H,G,S,w)
%
%   create object for system from Gutman1987** for horizon 2
%     H = [9.18 4.15; 4.15 2.63];
%     F = [4.18 1.62; 11.17 4.96];
%     Y = [4.06 5.71; 5.71 15.03];
%     G = [1 0; -1 0; zeros(4,2); 0 1; 0 -1; 1 0; -1 0; .5 0; -.5 0; -2.2 -1.58; 2.2 1.58; .65 .3; -.65 -.3];
%     E = [zeros(2,2); 0 -1; 0 1; -1 0; 1 0; zeros(2,2); 0 -1; 0 1; -1 -1; 1 1; .62 2.5; -.62 -2.5; -.36 -.83; .36 .83];
%     w = [1 1 5 5 25 25 1 1 5 5 25 25 1 1 1 1]';
%     myQP = pQP(H,F,Y,G,E,w)
%
%
% *P. Tondel, T. A. Johansen, and A. Bemporad, "Further results on
%  multiparametric quadratic programming," 42nd IEEE International
%  Conference on Decision and Control (IEEE Cat. No. 03CH37475), vol. 3,
%  pp. 3173-3178, 2003.
%
% **P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
%  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
%  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
%  pp. 251–254, 1987.
%
% 
% See also pCLQOCP and pCLQOCPsymm

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
    
    properties (SetAccess = protected)
        H        % cost function 1/2*U'*H*U + x0'*F*U + 1/2*x0'*Y*x0 or 1/2*z'*H*z with z=U+inv(H)*F'*x0
        F        % cost function 1/2*U'*H*U + x0'*F*U + 1/2*x0'*Y*x0
        Y        % cost function 1/2*U'*H*U + x0'*F*U + 1/2*x0'*Y*x0
        G        % constraints G*U <= w + E*x0 or G*z <= w + S*x0 with z=U+inv(H)*F'*x0
        E        % constraints G*U <= w + E*x0
        S        % constraints G*z <= w + S*x0 with z=U+inv(H)*F'*x0
        w        % constraints G*U <= w + E*x0 or G*z <= w + S*x0 with z=U+inv(H)*F'*x0
        
        q        % number of constraints
        solution % solution of pQP expressed as set of active sets. Each line of solution contains an active set expressed with bit tuples.
    end
        
    methods
        function obj = pQP(varargin)
            %pQP constructs an instance of class pQP
            %
            % 
            % See also pQP
            
            switch nargin
                % case no inputs
                case 0
                    return
                    
                % case inputs are (H,G,S,w)
                case 4
                    % test if input H is correct
                    if  ~isnumeric(varargin{1}) || size(varargin{1},1)~=size(varargin{1},2) || ~isreal(varargin{1}) || ~all(eig(varargin{1})>0)
                        error('Input H must be a numerical positive definite square matrix')
                    end
                    obj.H = varargin{1};
                    
                    % test if input G is correct
                    if ~isnumeric(varargin{2}) || ~isreal(varargin{2}) || size(obj.H,2)~=size(varargin{2},2)
                        error('Input G must be a numerical matrix with the same number of columns as H')
                    end
                    obj.G = varargin{2};
                    obj.q = size(obj.G,1);
                    
                    % test if input S is correct
                    if ~isnumeric(varargin{3}) || ~isreal(varargin{3}) || obj.q~=size(varargin{3},1)
                        error('Input S must be a numerical matrix with the same number of lines as G')
                    end
                    obj.S = varargin{3};
                    
                    % test if input w is correct
                    if ~isnumeric(varargin{4}) || ~isreal(varargin{4}) || obj.q~=size(varargin{4},1) || size(varargin{4},2)~=1
                        error('Input w must be a numerical column vector with the same number of lines as G')
                    end
                    obj.w = varargin{4};

                % case inputs are (H,F,Y,G,E,w)
                case 6
                    % test if input H is correct
                    if  ~isnumeric(varargin{1}) || size(varargin{1},1)~=size(varargin{1},2) || ~isreal(varargin{1}) || ~all(eig(varargin{1})>0)
                        error('Input H must be a numerical positive definite square matrix')
                    end
                    obj.H = varargin{1};
                    
                    % test if input F is correct
                    if ~isnumeric(varargin{2}) || ~isreal(varargin{2}) || size(obj.H,2)~=size(varargin{2},2) || size(varargin{3},1)~=size(varargin{2},1)
                        error('Input F must be a numerical matrix with the same number of lines as Y and columns as H')
                    end
                    obj.F = varargin{2};
                    
                    % test if input Y is correct
                    if ~isnumeric(varargin{3}) || size(varargin{3},1)~=size(varargin{3},2) || ~isreal(varargin{3})
                        error('Input Y must be a numerical square matrix')
                    end
                    obj.Y = varargin{3};
                    
                    % test if input G is correct
                    if ~isnumeric(varargin{4}) || ~isreal(varargin{4}) || size(obj.H,2)~=size(varargin{4},2)
                        error('Input G must be a numerical matrix with the same number of columns as H')
                    end
                    obj.G = varargin{4};
                    obj.q = size(obj.G,1);
                    
                    % test if input E is correct
                    if ~isnumeric(varargin{5}) || ~isreal(varargin{5}) || size(obj.G,1)~=size(varargin{5},1) || size(obj.Y,2)~=size(varargin{5},2)
                        error('Input E must be a numerical matrix with the same number of lines as G and columns as Y')
                    end
                    obj.E = varargin{5};
                    obj.S = obj.E+obj.G*inv(obj.H)*obj.F';
                    
                    % test if input w is correct
                    if ~isnumeric(varargin{6}) || ~isreal(varargin{6}) || obj.q~=size(varargin{6},1) || size(varargin{6},2)~=1
                        error('Input w must be a numerical column vector with the same number of lines as G')
                    end
                    obj.w = varargin{6};
                        

                % case invalid number of inputs
                otherwise
                    error('inputs must be (H,G,S,w) or (H,F,Y,G,E,w)')
                    
            end
        end
        
        
        function obj = solveMPT(obj)
            %SOLVEMPT determines the property solution with the mpt
            %   toolbox* using the method solve from class Opt.
            % 
            % Syntax:
            %   myQP = myQP.SOLVEMPT
            %
            %
            % Examples:
            %
            %   create object of problem given in Tondel2003**
            %     H = eye(3);
            %     G = [1 0 -1; -1 0 -1; 0 1 -1; 0 -1 -1; zeros(4,3)];
            %     S = [1 0; -1 0; 0 -1; 0 1; eye(2); -eye(2)];
            %     w = [-ones(4,1); 5*ones(4,1)];
            %     myQP = pQP(H,G,S,w);
            %
            %   compute solution
            %     myQP = myQP.solveMPT;
            %     myQP.solution
            %
            %
            % *M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
            %  Parametric Toolbox 3.0," Proc. of the European Control
            %  Conference, 2013, pp. 502-510.
            %
            % **P. Tondel, T. A. Johansen, and A. Bemporad, "Further
            %  results on multiparametric quadratic programming," 42nd IEEE
            %  International Conference on Decision and Control (IEEE Cat.
            %  No. 03CH37475), vol. 3, pp. 3173-3178, 2003.
            %
            % 
            % See also Opt
            
            % test if MPT is in path
            if ~exist('Opt.m','file') && ~exist('Opt.m','class')
                error('solveMPT uses functions from Multi-Parametric Toolbox. Toolbox is not in path.')
            end
            
            fprintf('\ndetermine solution with multi-parametric toolbox\n')
            
            % compute explicit QP solution with MPT
            system       = Opt('H',obj.H,'f',[],'A',obj.G,'b',obj.w,'pB',obj.S);
            explSolution = system.solve;

            % compute Chebychev center of all polytopes of solution
            xc = zeros(size(obj.S,2),explSolution.xopt.Num);
            for i=1:explSolution.xopt.Num
                temp    = chebyCenter(explSolution.xopt.Set(i,1));
                xc(:,i) = temp.x;
            end

            % find active sets corresponding to Chebychev centers
            obj.solution = false(size(xc,2),size(obj.G,1));
            for i=1:size(xc,2)
                [~,obj.solution(i,:)] = solvePointwise(obj,xc(:,i));
            end
        end
        
        function obj = solveCombinatorial(obj)
            %SOLVECOMBINATORIAL determines the property solution with
            %   combinatorial mpQP (algorithm from Gupta2011*).
            % 
            % Syntax:
            %   myQP = myQP.SOLVECOMBINATORIAL
            %
            %
            % Examples:
            %
            %   create object of problem given in Tondel2003**
            %     H = eye(3);
            %     G = [1 0 -1; -1 0 -1; 0 1 -1; 0 -1 -1; zeros(4,3)];
            %     S = [1 0; -1 0; 0 -1; 0 1; eye(2); -eye(2)];
            %     w = [-ones(4,1); 5*ones(4,1)];
            %     myQP = pQP(H,G,S,w);
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
            % **P. Tondel, T. A. Johansen, and A. Bemporad, "Further
            %  results on multiparametric quadratic programming," 42nd IEEE
            %  International Conference on Decision and Control (IEEE Cat.
            %  No. 03CH37475), vol. 3, pp. 3173-3178, 2003.
            
            fprintf('\ndetermine solution with combinatorial mpQP\n')

            % maximal cardinality for processed active sets
            maxCardinality = min(size(obj.G,2),obj.q);

            % initialization
            chunkSize = 1000; % number of preallocated lines (overestimated guess)
            obj.solution         = false(chunkSize,obj.q);
            infeasibleActiveSets = false(chunkSize,obj.q);
            counter_solution             = 0;
            counter_infeasibleActiveSets = 0;
            
            % for every active set in combinatorial tree
            candidateActiveSet = false(1,obj.q);
            while true
                
                % if obj.G(active set,:) has full row rank
                if sum(candidateActiveSet)==rank(obj.G(candidateActiveSet,:))

                    % if no subset is a detected infeasible active set
                    [obj,exitflag] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),candidateActiveSet);
                    if ~exitflag
                        
                        % solve optimality LP
                        [obj,t,exitflagOpt] = isOptimal(obj,candidateActiveSet);
                        if exitflagOpt==1
                            
                            % store if optimal and t>0
                            if t>eps
                                counter_solution = counter_solution+1;
                                if ~mod(counter_solution,chunkSize)
                                    obj.solution = [obj.solution; false(chunkSize,obj.q)]; % allocate more lines
                                end
                                obj.solution(counter_solution,:) = candidateActiveSet;
                            
                            % compute polytope if optimal and t=0
                            else
                                [T,d] = computePolytope(obj,candidateActiveSet);
                                if isFullDim(Polyhedron(T,d))
                                    % store if polytope is full-dimensional
                                    counter_solution = counter_solution+1;
                                    if ~mod(counter_solution,chunkSize)
                                        obj.solution = [obj.solution; false(chunkSize,obj.q)]; % allocate more lines
                                    end
                                    obj.solution(counter_solution,:) = candidateActiveSet;
                                end
                            end

                        % solve feasibility LP if not optimal
                        elseif exitflagOpt==-2
                            [obj,exitflagFeas] = isFeasible(obj,candidateActiveSet);

                            % add to set of infeasible active sets if infeasible
                            if exitflagFeas==-2
                                counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                                if ~mod(counter_infeasibleActiveSets,chunkSize)
                                    infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                end
                                infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidateActiveSet;
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
            obj.solution = obj.solution(1:counter_solution,:); % delete unnecessary lines
        end
                
        
        function [U,activeSet] = solvePointwise(obj,x0)
            %SOLVEPOINTWISE solves quadratic program pointwise at given
            % initial state. Note that the property F must be defined.
            % 
            % Syntax:
            %   [U,activeSet] = myQP.SOLVEPOINTWISE(x0)
            %
            %
            % Examples:
            %
            %   create object for system from Gutman1987* for horizon 2
            %     H = [9.18 4.15; 4.15 2.63];
            %     F = [4.18 1.62; 11.17 4.96];
            %     Y = [4.06 5.71; 5.71 15.03];
            %     G = [1 0; -1 0; zeros(4,2); 0 1; 0 -1; 1 0; -1 0; .5 0; -.5 0; -2.2 -1.58; 2.2 1.58; .65 .3; -.65 -.3];
            %     E = [zeros(2,2); 0 -1; 0 1; -1 0; 1 0; zeros(2,2); 0 -1; 0 1; -1 -1; 1 1; .62 2.5; -.62 -2.5; -.36 -.83; .36 .83];
            %     w = [1 1 5 5 25 25 1 1 5 5 25 25 1 1 1 1]';
            %     myQP = pQP(H,F,Y,G,E,w)
            %
            %   compute optimal active set
            %     x0 = [8; -2.5];
            %     [U,activeSet] = myQP.solvePointwise(x0)
            % 
            %
            % Input Arguments:
            %   x0: initial state as numerical column vector with length
            %   obj.n
            % 
            % Output Arguments:
            %   U:  optimal input sequence at x0 is U = [u(0); ...; u(N-1)]
            %
            %   activeSet:  optimal active set at x0 expressed as logical
            %   line vector of length obj.q (inactive constraint
            %   corresponds false, active constraint corresponds true)
            %
            %
            % *P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
            %  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
            %  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
            %  pp. 251–254, 1987.
            %
            % 
            % See also quadprog
            
            % test if input x0 is correct
            if ~isnumeric(x0) || ~isreal(x0) || size(x0,1)~=size(obj.S,2) || size(x0,2)~=1
                error('x0 must be a numerical column vector with length obj.n')
            end
            
            % solve QP at x0
            options             = optimoptions('quadprog','Display','none');          
            [zOpt,~,exitflag,~] = quadprog(obj.H,zeros(size(obj.H,1),1),obj.G,obj.w+obj.S*x0,[],[],[],[],[],options);
            
            % if solution to QP exists
            if exitflag==1
                % determine optimal input sequence
                U = zOpt-obj.H\obj.F'*x0;
                % determine optimal active set
                activeSet = (abs(obj.G*zOpt-(obj.w+obj.S*x0))<=1e-6)'; % !! tolerance manually chosen !!
            else
                error(['no solution to quadratic program at x0 = [' num2str(x0') ']^T'])
            end
        end
        
        function [T,d] = computePolytope(obj,activeSet)
            %COMPUTEPOLYTOPE computes the polytope corresponding to an
            %   active set. Equations see, e.g., Lemma 2 in Jost2015*.
            % 
            % Syntax:
            %   [T,d] = myQP.COMPUTEPOLYTOPE(activeSet)
            %
            %
            % Examples:
            %
            %   create object from Tondel2003**
            %     H = eye(3);
            %     G = [1 0 -1; -1 0 -1; 0 1 -1; 0 -1 -1; zeros(4,3)];
            %     S = [1 0; -1 0; 0 -1; 0 1; eye(2); -eye(2)];
            %     w = [-ones(4,1); 5*ones(4,1)];
            %     myQP = pQP(H,G,S,w);
            %
            %   compute polytope
            %     activeSet = false(1,8);
            %     activeSet([2,4]) = true;
            %     [T,d] = myQP.computePolytope(activeSet)
            % 
            %
            % Input Arguments:
            %   active set must be a logical line vector of length obj.q
            %   (inactive constraint corresponds false, active constraint
            %   corresponds true)
            % 
            % Output Arguments:
            %   T and d such that Polytope = {x(0) | T*x(0)<=d}
            %
            %
            % *M. Jost, M. Schulze Darup, and M. Mönnigmann, "Optimal and
            %  suboptimal event-triggering in linear model predicitive 
            %  control," Proc. of the European Control Conference, 2015,
            %  pp. 1147-1152.
            %
            % **P. Tondel, T. A. Johansen, and A. Bemporad, "Further
            %  results on multiparametric quadratic programming," 42nd IEEE
            %  International Conference on Decision and Control (IEEE Cat.
            %  No. 03CH37475), vol. 3, pp. 3173-3178, 2003.
            
            % test if input activeSet is correct
            if ~islogical(activeSet) || size(activeSet,1)>1 || size(activeSet,2)~=obj.q
                error('active set must be a binary 1xq vector')
            end
            
            % determine T and d
            T = [obj.G(~activeSet,:)/obj.H*obj.G(activeSet,:)'/(obj.G(activeSet,:)/obj.H*obj.G(activeSet,:)')*obj.S(activeSet,:)-obj.S(~activeSet,:);...
                (obj.G(activeSet,:)/obj.H*obj.G(activeSet,:)')\obj.S(activeSet,:)];
            d = -[obj.G(~activeSet,:)/obj.H*obj.G(activeSet,:)'/(obj.G(activeSet,:)/obj.H*obj.G(activeSet,:)')*obj.w(activeSet,:)-obj.w(~activeSet,:);...
                (obj.G(activeSet,:)/obj.H*obj.G(activeSet,:)')\obj.w(activeSet,:)];
        end
        
        function [K,b] = computeControlLaw(obj,activeSet)
            %COMPUTECONTROLLAW computes the control law corresponding the
            %   activeSet.  Equations see, e.g., Lemma 2 in Jost2015*. Note
            %   that the property F must be defined.
            % 
            % Syntax:
            %   [K,b] = myQP.COMPUTECONTROLLAW(activeSet)
            %
            %
            % Examples:
            %
            %   create object for system from Gutman1987** for horizon 2
            %     H = [9.18 4.15; 4.15 2.63];
            %     F = [4.18 1.62; 11.17 4.96];
            %     Y = [4.06 5.71; 5.71 15.03];
            %     G = [1 0; -1 0; zeros(4,2); 0 1; 0 -1; 1 0; -1 0; .5 0; -.5 0; -2.2 -1.58; 2.2 1.58; .65 .3; -.65 -.3];
            %     E = [zeros(2,2); 0 -1; 0 1; -1 0; 1 0; zeros(2,2); 0 -1; 0 1; -1 -1; 1 1; .62 2.5; -.62 -2.5; -.36 -.83; .36 .83];
            %     w = [1 1 5 5 25 25 1 1 5 5 25 25 1 1 1 1]';
            %     myQP = pQP(H,F,Y,G,E,w)
            %
            %   compute control law
            %     activeSet = false(1,16);
            %     activeSet([7,13]) = true;
            %     [K,b] = myQP.computeControlLaw(activeSet)
            %
            %
            % Input Arguments:
            %   active set must be a logical line vector of length obj.q
            %   (inactive constraint corresponds false, active constraint
            %   corresponds true)
            % 
            % Output Arguments:
            %   K and b such that optimal control law U = is K*x(0)+b
            %
            %
            % *M. Jost, M. Schulze Darup, and M. Mönnigmann, "Optimal and
            %  suboptimal event-triggering in linear model predicitive 
            %  control," Proc. of the European Control Conference, 2015,
            %  pp. 1147-1152.
            %
            % **P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
            %  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
            %  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
            %  pp. 251–254, 1987.
            
            % test if input activeSet is correct and property obj.F is defined
            if ~islogical(activeSet) || size(activeSet,1)>1 || size(activeSet,2)~=obj.q
                error('active set must be a binary line vector of length obj.q')
            elseif isempty(obj.F)
                error('to compute the control law, the property F must be defined')
            end
            
            % determine K and b
            K = obj.H\obj.G(activeSet,:)'/(obj.G(activeSet,:)/obj.H*obj.G(activeSet,:)')*obj.S(activeSet,:)-obj.H\obj.F';
            b = obj.H\obj.G(activeSet,:)'/(obj.G(activeSet,:)/obj.H*obj.G(activeSet,:)')*obj.w(activeSet);
        end
        
        function plotSolution(obj)
            %PLOTSOLUTION illustrates state space partition for problems
            %   with two states in figure using Polyhedron from the
            %   multi-parametric toolbox*.
            % 
            % Syntax:
            %   myQP.PLOTSOLUTION
            %
            %
            % Examples:
            %
            %   create object of problem given in Tondel2003** and compute
            %   its solution
            %     H = eye(3);
            %     G = [1 0 -1; -1 0 -1; 0 1 -1; 0 -1 -1; zeros(4,3)];
            %     S = [1 0; -1 0; 0 -1; 0 1; eye(2); -eye(2)];
            %     w = [-ones(4,1); 5*ones(4,1)];
            %     myQP = pQP(H,G,S,w);
            %     myQP = myQP.solveCombinatorial;
            %
            %   plot partition
            %     myQP.plotSolution
            %
            %
            % *M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
            %  Parametric Toolbox 3.0," Proc. of the European Control
            %  Conference, 2013, pp. 502-510.
            %
            % **P. Tondel, T. A. Johansen, and A. Bemporad, "Further
            %  results on multiparametric quadratic programming," 42nd IEEE
            %  International Conference on Decision and Control (IEEE Cat.
            %  No. 03CH37475), vol. 3, pp. 3173-3178, 2003.
            %
            %
            % See also Polyhedron
            
            % test if MPT is in path
            if ~exist('Polyhedron.m','file') && ~exist('Polyhedron.m','class')
                error('solveMPT uses functions from Multi-Parametric Toolbox. Toolbox is not in path.')
            end
            
            % test if state space is two dimensional
            if size(obj.S,2)~=2
                error('plot only possible if number of states =2')
            % test if property obj.solution has been determined
            elseif isempty(obj.solution)
                error('compute property solution with obj = obj.solveMPT, obj = obj.solveCombinatorial or equivalent first')
            end

            fprintf('\nplot state space partition\n')
            
            % prepare figure
            hold on
            grid off
            xlabel('x_1')
            ylabel('x_2')

            plotPolytopes(obj)
        end
    end
    
    methods (Access = protected)        
        function [obj,nextActiveSet,finished] = nextActiveSet(obj,currentActiveSet,maxCardinality)
            %NEXTACTIVESET computes next active set in combinatorial tree*
            %   in the order of increasing cardinality and increasing
            %   constraint indices (in the combinatorial tree, this
            %   strategy proceeds from top to bottom, and in each level
            %   from left to right).
            %
            %
            % Input Arguments:
            %   currentActiveSet: current active set in the tree (logical
            %   line vector of length obj.q, inactive constraint = false,
            %   active constraint = true)
            %
            %   maxCardinality: maximal cardinality of the active sets in
            %   the combinatorial tree
            %
            % Output Arguments:
            %   nextActiveSet: next active set in the tree (inactive 
            %   constraint = false, active constraint = true)
            %
            %   finished: logical scalar, returns true if last active set
            %   in combinatorial tree is reached, and false otherwise
            %
            %
            % *R. Rymon, "Search through systematic set enumeration," 
            %  Proc. of the Third International Conf. on Principles of
            %  Knowledge Representation and Reasoning, pp. 539–550, 1992.

            % initialize
            nextActiveSet = currentActiveSet;
            finished = false;
            
            % root
            if ~max(nextActiveSet) 
                nextActiveSet(1) = true;
            
            % go to next child node of parent node in the same level
            elseif ~nextActiveSet(end) 
                temp = find(nextActiveSet,1,'last');
                nextActiveSet(temp)   = false;
                nextActiveSet(temp+1) = true;
            
            % go to next branch in the same level
            elseif ~min(nextActiveSet(find(nextActiveSet,1,'first'):end)) 
                temp  = find(~nextActiveSet,1,'last');
                temp1 = find(nextActiveSet(1:temp),1,'last');
                nextActiveSet(temp1:end) = [false true(1,length(nextActiveSet)-temp+1) false(1,temp-temp1-1)];
                
            % increase level of the tree
            elseif sum(nextActiveSet)<maxCardinality 
                temp = sum(nextActiveSet);
                nextActiveSet = [true(1,temp+1) false(1,length(nextActiveSet)-temp-1)];
            
            % reach last element of the tree
            else 
                finished = true;
            end 
        end
        
        function [obj,exitflag] = isPruned(obj,infeasibleActiveSets,activeSet)
            %ISPRUNED tests if an activeSet is a superset of a known
            %   infeasible active set.
            %
            %
            % Input Arguments:
            %   activeSet: active set in form of a logical line vector of length
            %   obj.q (inactive constraint = false, active constraint = true)
            %
            %   infeasibleActiveSets: matrix containing all known infeasible
            %   active sets
            %
            % Output Arguments:
            %   exitflag: returns true if activeSet is superset of a known
            %   infeasible active set, and false otherwise
            
            % initialize
            exitflag = false;
            
            % test if active set is superset of a set in infeasibleActiveSets
            for i = 1:size(infeasibleActiveSets,1)
                if ~any(infeasibleActiveSets(i,:)-activeSet>0)
                    exitflag = true;
                    break;
                end
            end 
        end
        
        function [obj,t,exitflag] = isOptimal(obj,activeSet)
            %ISOPTIMAL tests if there exists an x0 for that the activeSet
            %   corresponds to the solution to the QP at x0 by solving
            %   linear program (16) described in Gupta2011*.
            %
            %
            % Input Arguments:
            %   activeSet: active set in form of a logical line vector of length
            %   obj.q (inactive constraint = false, active constraint = true)
            %
            % Output Arguments:
            %   t: maximum value of cost function of linear program
            %   
            %   exitflag: exitflag of linear program (feasible returns 1,
            %   infeasible returns -2)
            %
            %
            % *A. Gupta, S. Bhartiya, and P. Nataraj, "A novel approach to
            %  multiparametric quadratic programming," Automatica, vol. 47,
            %  pp. 2112-2117, 2011.
            %
            % 
            % See also linprog

            % set options for linprog
            options = optimoptions('linprog','Display','none');

            % cost function is -t
            J = [zeros(size(obj.H,1)+size(obj.S,2)+length(obj.w),1); -1];

            % set linear equality constraints: active constraints, inactive constraints plus slack variables and stationarity conditions are fulfilled with equality
            Aeq = [obj.G(activeSet,:) -obj.S(activeSet,:) zeros(sum(activeSet),length(obj.w)+1);...
                obj.G(~activeSet,:) -obj.S(~activeSet,:) eye(sum(~activeSet)) zeros(sum(~activeSet),sum(activeSet)+1);...
                obj.H zeros(size(obj.H,1),size(obj.S,2)) zeros(size(obj.H,1),sum(~activeSet)) obj.G(activeSet,:)' zeros(size(obj.H,1),1)];
            Beq = [obj.w(activeSet,:);...
                obj.w(~activeSet,:);...
                zeros(size(obj.H,1),1)];

            % linear inequality constraints: t minus slack variables and t minus Langrangian multipliers result less or equal zero
            Aineq = [zeros(sum(~activeSet),size(obj.H,1)+size(obj.S,2)) -eye(sum(~activeSet)) zeros(sum(~activeSet),sum(activeSet)) ones(sum(~activeSet),1);...
                zeros(sum(activeSet),size(obj.H,1)+size(obj.S,2)+sum(~activeSet)) -eye(sum(activeSet)) ones(sum(activeSet),1)];
            Bineq = [zeros(sum(~activeSet),1);...
                zeros(sum(activeSet),1)];

            % bounds t, slack variables and Langrangian multiplies have lower bound set to zero
            lb = [-inf*ones(size(obj.S,2)+size(obj.H,1),1); zeros(length(obj.w)+1,1)];

            % solve LP (optimization variable is [z0 ... zNm1,x0,slack variables,Langrangian multiplies,t])
            [sol,~,exitflag,~] = linprog(J,Aineq,Bineq,Aeq,Beq,lb,[],options);

            if exitflag==1
                t = sol(end);
            else
                t = [];
                if exitflag~=-2
                    fprintf(['\nWarning: exitflag of LP in isOptimal for active set [' num2str(find(activeSet)) '] is ' num2str(exitflag) '\n']) 
                end
            end
        end
        
        function [obj,exitflag] = isFeasible(obj,activeSet)
            %ISFEASIBLE tests if there exists an x0 for that the active
            %   and inactive constraints of the active set are satisfied by
            %   solving linear program (18) described in Gupta2011*.
            %
            %
            % Input Arguments:
            %   activeSet: active set in form of a logical line vector of
            %   length obj.q (inactive constraint corresponds false, 
            %   active constraint corresponds true)
            %
            % Output Arguments:
            %   exitflag: exitflag of linear program (feasible returns 1,
            %   infeasible returns -2)
            %
            %
            % *A. Gupta, S. Bhartiya, and P. Nataraj, "A novel approach to
            %  multiparametric quadratic programming," Automatica, vol. 47,
            %  pp. 2112-2117, 2011.
            %
            % 
            % See also linprog

            % set options for linprog
            options = optimoptions('linprog','Display','none');
    
            % cost function is -t
            J = [zeros(size(obj.H,1)+size(obj.S,2)+sum(~activeSet),1); -1];

            % linear equality constraints: active constraints and inactive constraints plus slack variables are fulfilled with equality
            Aeq = [obj.G(activeSet,:) -obj.S(activeSet,:) zeros(sum(activeSet),sum(~activeSet)+1);...
                obj.G(~activeSet,:) -obj.S(~activeSet,:) eye(sum(~activeSet)) zeros(sum(~activeSet),1)];
            Beq = [obj.w(activeSet,:);...
                obj.w(~activeSet,:)];

            % linear inequality constraints: t minus slack variables is less or equal zero
            Aineq = [zeros(sum(~activeSet),size(obj.H,1)+size(obj.S,2)) -eye(sum(~activeSet)) ones(sum(~activeSet),1)];
            Bineq = zeros(sum(~activeSet),1);

            % lower bound to t and slack variables is zero
            lb = [-inf*ones(size(obj.S,2)+size(obj.H,1),1); zeros(sum(~activeSet)+1,1)];

            % solve LP (optimization variable is [z0 ... zNm1,x0,slack variables,t])
            [~,~,exitflag,~] = linprog(J,Aineq,Bineq,Aeq,Beq,lb,[],options);

            if exitflag~=-2 && exitflag~=1
                fprintf(['\nWarning: exitflag of LP in isFeasible for active set [' num2str(find(activeSet)) '] is ' num2str(exitflag) '\n'])
            end
        end
        
        function plotPolytopes(obj)
            %PLOTPOLYTOPES plots polytopes defined by the active sets in
            %   obj.solution in current figure using Polyhedron from the
            %   multi-parametric toolbox*.
            %
            %
            % *M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
            %  Parametric Toolbox 3.0," Proc. of the European Control
            %  Conference, 2013, pp. 502-510.
            %
            % 
            % See also Polyhedron
            
            % plot polytopes defined by active sets in solution in white color
            for activeSet_i=1:size(obj.solution,1)  
                [T,d] = computePolytope(obj,obj.solution(activeSet_i,:));
                Poly  = Polyhedron(T,d);
                plot(Poly,'Color',[1 1 1])
            end
        end
    end
        
end

