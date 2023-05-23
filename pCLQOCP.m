classdef pCLQOCP < pQP
%PCLQOCP represents constrained linear-quadratic optimal control
%  problems (OCPs) that can be formulated as parametric quadratic programs
%  (pQPs).
%
% Syntax:
%   myQP = PCLQOCP(A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N)
%   myQP = PCLQOCP(myLTISystem,N)
%
%
% Description:
%
%   myQP = pCLQOCP(A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N) creates an object for
%       a linear-quatratic OCP of the form
%
%       min(U,X) x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
%       s.t.     x(k+1) = A*x(k) + B*u(k), k=0,...,N-1
%                Ou*u(k) <= cu, k=0,...,N-1
%                Ox*x(k) <= cx, k=0,...,N-1
%                Ot*x(N) <= ct.
%
%       This optimal control problem is converted into a pQP with the
%       constraints ordered stagewise.
%
%       If the inputs P and (Ot,ct) are entered as empty, the terminal 
%       penalty P is set to the LQR penalty and the matrices for the 
%       terminal set (Ot,ct) are set to the LQR invariant set, respectively. 
%
%       If no terminal penalty and no terminal set must be considered, set
%       P to the zero matrix of appropriate size and ct to inf, respectively.
%
%   myQP = pCLQOCP(myLTISystem,N) creates an object for an object of
%       class LTISystem from the multi-parametric toolbox* and a given
%       horizon N,
%
%       min(U,X) x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
%       s.t.     x(k+1) = A*x(k) + B*u(k), k=0,...,N-1
%                umin <= u(k) <= umax, k=0,...,N-1
%                xmin <= x(k) <= xmax, k=0,...,N-1
%                Ot*x(N) <= ct.
%
%       This optimal control problem is converted into a pQP with the
%       constraints ordered stagewise.
%
%       If the input myLTISystem does not define a terminal penalty and a
%       terminal set, i.e., has no properties myLTISystem.x,'terminalPenalty'
%       and myLTISystem.x,'terminalSet', the terminal penalty is set to
%       the LQR penalty and the terminal set is set to the LQR set, respectively.
%
%       If no terminal penalty and no terminal set must be considered, set
%       myLTISystem.x.with('terminalPenalty'); and
%       myLTISystem.x.with('terminalSet');, respectively, without further 
%       specifying them.
%
%
% Examples:
%
%   create object for a double integrator system from Gutman1987**
%      A    = [1 1; 0 1];
%      B    = [.5; 1];
%      Q    = [1 0; 0 1];
%      R    = .1;
%      N    = 5; 
%      Ox   = [0 1;0 -1;1 0; -1 0];
%      cx   = [5 5 25 25]';
%      Ou   = [1; -1];
%      cu   = [1 1]';
%      myQP = pCLQOCP(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N)
%
%   create object for a double integrator system from Gutman1987**
%     myLTISystem           = LTISystem('A',[1 1; 0 1],'B',[.5; 1]);
%     myLTISystem.x.min     = [-25; -5];
%     myLTISystem.x.max     = [25; 5];
%     myLTISystem.u.min     = -1;
%     myLTISystem.u.max     = 1;
%     myLTISystem.x.penalty = QuadFunction([1 0; 0 1]);
%     myLTISystem.u.penalty = QuadFunction(.1);
%     N = 5;
%     myQP = pCLQOCP(myLTISystem,N)
%
%
% *M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
%  Parametric Toolbox 3.0," Proc. of the European Control
%  Conference, 2013, pp. 502-510.
%
% **P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
%  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
%  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
%  pp. 251–254, 1987.
%
% 
% See also LTISystem and pQP and pCLQOCPsymm

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
    
    properties
        N   % prediction horizon
    end
    
    properties (SetAccess = private)
        A   % system matrix of linear discrete-time system x(k+1) = A*x(k) + B*u(k)     
        B   % input matrix of linear discrete-time system x(k+1) = A*x(k) + B*u(k)
        Q   % state penalty matrix in cost function x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
        R   % input penalty matrix in cost function x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
        P   % terminal state penalty matrix in cost function x(N)'*P*x(N) + sum(k=0,...,N-1) x(k)'*Q*x(k) + u(k)'*R*u(k)
        Ox  % state constraints are Ox*x(k) <= cx, k=0,...,N-1
        cx  % state constraints are Ox*x(k) <= cx, k=0,...,N-1
        Ou  % input constraints are Ou*u(k) <= cu, k=0,...,N-1
        cu  % input constraints are Ou*u(k) <= cu, k=0,...,N-1
        Ot  % terminal state constraints are Ot*x(N) <= ct
        ct  % terminal state constraints are Ot*x(N) <= ct
        
        n   % number of states
        m   % number of inputs
        qUX % number of state- and input constraints per stage
        qT  % number of terminal constraints
    end
    
    methods
        function obj = set.N(obj,new_N)
            % test if input N is correct
            if ~isscalar(new_N) || ~isnumeric(new_N) || ~isreal(new_N) || ceil(new_N)~=floor(new_N) || new_N<0
                error('horizon N must be a positive integer')
            else
                obj.N = new_N;
                % build matrices for pQP
                obj   = ocp2mpqp(obj);
            end
        end
        
        
        function obj = pCLQOCP(varargin)
            %PCLQOCP constructs an instance of class pCLQOCP
            %
            % 
            % See also pCLQOCP
            
            obj = obj@pQP();
            
            switch nargin                    
                % case inputs are (myLTISystem,N)
                case 2
            
                    % test if MPT is in path
                    if ~exist('LTISystem.m','file') && ~exist('LTISystem.m','class')
                        error('LTISystem is a class from Multi-Parametric Toolbox. Toolbox is not in path.')
                    end
                    
                    % test if input myLTISystem is correct
                    if  ~isa(varargin{1},'LTISystem')
                        error('Input myLTISystem must be an object of class LTISystem')
                    end
                    myLTISystem = varargin{1};
                    
                    % test if A,B is controllable
                    if rank(ctrb(myLTISystem.A,myLTISystem.B))~=myLTISystem.nx
                        error('System resulting from inputs (A,B) must be controllable')
                    end
                    
                    obj.A  = myLTISystem.A;
                    obj.B  = myLTISystem.B;
                    obj.Q  = myLTISystem.x.penalty.H;
                    obj.R  = myLTISystem.u.penalty.H;
                    obj.n  = myLTISystem.nx;
                    obj.m  = myLTISystem.nu;
                    
                    % test if inputs Ox and cx are correct
                    obj.Ox = [-eye(obj.n); eye(obj.n)];
                    if all([-myLTISystem.x.min; myLTISystem.x.max]>0)
                        obj.cx = [-myLTISystem.x.min; myLTISystem.x.max];
                    else
                        error('Origin of state space does not satisfy the state constraints')
                    end
                    
                    % test if inputs Ou and cu are correct
                    obj.Ou = [-eye(obj.m); eye(obj.m)];
                    if all([-myLTISystem.u.min; myLTISystem.u.max]>0)
                        obj.cu = [-myLTISystem.u.min; myLTISystem.u.max];
                    else
                        error('Origin of input space does not satisfy the input constraints')
                    end
                    obj.qUX = length(find(abs([obj.cx; obj.cu])<inf));

                    % case no terminal penalty is to be considered
                    if isprop(myLTISystem.x,'terminalPenalty') && isempty(myLTISystem.x.terminalPenalty)
                        obj.P = zeros(obj.n,obj.n);
                    else
                        % case terminal penalty is set to LQR penalty
                        if ~isprop(myLTISystem.x,'terminalPenalty')
                            fprintf('\nNote: Terminal penalty is set to LQR penalty (see doc pCLQOCP for more information)\n')
                            myLTISystem.x.with('terminalPenalty');
                            myLTISystem.x.terminalPenalty = myLTISystem.LQRPenalty;
                        end
                        obj.P = myLTISystem.x.terminalPenalty.H;
                    end

                    % case no terminal set is to be considered
                    if isprop(myLTISystem.x,'terminalSet') && isempty(myLTISystem.x.terminalSet)
                        obj.Ot = [];
                        obj.ct = [];
                    else
                        % case terminal set is set to LQR set
                        if ~isprop(myLTISystem.x,'terminalSet')
                            fprintf('\nNote: Terminal set is set to LQR invariant set (see doc pCLQOCP for more information)\n')
                            myLTISystem.x.with('terminalSet');
                            myLTISystem.x.terminalSet = myLTISystem.LQRSet;
                        end
                        obj.Ot = myLTISystem.x.terminalSet.H(:,1:end-1);
                        obj.ct = myLTISystem.x.terminalSet.H(:,end);
                    end
                    obj.qT = length(find(abs(obj.ct)<inf));
                    
                    % test if input N is correct
                    if ~isscalar(varargin{2}) || ~isnumeric(varargin{2}) || ~isreal(varargin{2}) || ceil(varargin{2})~=floor(varargin{2}) || varargin{2}<0
                        error('Input N must be a positive integer')
                    end
                    obj.N = varargin{2}; % setting N triggers computation of QP matrices

                % case inputs are (A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N)
                case 12
                    % test if input A is correct
                    if ~isnumeric(varargin{1}) || ~isreal(varargin{1}) || size(varargin{1},1)~=size(varargin{1},2)
                        error('Input A must be a numerical square matrix')
                    end
                    obj.A  = varargin{1};
                    obj.n = size(obj.A,1);
                    
                    % test if input B is correct
                    if ~isnumeric(varargin{2}) || ~isreal(varargin{2}) || size(varargin{2},1)~=obj.n
                        error('Input B must be a numerical matrix with the same number of lines as G')
                    end
                    obj.B  = varargin{2};
                    obj.m = size(obj.B,2);
                    
                    % test if A,B is controllable
                    if rank(ctrb(obj.A,obj.B))~=obj.n
                        error('System resulting from inputs (A,B) must be controllable')
                    end
                    
                    % test if input Q is correct
                    if ~isnumeric(varargin{3}) || ~isreal(varargin{3}) || size(varargin{3},1)~=size(varargin{3},2) || size(varargin{3},1)~=obj.n || any(sign(eig((varargin{3}+varargin{3}')/2))==-1)
                        error('Input Q must be a numerical positive semidefinit square matrix with the same size as A')
                    end
                    obj.Q  = varargin{3};
                    
                    % test if input R is correct
                    if ~isnumeric(varargin{4}) || ~isreal(varargin{4}) || size(varargin{4},1)~=size(varargin{4},2) || size(varargin{4},1)~=obj.m || ~all(sign(eig((varargin{3}+varargin{3}')/2))==1)
                        error('Input R must be a numerical positive definit square matrix with the same number of columns as B')
                    end
                    obj.R  = varargin{4};
                    
                    % test if input P is correct
                    if ~isempty(varargin{5}) && (~isnumeric(varargin{5}) || ~isreal(varargin{5}) || size(varargin{5},1)~=size(varargin{5},2) || size(varargin{5},1)~=obj.n)
                        error('Input P must be empty or a numerical square matrix with the same size as Q')
                    end 
                    if isempty(varargin{5})
                        % case terminal penalty is set to LQR penalty
                        fprintf('\nNote: Terminal penalty P is set to LQR penalty (see doc pCLQOCP for more information)\n')
                        [~,obj.P,~] = dlqr(obj.A,obj.B,obj.Q,obj.R,zeros(obj.n,obj.m));
                    else
                        obj.P  = varargin{5};
                    end
                    
                    % test if input Ox is correct
                    if ~isnumeric(varargin{6}) || ~isreal(varargin{6}) || size(varargin{6},2)~=obj.n
                        error('Input Ox must be a numerical matrix with the same number of columns as A')
                    end 
                    obj.Ox = varargin{6};
                    
                    % test if input cx is correct
                    if ~isnumeric(varargin{7}) || ~isreal(varargin{7}) || size(varargin{7},2)~=1 || size(varargin{7},1)~=size(obj.Ox,1)
                        error('Input cx must be a numerical column vector with the same number of lines as Ox')
                    end
                    if any(varargin{7}<=0)
                        error('Inputs Ox and cx must define a full-dimensional polytope Ox*x<=cx that contains the origin in its interior')
                    end 
                    obj.cx = varargin{7};
                    
                    % test if input Ou is correct
                    if ~isnumeric(varargin{8}) || ~isreal(varargin{8}) || size(varargin{8},2)~=obj.m
                        error('Input Ou must be a numerical matrix with the same number of columns as B')
                    end 
                    obj.Ou = varargin{8};
                    
                    % test if input cu is correct
                    if ~isnumeric(varargin{9}) || ~isreal(varargin{9}) || size(varargin{9},2)~=1 || size(varargin{9},1)~=size(obj.Ou,1)
                        error('Input cu must be a numerical column vector with the same number of lines as Ou')
                    end
                    if any(varargin{9}<=0)
                        error('Inputs Ou and cu must define a full-dimensional polytope Ou*u<=cu that contains the origin in its interior')
                    end 
                    obj.cu = varargin{9};
                    obj.qUX = length(find(abs([obj.cx; obj.cu])<inf));
                    
                    % case no terminal set is to be considered (input is ct=inf)
                    if ~isempty(varargin{11}) && all(varargin{11}==inf)
                        obj.Ot = [];
                        obj.ct = [];
                    % case terminal set is set to LQR set (inputs are Ot=ct=[])
                    elseif isempty(varargin{10}) && isempty(varargin{11})
                        fprintf('\nNote: Terminal set Ot*x<=ct is set to LQR invariant set (see doc pCLQOCP for more information)\n')
                        [obj.Ot,obj.ct] = determineTerminalSet(obj);
                    % case user specified terminal set
                    else
                        % test if input Ot is correct
                        if ~isnumeric(varargin{10}) || ~isreal(varargin{10}) || size(varargin{10},2)~=obj.n
                            error('Input Ot must be a numerical matrix with the same number of columns as A. Enter Ot=ct=[] to set the terminal set as the LQR invariant set.')
                        end
                        obj.Ot = varargin{10};
                        
                        % test if input ct is correct
                        if ~isnumeric(varargin{11}) || ~isreal(varargin{11}) || size(varargin{11},2)~=1 || size(varargin{11},1)~=size(obj.Ot,1)
                            error('Input ct must be a numerical column vector with the same number of lines as Ot. Enter Ot=ct=[] to set the terminal set as the LQR invariant set.')
                        end
                        if any(varargin{11}<=0)
                            error('Inputs Ot and ct must define a full-dimensional polytope Ot*x<=ct that contains the origin in its interior')
                        end 
                        obj.ct = varargin{11};
                    end
                    obj.qT  = length(obj.ct);
                    
                    % test if input N is correct
                    if ~isscalar(varargin{12}) || ~isnumeric(varargin{12}) || ~isreal(varargin{12}) || ceil(varargin{12})~=floor(varargin{12}) || varargin{12}<0
                        error('Input N must be a positive integer')
                    end
                    obj.N  = varargin{12};

                % case invalid inputs
                otherwise
                    error('Inputs musst follow either (myLTISystem,N) or (A,B,Q,R,P,Ox,cx,Ou,cu,Ot,ct,N)'); 

            end
        end
        
        
        function obj = solveDP(obj)
            %SOLVEDP determines the property solution using the dynamic  
            %   programming algorithm from Mitze2020* und the cardinality 
            %   limit for processed active sets introduced in Mitze2022**.
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
            %   Double integrator from Gutman1987***
            %      A    = [1 1; 0 1];
            %      B    = [.5; 1];
            %      Q    = [1 0; 0 1];
            %      R    = .1;
            %      N    = 20; 
            %      Ox   = [0 1;0 -1;1 0; -1 0];
            %      cx   = [5 5 25 25]';
            %      Ou   = [1; -1];
            %      cu   = [1 1]';
            %      myQP = pCLQOCP(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N)
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
            % **R. Mitze and M. Mönnigmann, "Improved active set dynamic
            %  programming for solving linear-quadratic optimal control
            %  problems", 2022 IEEE 61st Conference on Decision and Control
            %  (CDC), pp. 1764-1769, 2022.
            %
            % ***P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
            %  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
            %  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
            %  pp. 251–254, 1987.
            
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

            % determine solution for initial horizon N=1
            obj.N = 1;
            fprintf('\ndetermine solution for initial horizon N=1\n')
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
                        
            obj.solution = solution;
        end
        
        
        function [X,U,activeSet] = computeOpenLoopTrajectory(obj,myInput)
            % COMPUTEOPENLOOPTRAJECTORY computes open loop trajectory from
            %  given initial state or Chebychev center of a given active 
            %  set.
            % 
            % Syntax:
            %   [X,U,activeSet] = myQP.COMPUTEOPENLOOPTRAJECTORY(x0)
            %   [X,U,activeSet] = myQP.COMPUTEOPENLOOPTRAJECTORY(activeSet)
            %
            %
            % Examples:
            %
            %   Double integrator from Gutman1987*
            %      A    = [1 1; 0 1];
            %      B    = [.5; 1];
            %      Q    = [1 0; 0 1];
            %      R    = .1;
            %      N    = 5; 
            %      Ox   = [0 1;0 -1;1 0; -1 0];
            %      cx   = [5 5 25 25]';
            %      Ou   = [1; -1];
            %      cu   = [1 1]';
            %      myQP = pCLQOCP(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N);
            %
            %   compute trajectory (Syntax 1)
            %     x0 = [-2; 3];
            %     X = myQP.computeOpenLoopTrajectory(x0)
            %
            %   compute trajectory (Syntax 2)
            %     activeSet = false(1,34);
            %     activeSet([1,7,13,19]) = true;
            %     X = myQP.computeOpenLoopTrajectory(activeSet)
            % 
            %
            % Input Arguments:
            %   activeSet: logical line vector with length obj.n (inactive
            %   constraint = false, active constraint = true)
            %
            %   x0: numerical column vector with length obj.n
            % 
            % Output Arguments:
            %   X: open loop state trajectory is X = [x(0) x(1) ... x(N)]
            %
            %   U: open loop optimal input sequence is U = [u(0) u(1) ... u(N)]
            %
            %   activeSet: optimal active set at state x(0) as logical line
            %   vector with length obj.n (inactive constraint = false, 
            %   active constraint = true)
            %
            %
            % *P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
            %  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
            %  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
            %  pp. 251–254, 1987.
            
            % if input is x0
            if isnumeric(myInput) && size(myInput,2)==1 && length(myInput)==obj.n
                x0 = myInput;
                
                % compute optimal active set and optimal input sequence at x0
                [U,activeSet] = solvePointwise(obj,x0);
                
            % if input is activeSet
            elseif min(islogical(myInput)) && size(myInput,1)==1 && length(myInput)==obj.q
                activeSet = myInput;
                
                [~,~,exitflagOpt] = isOptimal(obj,activeSet);
                if exitflagOpt==-2
                    error('input active set is not element of the solution')
                end
                
                % compute Chebychev center of polytope defined by active set 
                [T,d]  = computePolytope(obj,activeSet);
                Poly   = Polyhedron(T,d);
                Center = chebyCenter(Poly);
                x0     = Center.x;
            
                % compute control law and optimal input sequence defined by active set
                [K,b] = computeControlLaw(obj,activeSet);
                U     = K*x0+b;
                
            % otherwise
            else
                error('input must be a state x0 (numerical column vector of length obj.n) or an active set (logical line vector of length obj.q)')
            end
                
            % compute trajectory
            X = [x0 zeros(obj.n,obj.N)];
            for k_i=1:obj.N
                X(:,k_i+1)=obj.A*X(:,k_i)+obj.B*U((k_i-1)*obj.m+1:k_i*obj.m);
            end
        end
        
        function [xTraj,uTraj,activeSets,costTraj] = computeClosedLoopTrajectory(obj,myInput,numberSteps)
            % COMPUTECLOSEDLOOPTRAJECTORY computes closed loop trajectory
            %  from initial state or Chebychev center of an active set for
            %  a discrete number of time steps.
            % 
            % Syntax:
            %   [xTraj,uTraj,activeSets,costTraj] = myQP.COMPUTECLOSEDLOOPTRAJECTORY(x0,numberSteps)
            %   [xTraj,uTraj,activeSets,costTraj] = myQP.COMPUTECLOSEDLOOPTRAJECTORY(activeSet,numberSteps)
            %
            %
            % Examples:
            %
            %   Double integrator from Gutman1987*
            %      A    = [1 1; 0 1];
            %      B    = [.5; 1];
            %      Q    = [1 0; 0 1];
            %      R    = .1;
            %      N    = 5; 
            %      Ox   = [0 1;0 -1;1 0; -1 0];
            %      cx   = [5 5 25 25]';
            %      Ou   = [1; -1];
            %      cu   = [1 1]';
            %      myQP = pCLQOCP(A,B,Q,R,[],Ox,cx,Ou,cu,[],[],N);
            %
            %   compute trajectory (Syntax 1)
            %     x0 = [-2; 3];
            %     [xTraj,uTraj,activeSets,costTraj] = myQP.computeClosedLoopTrajectory(x0,10)
            %
            %   compute trajectory (Syntax 2)
            %     activeSet = false(1,34);
            %     activeSet([1,7,13,19]) = true;
            %     [xTraj,uTraj,activeSets,costTraj] = myQP.computeClosedLoopTrajectory(activeSet,10)
            % 
            %
            % Input Arguments:
            %   activeSet: active set as logical line vector of length
            %   obj.n (inactive constraint = false, active constraint = 
            %   true)
            %
            %   x0: initial state as numerical column vector with length 
            %   obj.n
            %
            %   numberSteps: number of time steps
            % 
            % Output Arguments:
            %   xTraj: states x(0) on closed loop trajectory are xTraj = 
            %   [x(0,timeStep=0) ... x(0,timeStep=numberSteps)]
            %
            %   uTraj: first optimal input for states on xTraj are uTraj =
            %   [u(0,timeStep=0) ... u(0,timeStep=numberSteps-1)]
            %
            %   activeSet: ith line contains optimal active set at ith 
            %   state of xTraj (inactive constraint = false, active 
            %   constraint = true)
            %
            %   costTraj: ith entry contains minimal value of cost function
            %   at ith state of xTraj
            %
            %
            % *P. O. Gutman and M. Cwikel, "An Algorithm to Find Maximal State
            %  Constraint Sets for Discrete-Time Linear Dynamical Systems with Bounded
            %  Controls and States", IEEE Transactions on Automatic Control, vol. 32,
            %  pp. 251–254, 1987.
            
            % if input is x0
            if isnumeric(myInput) && size(myInput,2)==1 && length(myInput)==obj.n
                x0 = myInput;
                
                % compute optimal active set and optimal input sequence at x0
                [U,activeSet] = solvePointwise(obj,x0);
                
            % if input is activeSet
            elseif min(islogical(myInput)) && size(myInput,1)==1 && length(myInput)==obj.q
                activeSet = myInput;
                
                [~,~,exitflagOpt] = isOptimal(obj,activeSet);
                if exitflagOpt==-2
                    error('input active set is not element of the solution')
                end
                
                % compute Chebychev center of polytope defined by active set 
                [T,d]  = computePolytope(obj,activeSet);
                Poly   = Polyhedron(T,d);
                Center = chebyCenter(Poly);
                x0     = Center.x;
                
                % compute control law and optimal input defined by active set
                [K,b] = computeControlLaw(obj,activeSet);
                U = K*x0+b;
                
            % otherwise
            else
                error('Inputs musst follow either (x0,numberSteps) or (activeSet,numberSteps)')
            end
            
            % initialize
            xTraj = [x0 zeros(obj.n,numberSteps)];
            uTraj = [U(1:obj.m) zeros(obj.m,numberSteps-1)];
            activeSets = [activeSet; false(numberSteps-1,obj.q)];
            costTraj = [1/2*U'*obj.H*U+x0'*obj.F*U+1/2*x0'*obj.Y*x0 zeros(1,numberSteps-1)];
            
            % for increasing time steps
            for i=1:numberSteps
            
                % compute subsequent state by applying the first optimal input
                xTraj(:,i+1) = obj.A*xTraj(:,i)+obj.B*U(1:obj.m);
                
                if i<numberSteps
                    % compute optimal input and optimal active set at new state
                    [U,activeSets(i+1,:)] = solvePointwise(obj,xTraj(:,i+1));
                    uTraj(:,i+1) = U(1:obj.m);

                    % compute value of cost function at new state
                    costTraj(:,i+1) = 1/2*U'*obj.H*U+xTraj(:,i+1)'*obj.F*U+1/2*xTraj(:,i+1)'*obj.Y*xTraj(:,i+1);
                end
            end
        end
    end
    
    methods (Access = protected)
        function obj = ocp2mpqp(obj)
            %OCP2MPQP converts a constrained linear-quadratic optimal
            %  control problem into a quadratic program and sets properties
            %  H, F, Y, G, w, S, and E. Equations see, e.g., Appendix C in
            %  dissertation Jost2015*. The constraints of the quadratic
            %  program are ordered by increasing stage.
            %
            %
            % *M. Jost, "Accelerating the calculation of model predictive
            %  control laws for constrained linear systems: From the
            %  explicit solution to fast online calculations,"
            %  dissertation, Ruhr-Universität Bochum, 2015.

            qU = length(obj.cu);
            qX = length(obj.cx); 

            % insert linear system
            calA = zeros(obj.n*(obj.N+1),obj.n);
            calB = zeros(obj.n*(obj.N+1),obj.m*obj.N);
            calQ = zeros(obj.n*(obj.N+1),obj.n*(obj.N+1));
            calR = zeros(obj.m*obj.N,obj.m*obj.N);
            temp = zeros(2*obj.m*obj.N,obj.m*obj.N); 
            for i = 1:obj.N
                calA(1+obj.n*(i-1):obj.n*i,:) = obj.A^(i-1);
                for j = 1:obj.N
                    if j<=i
                        calB(1+obj.n*i:(i+1)*obj.n,1+obj.m*(j-1):obj.m*j) = obj.A^(i-j)*obj.B;
                    end
                end
                calQ(1+obj.n*(i-1):obj.n*i,1+obj.n*(i-1):obj.n*i) = obj.Q; 
                calR(1+obj.m*(i-1):obj.m*i,1+obj.m*(i-1):obj.m*i) = obj.R; 
                temp(1+obj.m*(2*i-2):obj.m*2*i,1+obj.m*(i-1):obj.m*i) = [-eye(obj.m); eye(obj.m)];
            end
            calA(1+obj.n*obj.N:obj.n*(obj.N+1),:) = obj.A^obj.N; 
            calQ(1+obj.n*obj.N:obj.n*(obj.N+1),1+obj.n*obj.N:obj.n*(obj.N+1)) = obj.P;

            % weighting matrices
            obj.Y = calA'*calQ*calA;
            obj.F = calA'*calQ*calB;
            obj.H = calB'*calQ*calB+calR;

            % input constraints
            GU = obj.Ou;
            for i=1:obj.N-1
                GU = blkdiag(GU,obj.Ou);
            end
            EU = zeros(qU*obj.N,obj.n);           
            WU = repmat(obj.cu,obj.N,1);

            % state constraints             
            calJ = obj.Ox;
            for i=1:obj.N-1
                calJ = blkdiag(calJ,obj.Ox);
            end
            calJ = [calJ zeros(size(calJ,1),obj.n)];
            GX = calJ*calB;
            EX = -calJ*calA;
            WX = repmat(obj.cx,obj.N,1);

            % order state and input constraints by stage 
            obj.G = zeros(obj.N*(qU+qX)+obj.qT,obj.N*obj.m);
            obj.E = zeros(obj.N*(qU+qX)+obj.qT,obj.n);
            obj.w = zeros(obj.N*(qU+qX)+obj.qT,1);
            for i=1:obj.N
                % stage i: u constraints
                obj.G((i-1)*(qU+qX)+(1:qU),:) = GU((i-1)*qU+(1:qU),:); 
                obj.E((i-1)*(qU+qX)+(1:qU),:) = EU((i-1)*qU+(1:qU),:); 
                obj.w((i-1)*(qU+qX)+(1:qU),:) = WU((i-1)*qU+(1:qU),:); 

                % stage i: x constraints
                obj.G((i-1)*(qU+qX)+(qU+1:(qU+qX)),:) = GX((i-1)*qX+(1:qX),:); 
                obj.E((i-1)*(qU+qX)+(qU+1:(qU+qX)),:) = EX((i-1)*qX+(1:qX),:); 
                obj.w((i-1)*(qU+qX)+(qU+1:(qU+qX)),:) = WX((i-1)*qX+(1:qX),:); 
            end 

            % add terminal constraints
            if obj.qT>0
                obj.G(end-obj.qT+1:end,:) = obj.Ot*[zeros(obj.n,obj.N*obj.n) eye(obj.n)]*calB;
                obj.E(end-obj.qT+1:end,:) = -obj.Ot*[zeros(obj.n,obj.N*obj.n) eye(obj.n)]*calA;
                obj.w(end-obj.qT+1:end,:) = obj.ct;
            end
            
            % delete unconstrained states and inputs
            unconstrained = find(abs(obj.w)==inf);
            obj.G(unconstrained,:) = [];
            obj.E(unconstrained,:) = [];
            obj.w(unconstrained,:) = [];
                
            % introduce augmented parameter z:
            obj.S = obj.E+obj.G/obj.H*obj.F';
            
            % total number of constraints
            obj.q = size(obj.G,1);
            
            obj.solution = [];
        end
        
        function [obj,solution,degenerate] = solveCombinatorialConsiderDegeneracy(obj,maxCardinality)
            %SOLVECOMBINATORIALCONSIDERDEGENERACY determines all optimal
            %   active sets with combinatorial mpQP (Alg. 3 from Mitze2020*)
            %   and stores them in solution. Active sets such that
            %   obj.G(active set,:) is not of full row rank are marked.
            % 
            %
            % Input Arguments:
            %   maxCardinality: maximal cardinality of active sets
            %
            % Output Arguments:
            %   solution: matrix containing all optimal active sets 
            %   (expressed as logical line vectors of length obj.q,
            %   inactive constraint = false, active constraint = true)
            %
            %   degenerate: logical vector. Every entry corresponds to a
            %   line in solution. If the active set in that line has the
            %   solution t=0 to the optimality LP, then the entry is true
            %
            %
            % *R. Mitze and M. Mönnigmann, "A dynamic programming
            %  approach to solving constrained linear-quadratic optimal
            %  control problems", Automatica, vol. 120, p. 109132, 2020.

            % save maximal cardinality for processed active sets
            maxCardinality = min(maxCardinality,obj.q);
            
            % initialize
            chunkSize = 1000; % number of preallocated lines (overestimated guess)
            solution             = false(chunkSize,obj.q);
            infeasibleActiveSets = false(chunkSize,obj.q);
            degenerate           = false(chunkSize,1);
            counter_solution             = 0;
            counter_infeasibleActiveSets = 0;
            
            % for every active set in the combinatorial tree
            candidateActiveSet = false(1,obj.q);
            while true

                % if no subset is a detected infeasible active set
                [obj,exitflag] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),candidateActiveSet);
                if ~exitflag

                    % solve optimality LP
                    [obj,t,exitflagOpt] = isOptimal(obj,candidateActiveSet);

                    if exitflagOpt==1
                        % store in solution if optimal
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

                    % solve feasibility LP if not optimal
                    elseif exitflagOpt==-2
                        [obj,exitflagFeas] = isFeasible(obj,candidateActiveSet);

                        % add to set of set of infeasible active sets if infeasible
                        if exitflagFeas==-2
                            counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                            if ~mod(counter_infeasibleActiveSets,chunkSize)
                                infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                            end
                            infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidateActiveSet;
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

        function [obj,currentSolution,currentDegenerate] = solveWithPreviousSolution(obj,previousSolution,previousDegenerate,maxCardinality)
            %SOLVEWITHPREVIOUSSOLUTION determines the active sets in the
            % solution based on the solution for the previous horizon N-1
            % (Alg. 2 Mitze2020*, see also Mönnigmann2019**).
            %
            %
            % Input Arguments:
            %   previousSolution: matrix containing active sets of the
            %       solution for horizon N-1 (logical line vectors of 
            %       length obj.q-obj.qUX, inactive constraint = false, 
            %       active constraint = true)
            %
            %   previousDegenerate: logical vector. Every entry corresponds
            %       to a line in previousSolution and indicates true, if 
            %       the active set in that line has the solution t=0 to the
            %       LP implemented in isOptimal
            %
            %   maxCardinality: maximal cardinality of active sets.
            % 
            % Output Arguments:
            %   currentSolution: matrix containing active sets of the
            %   solution for horizon N (logical line vectors of length
            %   obj.q, inactive constraint = false, active constraint =
            %    true)
            %
            %   currentDegenerate: logical vector. Every entry corresponds
            %   to a line in currentSolution and indicates true, if the
            %   active set in that line has the solution t=0 to the
            %   optimality LP
            %
            %
            % *R. Mitze and M. Mönnigmann, "A dynamic programming
            %  approach to solving constrained linear-quadratic optimal
            %  control problems", Automatica, vol. 120, p. 109132, 2020.
            %
            % **M. Mönnigmann, "On the Structure of the Set of Active
            %  Sets," Automatica, vol. 106, pp. 61-69, 2019.
            %  doi.org/10.1016/j.automatica.2019.04.017. 
            
            % Initialize
            chunkSize = 1000; % number of preallocated lines (overestimated guess)
            currentSolution      = false(chunkSize,obj.q);        
            infeasibleActiveSets = false(chunkSize,obj.q);
            currentDegenerate    = false(chunkSize,1);
            counter_currentSolution      = 0;
            counter_infeasibleActiveSets = 0;

            % for every active set in solution for previous horizon
            for previousSolution_i=1:size(previousSolution,1)
                activeSet = previousSolution(previousSolution_i,:);
  
                % if all constraints corresponding to terminal stage are inactive
                if ~max(activeSet((obj.N-1)*obj.qUX+1:end))
                    % add to solution
                    counter_currentSolution = counter_currentSolution+1;
                    if ~mod(counter_currentSolution,chunkSize)
                        currentSolution   = [currentSolution; false(chunkSize,obj.q)]; % allocate more lines
                        currentDegenerate = [currentDegenerate; false(chunkSize,1)];   % allocate more lines
                    end
                    currentSolution(counter_currentSolution,:) = [activeSet false(1,obj.qUX)]; 
                    if previousDegenerate(previousSolution_i,:)
                        % set degeneracy flag according to degeneracy flag for previous horizon
                        currentDegenerate(counter_currentSolution,:) = true;
                    end           
                end

                % if any constraint corresponding to the terminal stage or the previous one is active
                if max(activeSet((obj.N-2)*obj.qUX+1:end))
                    
                    % extend active set by every combination of active/inactive constraints in stage 0
                    candidateStage0 = false(1,obj.qUX);
                    while true
                        candidate = [candidateStage0 activeSet];

                            % if no subset is detected infeasible active set
                            [obj,is_infeasible] = isPruned(obj,infeasibleActiveSets(1:counter_infeasibleActiveSets,:),candidate);
                            if ~is_infeasible
                                
                                % solve optimality LP
                                [obj,t,exitflagOpt] = isOptimal(obj,candidate);
                                
                                % store in solution if optimal
                                if exitflagOpt==1
                                    counter_currentSolution = counter_currentSolution+1;
                                    if ~mod(counter_currentSolution,chunkSize)
                                        currentSolution   = [currentSolution; false(chunkSize,obj.q)]; % allocate more lines
                                        currentDegenerate = [currentDegenerate; false(chunkSize,1)];   % allocate more lines
                                    end
                                    currentSolution(counter_currentSolution,:) = candidate;
                                    if t<=eps
                                        % set degeneracy flag according to degeneracy flag from previous horizon
                                        currentDegenerate(counter_currentSolution,:) = true;
                                    end

                                elseif exitflagOpt==-2
                                    % solve feasibility LP if not optimal
                                    [obj,exitflagFeas] = isFeasible(obj,candidate);
                                    if exitflagFeas==-2
                                        % add to set of infeasible active sets if infeasible
                                        counter_infeasibleActiveSets = counter_infeasibleActiveSets+1;
                                        if ~mod(counter_infeasibleActiveSets,chunkSize)
                                            infeasibleActiveSets = [infeasibleActiveSets; false(chunkSize,obj.q)]; % allocate more lines
                                        end
                                        infeasibleActiveSets(counter_infeasibleActiveSets,:) = candidate;
                                    end
                                end
                            end

                        % compute next active set in combinatorial tree for the extension in stage 0
                        maxCardinality_candidate = min(maxCardinality-sum(activeSet),obj.qUX);
                        [obj,candidateStage0,finished] = nextActiveSet(obj,candidateStage0,maxCardinality_candidate);
                        if finished
                            break
                        end
                    end
                end
            end
            currentSolution   = currentSolution(1:counter_currentSolution,:);   % delete unnecessary lines
            currentDegenerate = currentDegenerate(1:counter_currentSolution,:); % delete unnecessary lines
        end
                
        function [Ot,ct] = determineTerminalSet(obj)
            %DETERMINETERMINALSET determines the largest possible set such
            %  that the optimal feedback for the unconstrained infinite
            %  horizon problem stabilizes the system without violating the
            %  constraints with the approach from Gilbert1991*. The 
            %  terminal set is given in minimal H-representation using
            %  Polyhedron from the multi-parametric toolbox**.
            %
            %
            % Output Arguments:
            %   Ot and ct such that terminal set is {x | Ot*x<=ct}
            %
            %
            % *E. G. Gilbert and K. T. Tan, "Linear Systems with State and
            %  Constrol Constraints: The Theory and Application of Maximal
            %  Output Admissible Sets", IEEE Transactions on Automatic
            %  Control, vol. 36 (9), pp. 1008-1020, 1991.
            %
            % **M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
            %  Parametric Toolbox 3.0," Proc. of the European Control
            %  Conference, 2013, pp. 502-510.
            %
            %
            % See also Polyhedron and linprog and dlqr
            
            % test if MPT is in path
            if ~exist('Polyhedron.m','file') && ~exist('Polyhedron.m','class')
                error('determineTerminalSet uses functions from Multi-Parametric Toolbox. Toolbox is not in path.')
            end
            
            % LQR
            K = (-1)*dlqr(obj.A,obj.B,obj.Q,obj.R,zeros(obj.n,obj.m));

            % set options for linprog
            options = optimoptions('linprog','Display','none');

            % combine state and input constraints for x and K*x, respectively
            OUX = [obj.Ox(find(abs(obj.cx)<inf),:); obj.Ou(find(abs(obj.cu)<inf),:)*K]; 
            cUX = [obj.cx(find(abs(obj.cx)<inf)); obj.cu(find(abs(obj.cu)<inf))];

            % initialize
            k  = 0;
            Ot = OUX;
            ct = cUX;

            % increase time step
            while true
                k = k+1;
                flag = false;
                
                % insert linear system and LQR such that xk = calA*x0
                calA = (obj.A+obj.B*K)^k;
                
                % for each constraint in the set of the state and input constraints
                for i=1:obj.qUX
                    
                    % determine constraint for current time step (with xk)
                    Oi = OUX(i,:)*calA;
                    ci = cUX(i);
                    
                    % test if constraint is redundant to current terminal set
                    % (does there exist a point xk, where Ot*x<=ct, that violates Oi*xk<=ci?)
                    [xk,~,exitflag,~] = linprog(-Oi',Ot,ct,[],[],[],[],options);
                    if exitflag<0
                        error(['Exitflag of LP in determineTerminalSet is ' num2str(exitflag)])
                        
                    % if is not redundant, add constraint to current terminal set
                    elseif Oi*xk>ci-1e-6  % !! tolerance manually chosen !!
                        flag = true;
                        Ot   = [Ot; Oi];
                        ct   = [ct; ci];
                    end
                end
                
                % break if no constraint was added in current time step
                if ~flag
                    break
                end
            end
            
            % compute minimal representation of terminal set
            PolyT   = Polyhedron(Ot,ct);
            [~,sol] = PolyT.minHRep();
            Ot      = sol.H(:,1:end-1);
            ct      = sol.H(:,end);
        end
        
        function plotPolytopes(obj)
            %PLOTPOLYTOPES plots polytopes defined by the active sets in
            %  obj.solution in current figure using Polyhedron from the
            %  multi-parametric toolbox*. Polytopes definded by active
            %  sets with active terminal constraints are shown in gray,
            %  polytopes with inactive constraints in stages N and N-1 are
            %  shown in white, other polytopes are shown in blue.
            %
            %
            % *M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari "Multi-
            %  Parametric Toolbox 3.0," Proc. of the European Control
            %  Conference, 2013, pp. 502-510.
            %
            % 
            % See also Polyhedron
            
            fprintf('\nPolytopes definded by active sets with active terminal constraints are shown in gray,\npolytopes definded by active sets with inactive constraints in stages N and N-1 are shown in white,\nother polytopes are shown in blue.\n')
            
            % for all active sets in solution
            for activeSet_i=1:size(obj.solution,1)
                [T,d] = computePolytope(obj,obj.solution(activeSet_i,:));
                Poly  = Polyhedron(T,d);
                
                % plot polytope defined by active set with only inactive constraints in stages N and N-1 in white color
                if ~any(obj.solution(activeSet_i,end-obj.qT-obj.qUX+1:end))
                    plot(Poly,'Color',[1 1 1])
                
                % plot polytope defined by active set with active terminal constraints in grey color
                elseif any(obj.solution(activeSet_i,end-obj.qT+1:end))
                    plot(Poly,'Color',[.8 .8 .8])
                
                % otherwise, plot polytope defined by active set in cyan color
                else
                    plot(Poly,'Color','cyan')
                end
            end
        end
    end
        
end