% Code by David Ungarish.

classdef Subplots < handle
    % A subplot handler, to be used instead of subplot().
    % Main features:
    %
    % - Subplot grid is defined once, instead of at each axis creation.
    %       e.g. sp = Subplots(ncol, nrows)
    %
    % - Option to auto-determine grid dimensions by giving only the total
    %   number of axes.
    %       e.g. sp = Subplots(n, -1) % n axis grid, with #rows >= #cols
    %            sp = Subplots(-1, n) % n axis grid, with #rows <= #cols
    %
    % - Active axis can be specified in several ways:
    %   1. "Flat" index, i.e. sp.axis(idx)
    %   2. Sub-index, i.e. sp.axis([row, col])
    %   3. Determined automatically according to axis order, i.e. sp.axis()
    %
    % - The axis handles are saved to an array for later access.
    %
    % Usage:
    %
    %	1. Create a Subplots instance (in current figure):
    %       
    %       figure();
    %       sp = Subplots(n, m)
    %       
    %       * If n == -1, m is regarded as the desired total number of
    %           axes, and a grid is created in size NxM such that N*M == m,
    %           and N>=M where N & M are as close as possible. I.e. the
    %           grid is as close as possible to being square, under the
    %           constraint that N*M must equal the total number of axis
    %           (or greater than, if number of axis is prime).
    %
    %       * If m == -1, the same applies, but with N<=M.
    %
    %   2. Use axis() method to navigate between axes. The axis to activate
    %       can be specified in three ways:
    %       
    %       sp.axis(i)          % activate i-th axis (flat index notation)
    %       sp.axis([j, k])     % activate axis j,k (sub index notation)
    %       sp.axis()           % activate next axis in the grid.
    %
    %   3. To span several axes, use axis() with a column array. i.e.
    %       sp.axis([i1;i2;i3]) will span axis i1,i2,i3
    %       similarly, sp.axis([j1,k1;j2,k2;j3,k3]) will span 3 axes, using
    %       sub-index notation.
    %
    %   * Note that by default, axes are created only when they are first
    %       activated. To create all axes in advance, use: sp.make()
    %
    % Examples:
    %   
    %   figure();
    %   sp = Subplots(5,3); % initialize 5 * 3 subplot axis
    %   sp.make()  % optional, create axis in advance (usually not needed)
    %
    %   sp.axis([1,2]); % focus on axis at row=1 col=2
    %   plot(1:10,'r')
    %   sp.axis(7); % focus on 7th axis (in this case row=3 col=1)
    %   plot(1:10, 'b')
    %   sp.axis() % focus on the next axis (in this case: row=3 col=2)
    %   plot(1:10,'g')
    %   sp.axis([2,3;2,4]) % span axis on axis (2,3) & (2,4)
    %
    %   Alternative ways to initalize:
    %
    %   sp = Subplots(12, -1)   % auto-determine the grid dimensions. in
    %                           % this case, the result will be a 3*4 grid.
    %
    %   sp = Subplots(-1, 12)   % auto-determine the grid dimensions. in
    %                           % this case, the result will be a 4*3 grid.
    %
    %   sp = Subplots(12)       % the same as Subplots(12, -1)
    %
    %   For loop example:
    %
    %       Draw 6 plots in a 2x3 grid, then set 'box' to 'on' in all of
    %       them.
    %
    %   figure();
    %   sp = Subplots(2,3);
    %   for i = 1 : 6
    %       sp.axis();
    %       plot(rand(1,10));
    %   
    %   set(sp.handles,'box','on');

    properties
        gridsize    % size of grid
        handles     % handles to axis
        parentfig   % handle to parent figure
        lastpos     % index of axis last plotted to
    end
    
    methods
    
        function obj = Subplots(nrows, ncols)
            % Constructor
            
            if ~exist('ncols','var'), ncols = -1; end
            assert(nrows > 0 || ncols > 0, "Cannot determine grid size")
            
            if ncols == -1 || nrows == -1
                % Determine grid size for a given number of axis. Try to
                % make the grid as sqaure as possible.

                n = max([nrows, ncols]) ;
                ms = 1:ceil(sqrt(n));
                err = n - ms.*ceil(n./ms);
                df = ms-ceil(n./ms);
                [~,mni] = min(2*abs(err)+abs(df));
                
                dim1 = ms(mni);
                dim2 = ceil(n/dim1);
                
                if ncols == -1
                    nrows = dim1;
                    ncols = dim2;
                else
                    nrows = dim2;
                    ncols = dim1;
                end
                
            end
            
            obj.parentfig = gcf ;
            obj.gridsize = [nrows,ncols];
            obj.handles = zeros(obj.gridsize);
            obj.lastpos = 0;
            
        end
        
        function [h,indx] = axis(obj,gridpos)
            % Create axis in grid position [gridpos]
            % Input:
            %   Indices of grid to plot to. Can be either a sub or an ind.
            %   Stack indices in a coloumn to provide multiple indices,i.e:
            %       [ i(1),j(1) ; i(2),j(2) ; i(3),j(3) , .. ]
            %       OR:
            %       [ ind(1) ; ind(2) .. ]
            %   If no input is given, the next grid in the list will be
            %   selected.
            %
            % Output:
            %   h - handle to axis
            %   indx - index of axis
            
            
            % **
            % convert input grid-pos to grid-pos thats supported by
            % matlab's subplot() :
            
            if ~exist('gridpos','var')
                % no position is given- increment last position:
                indx = obj.lastpos + 1;
                
            elseif size(gridpos,2)==1
                % position is given as an index (not subindex)
                indx = gridpos;
                
            else
                % position is given as a list of subindices,
                % convert to inds:
                gridsz = obj.gridsize([2,1]) ;
                indx = nan(size(gridpos,1),1);
                for i = 1 : size(gridpos,1)
                    subindx = gridpos(i,:);
                    indx(i) = sub2ind(gridsz,subindx(2),subindx(1));
                end
                
            end
            
            % **
            % get axis handle:
            
            figure(obj.parentfig);
            h = subplot(obj.gridsize(1),obj.gridsize(2),indx);

            obj.handles(indx) = h ; % append to list of axis handles
            obj.lastpos = max(indx); % keep track of position in grid

        end
        
        function r = nAxes(obj)
            % Number of axes
            r = numel(obj.handles);
        end
        
        function make(obj)
            % make the axes
            
            inds = find(obj.handles==0);
            for ind = inds(:)'
                obj.axis(ind);
            end
            obj.lastpos = 0 ;
        end
        
    end
    
end

