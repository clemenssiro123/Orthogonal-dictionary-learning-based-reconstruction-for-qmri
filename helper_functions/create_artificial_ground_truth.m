function [T1,T2,Density] = create_artificial_ground_truth(complexity,image_size)
%CREATE_ARTIFICIAL_GROUND_TRUTH Simple qMRI test images with regions and texture.
%   [T1,T2,Density] = create_artificial_ground_truth(complexity,image_size)
%   returns double arrays in physical units: T1/T2 in ms and Density in the
%   demo's intensity units. This is a synthetic phantom, not anatomy.
%
%   complexity: nonnegative integer controlling BOTH regions and texture.
%       0 = uniform object; 5 = a few regions with gentle intensity variation;
%       30 = many smaller regions with stronger, finer texture inside them.
%       The number of added ellipses equals complexity. Ellipses can overlap;
%       later ones overwrite earlier ones. Texture is added afterwards.
%       This setting is not a measured image-complexity score.
%   image_size: side length, e.g. 256, or [rows columns], e.g. [128 192].
%   Defaults: complexity=10, image_size=256. Small details may disappear
%   on coarse grids or after subsampling with the demo's pixel_stride.
%
%   Example (from the project folder, after adding helper_functions):
%       [T1,T2,Density] = create_artificial_ground_truth(10,256);
%       save('data/artificial_maps.mat','T1','T2','Density');
%   Then set the demo's input_file to data/artificial_maps.mat.
%   Use pixel_stride=2 to reconstruct these 256-by-256 maps at 128-by-128.
%   All values fit the demos' default bounds; no toolbox is required.
%   Texture belongs to the clean ground truth; measurement noise is added by the demos.

    %% 1. Choose the complexity and image grid
    if nargin<1, complexity=10; end
    if nargin<2, image_size=256; end
    validateattributes(complexity,{'numeric'}, ...
        {'real','finite','scalar','integer','nonnegative'});
    validateattributes(image_size,{'numeric'}, ...
        {'real','finite','vector','integer','>=',3});
    if isscalar(image_size), image_size=[image_size image_size]; end
    assert(numel(image_size)==2,'image_size must be a side length or [rows columns].');
    % Use double arithmetic even if settings were supplied as integer arrays.
    complexity=double(complexity); image_size=double(image_size(:).');
    % Normalized coordinates keep the same geometry when resolution changes.
    [x,y]=meshgrid(linspace(-1,1,image_size(2)),linspace(-1,1,image_size(1)));

    %% 2. Create one uniform elliptical object in a zero-density background
    body=(x/.85).^2+(y/.95).^2<=1;
    % Positive background times keep Bloch calculations valid; zero density
    % makes the background signal zero regardless of these relaxation times.
    T1=1000*ones(image_size); T2=100*ones(image_size); Density=zeros(image_size);
    T1(body)=900; T2(body)=80; Density(body)=80;
    % Edit this small table to change tissue contrasts: each row is [T1 T2 Density].
    tissues=[600 40 55; 1400 110 90; 2600 200 75; 4500 450 100];

    %% 3. Add more and smaller tissue regions as complexity increases
    % A private fixed seed makes images repeatable without changing the
    % demo's random-number sequence used later for measurement noise.
    stream=RandStream('mt19937ar','Seed',1);
    radius=.4/sqrt(1+complexity/5);
    for region=1:complexity
        % Place the center inside the object using a random angle and radius.
        angle=2*pi*rand(stream); distance=.7*sqrt(rand(stream));
        center_x=.85*distance*cos(angle); center_y=.95*distance*sin(angle);
        % Vary ellipse width, height and orientation while retaining sharp edges.
        a=radius*(.6+.4*rand(stream)); b=radius*(.4+.3*rand(stream));
        rotation=2*pi*rand(stream);
        xr=(x-center_x)*cos(rotation)+(y-center_y)*sin(rotation);
        yr=-(x-center_x)*sin(rotation)+(y-center_y)*cos(rotation);
        inside=body & ((xr/a).^2+(yr/b).^2<=1);
        % Use the same region in all maps, cycling through the tissue table.
        tissue=1+mod(region-1,size(tissues,1));
        T1(inside)=tissues(tissue,1);
        T2(inside)=tissues(tissue,2);
        Density(inside)=tissues(tissue,3);
    end

    %% 4. Add smooth texture within the tissue while keeping its sharp boundaries
    if complexity>0
        % Higher complexity increases the spatial frequency. Cap it so the
        % finest wave still has at least four pixels per cycle on this grid.
        frequency=min(2+complexity/2,(min(image_size)-1)/8);
        % Three differently oriented waves give a smooth, repeatable texture
        % in [-1,1]. These are spatial variations, not independent pixel noise.
        texture=(sin(pi*frequency*x).*cos(pi*frequency*y) ...
            +sin(pi*frequency*(.6*x+.8*y)+1) ...
            +sin(2*pi*frequency*(.8*x-.6*y)+2))/3;
        % Gradually increase contrast, with less than 8% variation around
        % each tissue value. This keeps the default tissue table within bounds.
        strength=.08*complexity/(complexity+5);
        factor=1+strength*texture(body);
        % Apply matching texture in all three maps; leave the background unchanged.
        T1(body)=T1(body).*factor;
        T2(body)=T2(body).*factor;
        Density(body)=Density(body).*factor;
    end
end
