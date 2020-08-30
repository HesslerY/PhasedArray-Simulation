%#ok<*DEFNU,*UNRCH>

%clear all; %#ok<CLALL>

f_Hz = 10e9;
lambda_m = Constants.speed_of_light/f_Hz;
NUM_ELEMENTS = 888;

PLOT_SINGLE = true;
PLOT_ANIMATE = true;
SAVE_ANIMATED_GIF = true;

%Makes a grid
% elements = createCurvedGrid(NUM_ELEMENTS,lambda_m);
elements = createFlatGrid(NUM_ELEMENTS,lambda_m);
elements = steerElements(elements,[-1;.3;0],lambda_m);

if ~exist('tmp','dir')
  mkdir('tmp')
end

if PLOT_SINGLE
  figure(1);clf;
  directivity = plotSpatialPattern(elements,lambda_m,215);
  fprintf('Directivity %.2f dB\n',10*log10(directivity));
end

if PLOT_ANIMATE
  H=figure(2);clf;
  if SAVE_ANIMATED_GIF
    filename = fullfile('tmp',datestr(now,'yyyymmdd-HHMMSS'));
    gs = GifSaver(H,filename);
    gs.delay_time = 1/24;
  end
  thetas = linspace(-pi/4,2*pi-pi/4,151);
  %thetas = linspace(-pi/3,pi/3,75);
  %thetas = [thetas,fliplr(thetas)];
  for n=1:numel(thetas)
    th = thetas(n);
    phi = 20*th;
    psi = .2;
    u = [-cos(th);sin(th);0]; %spins around up axis
    % u = [cos(th),sin(th),0;-sin(th),cos(th),0;0,0,1]*[1,0,0;0,cos(phi),sin(phi);0,-sin(phi),cos(phi)]*[-cos(psi);0;sin(psi)];
    elements = steerElements(elements,u,lambda_m);
    directivity = plotSpatialPattern(elements,lambda_m,215);
    fprintf('Steer angle %.1f deg, Directivity %.2f dB\n',th*180/pi,10*log10(directivity));
    if SAVE_ANIMATED_GIF
      gs.captureFrame();
    else
      pause(1/100);
    end
  end
  if SAVE_ANIMATED_GIF
    gs.save();
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [directivity] = plotSpatialPattern(elements,lambda_m,NUM_POINTS)
  if nargin<3 || isempty(NUM_POINTS)
    NUM_POINTS = 100;
  end
  az = linspace(-pi,pi,NUM_POINTS);
  el = asin(linspace(-1,1,NUM_POINTS));
  [azaz,elel] = meshgrid(az,el);
  ned_unit = getNedFromAzEl(azaz(:),elel(:));

  targets = ned_unit;
  pattern = calculatePattern(elements,targets,lambda_m);
  pattern = reshape(pattern,size(azaz));
  
  P = abs(pattern).^2;
  solid_angle = trapz(az,trapz(el,cos(elel).*P,1));
  max_value = max(P,[],'all');
  average_P = (solid_angle)/(4*pi);
  directivity = max_value/(average_P);
  
  FC = 10*log10(P/(average_P));
  FC = clamp(FC,max(FC,[],'all')-60,inf);

  axs = [];
  
  subplot(3,1,1);
  plotElementLocations([elements.position],[elements.phase],lambda_m);
  axs(end+1) = gca();

  subplot(3,1,2);
  plotPatternOnSphere(targets,FC);
  axs(end+1) = gca();

  subplot(3,1,3);
  plotPatternWithHeight(targets,FC);
  axs(end+1) = gca();
  
  Link = linkprop(axs,{'CameraUpVector', 'CameraPosition', 'CameraTarget'});
  setappdata(gcf, 'StoreTheLink', Link);
end

function plotElementLocations(positions,phases_rad,lambda_m)
  if nargin<3 || isempty(lambda_m)
    lambda_m = [];
  end
  if nargin<2 || isempty(phases_rad)
    phases_rad = [];
  end
  if ~isempty(lambda_m)
    positions = positions/lambda_m;
  end
  if isempty(phases_rad)
    plot3(positions(2,:),positions(1,:),-positions(3,:),'b.');
  else
    phases_rad = mod(phases_rad,2*pi);
    uniform_phases = linspace(0,2*pi,numel(phases_rad));
    idx = arrayfun(@(c)argmin(abs(c-uniform_phases)),phases_rad);
    colors = hsv(numel(phases_rad));
    colors = colors(idx,:);
    csize = clamp(15-log(numel(phases_rad)),5,20);
    scatter3(positions(2,:),positions(1,:),-positions(3,:),csize,colors,'filled');
    h=colorbar;
    h.Label.String='Element Phase';
  end
  axis('square');
  axis('equal');
  colorbar;
  grid('on');
  colorbar('off');
  s = tern(isempty(lambda_m),'','/\lambda');
  xlabel(['\bfEast [m',s,']']);
  ylabel(['\bfNorth [m',s,']']);
  zlabel(['\bfUp [m',s,']']);
  h=title(['Element Locations' tern(isempty(phases_rad),'',' with Phases')]);
  h.FontSize=14;
  view([7,14]);
  hold('off');
end

function plotPatternOnSphere(directions,pattern)
  nn = reshape(directions(1,:),size(pattern));
  ee = reshape(directions(2,:),size(pattern));
  dd = reshape(directions(3,:),size(pattern));
  surf(ee,nn,-dd,pattern);
  axis('square');
  h=colorbar;
  h.Label.String = '\bf20*log_{10}|G_{pattern}|';
  h.Label.FontSize=14;
  grid('on');
  xlabel('\bfEast');
  ylabel('\bfNorth');
  zlabel('\bfUp');
  shading('interp');
  h=title('Pattern Shape');h.FontSize=14;
  view([7,14]);
  hold('off');
end

function plotPatternWithHeight(directions,pattern)
  axs1 = axis();
  nn = reshape(directions(1,:),size(pattern));
  ee = reshape(directions(2,:),size(pattern));
  dd = reshape(directions(3,:),size(pattern));
  R = pattern-min(pattern,[],'all');
  surf(ee.*R,nn.*R,-dd.*R,pattern);
  axis('square');
  h=colorbar;
  h.Label.String = '\bf20*log_{10}|G_{pattern}|';
  h.Label.FontSize=14;
  grid('on');
  xlabel('\bfEast');
  ylabel('\bfNorth');
  zlabel('\bfUp');
  shading('interp');
  h=title('Pattern Shape');h.FontSize=14;
  view([7,14]);
  axs2 = axis();
  axs3 = combineAxes(axs1,axs2);
  axis(axs3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [element] = getElement(dim)
  if nargin<1 || isempty(dim)
    dim = 3;
  end
  element = struct('position',zeros(dim,1),'amplitude',1,'phase',0);
end

function [pattern] = calculatePattern(elements,targets,lambda_m)
  positions = [elements.position];
  %path_diff = nan(size(targets,2),size(positions,2));
  %for i=1:size(targets,2)
  %  for j=1:size(positions,2)
  %    path_diff(i,j) = sum(targets(:,i).*positions(:,j),1);
  %  end
  %end
  % same above as below, just fancier/quicker because it uses broadcasting
  path_diff = squeeze(sum(targets.*permute(positions,[1,3,2]),1));
  phase = 2*pi/lambda_m*path_diff+[elements.phase];
  phasor = sum([elements.amplitude].*exp(1j*phase),2);
  pattern = phasor.';
end

function [elements] = steerElements(elements,target,lambda_m)
  positions = [elements.position];
  target = target./rootSumSq(target);
  path_diff = squeeze(sum(target.*permute(positions,[1,3,2]),1));
  phase = 2*pi/lambda_m*path_diff;
  for n=1:numel(elements)
    elements(n).phase=-phase(n);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elements] = createFlatGrid(num_elements,lambda_m)
  de = lambda_m/sqrt(2);
  du = lambda_m/sqrt(2);
  Ne = floor(sqrt(num_elements));
  Nu = floor(num_elements/Ne);
  elements = reshape(arrayfun(@(n)ElementRadiator(),1:(Nu*Ne)),[Nu,Ne]);
  for ne=1:Ne
    for nu=1:Nu
      position = [...
        0;
        de*(ne-1);... %columns
        du*(nu-1)+mod(ne,2)*(du/2) ... %offset columns
      ];
      elements(nu,ne).position = position;
    end
  end
  center = mean([elements.position],2);
  for n=1:numel(elements)
    elements(n).position = elements(n).position-center;
  end
end

function [elements] = createCurvedGrid(num_elements,lambda_m)
  de = lambda_m/3;
  dn = lambda_m/2;
  du = lambda_m/2;
  Ne = floor(sqrt(num_elements));
  Nu = floor(num_elements/Ne);
  elements = reshape(arrayfun(@(n)ElementRadiator(),1:(Nu*Ne)),[Nu,Ne]);
  for ne=1:Ne
    for nu=1:Nu
      position = [...
        dn*sqrt(Ne*Nu)/10*(sin(pi*(ne-1)/Ne)+sin(pi*(nu-1)/Nu));... %curved thing
        de*(ne-1);... %columns
        du*(nu-1)+mod(ne,2)*(du/2) ... %offset columns
      ];
      elements(nu,ne).position = position;
    end
  end
  center = mean([elements.position],2);
  for n=1:numel(elements)
    elements(n).position = elements(n).position-center;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wuv] = getWuvFromNed(ned)
  R = sqrt(sum(abs(ned).^2,1));
  ned_unit = ned./R;
  wuv = [1;1;-1].*ned_unit;
end

function [ned_unit] = getNedFromWuv(wuv)
  ned_unit = [1;1;-1].*wuv;
end

function [az,el,varargout] = getAzElFromNed(ned) 
  wuv = getWuvFromNed(ned);
  az = atan2(wuv(2,:),wuv(1,:));
  el = asin(wuv(3,:));
  varargout = {};
  if nargout>=3
    varargout{1} = wuv(1,:)>=0;
  end
end

function [wuv] = getWuvFromAzEl(az,el,is_front)
  N = max([numel(az),numel(el)]);
  if nargin<2 || isempty(is_front)
    is_front = true(1,N);
  end
  el = el(:).';
  az = az(:).';
  az = az+pi.*(~is_front);
  wuv = [...
    cos(el).*cos(az);...
    cos(el).*sin(az);...
    sin(el)     ...
 ];
end

function [ned_unit] = getNedFromAzEl(az,el,is_front)
  if nargin<3 || isempty(is_front)
    is_front = [];
  end
  wuv = getWuvFromAzEl(az,el,is_front);
  ned_unit = getNedFromWuv(wuv);
end

function [axs3] = combineAxes(axs1,axs2)
  axs3=[];
  if isempty(axs2) || numel(axs1)>numel(axs2)
    axs3=axs1;
    return;
  end
  if isempty(axs1) || numel(axs2)>numel(axs1)
    axs3=axs2;
    return;
  end
  tmp = [reshape(axs1,2,[]);reshape(axs2,2,[])];
  axs3 = reshape([min(tmp,[],1);max(tmp,[],1)],1,[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = argmin(varargin)
  [~,out] = min(varargin{:});
end

function [out] = argmax(varargin)
  [~,out] = max(varargin{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
