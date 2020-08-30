classdef GifSaver < handle
  properties
    data;
    colrmap;
    filename;
    delay_time;
    H;%graphic handle
  end
  methods
    function [self] = GifSaver(H,filename)
      if nargin<2 || isempty(filename)
        filename = ['untitled_' datestr(now,'yyyy-mm-dd_HH-MM-SS')];
      end
      self.H = H;
      self.filename = filename;
      self.colrmap  =[];
      self.data = uint8([]);
      self.delay_time = .01;
    end
  
    function captureFrame(self)
      h = gcf();
      changed_figures = false;
      if h.Number~=self.H.Number
        figure(self.H);
        changed_figures = true;
      end
      rgb = frame2im(getframe(self.H));
      [imind,clrmap] = rgb2ind(rgb,256); %#ok<ASGLU>
      if false
        if isempty(self.colrmap) %#ok<UNRCH>
          self.colrmap = clrmap;
          self.data(:,:,:,1) = imind;
        else
          self.data(:,:,end+1,:) = imind;
        end
      else
        if isempty(self.data)
          self.data(:,:,:,1) = rgb;
        else
          self.data(:,:,:,end+1) = rgb;
        end
      end
      if changed_figures
        figure(h);
      end
    end
    
    function clear(self)
      self.colrmap = [];
      self.data = [];
    end
    
    function save(self)
      if isempty(self.data)
        warning('Attempting to save GIF with no data');
        return;
      end
      if ~isempty(self.colrmap)
        imwrite(self.data,self.colrmap,[self.filename '.gif'],'gif','DelayTime',self.delay_time,'Loopcount',inf);
      else
        %smush a subsampled all frames to get a cohesive colormap
        tmp = self.data(:,:,:,1);
        for i=floor(linspace(1,size(self.data,4),150))
          tmp = cat(1,tmp,self.data(:,:,:,i));
        end 
        [~,self.colrmap] = rgb2ind(tmp,240);
        iminds = repmat(uint8(0),[size(self.data,1),size(self.data,2),1,size(self.data,4)]);
        for i=1:size(self.data,4)
          [iminds(:,:,1,i),~] = rgb2ind(self.data(:,:,:,i),self.colrmap);
        end
        imwrite(iminds,self.colrmap,[self.filename '.gif'],'gif','DelayTime',self.delay_time,'Loopcount',inf);
      end
    end
  end
end