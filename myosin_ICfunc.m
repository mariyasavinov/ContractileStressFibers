function [m_func, dmdx_func, dmdy_func] = myosin_ICfunc(IC_myosin,La,Lb,Lc,Ld)

if IC_myosin==0
    %%
    m_func = @(x,y) zeros(size(x)).*zeros(size(y));
    dmdx_func = @(x,y) zeros(size(x)).*zeros(size(y));
    dmdy_func = @(x,y) zeros(size(x)).*zeros(size(y));
elseif IC_myosin == 1
    %%
    % peak location
    gx1 = (Lb-La)/2+La;   gy1 = (Ld-Lc)/2+Lc;
    % peak spread
    sx1 = (Lb-La)/10;     sy1 = (Ld-Lc)/10;
    
    m_func = @(x,y) exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2);
    dmdx_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x);
    dmdy_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y);

elseif IC_myosin == 2
    %%
    % peak location
    gx1 = (Lb-La)*1/4+La;    gy1 = (Ld-Lc)*1/4+Lc;
    % peak spread
    sx1 = (Lb-La)/10;      sy1 = (Ld-Lc)/10;
    
    m_func = @(x,y) exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2);
    dmdx_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x);
    dmdy_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y);

elseif IC_myosin == 3
    %%
    % peak 1 location and spread
    gx1 = (Lb-La)*1/4+La;    gy1 = (Ld-Lc)*3/4+Lc;
    sx1 = (Lb-La)/10;      sy1 = (Lb-La)/10;
    % peak 2 location and spread
    gx2 = (Lb-La)*3/4+La;  gy2 = (Ld-Lc)/4+Lc;
    sx2 = (Lb-La)/10;      sy2 = (Lb-La)/10;
    
    m_func = @(x,y) exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2)+exp(1).^((-1).*sx2.^(-2).*((-1).*gx2+x).^2+(-1).* ...
      sy2.^(-2).*((-1).*gy2+y).^2);

    dmdx_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x)+(-2).*exp(1).^((-1).* ...
      sx2.^(-2).*((-1).*gx2+x).^2+(-1).*sy2.^(-2).*((-1).*gy2+y).^2).* ...
      sx2.^(-2).*((-1).*gx2+x);

    dmdy_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y)+(-2).*exp(1).^((-1).* ...
      sx2.^(-2).*((-1).*gx2+x).^2+(-1).*sy2.^(-2).*((-1).*gy2+y).^2).* ...
      sy2.^(-2).*((-1).*gy2+y);
elseif IC_myosin == 4
    %%
    m_func = @(x,y) ones(size(x)).*ones(size(y));
    dmdx_func = @(x,y) zeros(size(x)).*zeros(size(y));
    dmdy_func = @(x,y) zeros(size(x)).*zeros(size(y));
elseif IC_myosin == 5
    %%
    % peak 1 location and spread
    gx1 = (Lb-La)*3/4+La;    gy1 = (Ld-Lc)*3/4+Lc;
    sx1 = (Lb-La)/10;      sy1 = (Lb-La)/10;
    
    % peak 2 location and spread
    gx2 = (Lb-La)*3/4+La;  gy2 = (Ld-Lc)/4+Lc;
    sx2 = (Lb-La)/10;      sy2 = (Lb-La)/10;
    
    m_func = @(x,y) exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2)+exp(1).^((-1).*sx2.^(-2).*((-1).*gx2+x).^2+(-1).* ...
      sy2.^(-2).*((-1).*gy2+y).^2);

    dmdx_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x)+(-2).*exp(1).^((-1).* ...
      sx2.^(-2).*((-1).*gx2+x).^2+(-1).*sy2.^(-2).*((-1).*gy2+y).^2).* ...
      sx2.^(-2).*((-1).*gx2+x);

    dmdy_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y)+(-2).*exp(1).^((-1).* ...
      sx2.^(-2).*((-1).*gx2+x).^2+(-1).*sy2.^(-2).*((-1).*gy2+y).^2).* ...
      sy2.^(-2).*((-1).*gy2+y);
  
elseif IC_myosin == 6
    %%
    p = 6;
    m_func = @(x,y) (sin(pi*(x-La)/(Lb-La)).*sin(pi*(y-Lc)/(Ld-Lc))).^p;
    dmdx_func = @(x,y) ((-1).*La+Lb).^(-1).*p.*pi.*cos(((-1).*La+Lb).^(-1).*pi.*((-1).* ...
      La+x)).*sin(((-1).*Lc+Ld).^(-1).*pi.*((-1).*Lc+y)).*(sin(((-1).* ...
      La+Lb).^(-1).*pi.*((-1).*La+x)).*sin(((-1).*Lc+Ld).^(-1).*pi.*(( ...
      -1).*Lc+y))).^((-1)+p);
    dmdy_func = @(x,y) ((-1).*Lc+Ld).^(-1).*p.*pi.*cos(((-1).*Lc+Ld).^(-1).*pi.*((-1).* ...
      Lc+y)).*sin(((-1).*La+Lb).^(-1).*pi.*((-1).*La+x)).*(sin(((-1).* ...
      La+Lb).^(-1).*pi.*((-1).*La+x)).*sin(((-1).*Lc+Ld).^(-1).*pi.*(( ...
      -1).*Lc+y))).^((-1)+p);


elseif IC_myosin == 7
    %%
    % peak location
    gx1 = (Lb-La)*3/4+La;    gy1 = (Ld-Lc)*1/4+Lc;
    % peak spread
    sx1 = (Lb-La)/10;      sy1 = (Ld-Lc)/10;
    
    m_func = @(x,y) exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2);
    dmdx_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x);
    dmdy_func = @(x,y) (-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y);
elseif IC_myosin == 8
    %%
    % peak location
    gx1 = (Lb-La)/2+La;   gy1 = (Ld-Lc)/2+Lc;
    % peak spread
    sx1 = (Lb-La)/10;     sy1 = (Ld-Lc)/10;
    
    m_func = @(x,y) 1/2*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2);
    dmdx_func = @(x,y) 1/2*(-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x);
    dmdy_func = @(x,y) 1/2*(-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y);


elseif IC_myosin == 9
    %%
    % peak location
    gx1 = (Lb-La)*1/4+La;    gy1 = (Ld-Lc)*1/4+Lc;
    % peak spread
    sx1 = (Lb-La)/10;      sy1 = (Ld-Lc)/10;
    
    m_func = @(x,y) 1/2*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).*((-1) ...
      .*gy1+y).^2);
    dmdx_func = @(x,y) 1/2*(-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sx1.^(-2).*((-1).*gx1+x);
    dmdy_func = @(x,y) 1/2*(-2).*exp(1).^((-1).*sx1.^(-2).*((-1).*gx1+x).^2+(-1).*sy1.^(-2).* ...
      ((-1).*gy1+y).^2).*sy1.^(-2).*((-1).*gy1+y);
end