## Copyright (C) 2019 Petros Lazaridis
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- This function solves numerically the diffusion partial differential equation
##(Fick's second law) in time and concrete depth domain. -*- 
## @deftypefn {} {@var{retval} =} clorideExposure (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: peter <peter@LAPTOP-OG3Q24A2>
## Created: 2019-10-19
% This function solves numerically the diffusion partial differential equation
% (Fick's second law) in time and concrete depth domain, accordind to Life-365
%  2.2.1 manual (pages 17-18).

function [U, D, r] = clorideExposure (c, te, dx, dt, D_28, m, Temp,...
   maxSurfConc, maxSurfConcFirstYear)
  
  i = 1;
  nmbofConcreteSlices = c/dx;
  t28 = 28/30;
  t = [0:dt:12*te]; % t in months
  Dref = D_28*(t28./t).^m ;
  Dref(Dref > D_28 | Dref == inf) = D_28;
  Dref(t > 25*12) = Dref(t == 25*12);
  TempperTime = interp1( [0:12], Temp([end, 1:end]) + 273, mod(t, 12) );
  SurfConcentr = interp1( [0 maxSurfConcFirstYear, te]*12,...
  [0, maxSurfConc, maxSurfConc], t ) ;
  D = Dref.*exp ( 3500/8.5*(1/293 - 1./TempperTime) );
  
  te = 12*30*24*3600*te; %conversion to seconds
  dt = 30*24*3600*dt; %conversion to seconds
  
  U(:, 1) = [ SurfConcentr(1), zeros(1, nmbofConcreteSlices) ]';
  
  for t = dt:dt:te
    i++;
  
    r = D(i-1)*dt/2/dx^2;
  
    helpA = [1, ones( 1, nmbofConcreteSlices -1 )*(1+2*r), 1]';
    help1 = [ ones( 1, nmbofConcreteSlices )*r, 0]';
    A = [ -help1, helpA, -help1(end:-1:1) ];
    A = spdiags(A, [-1 0 1], nmbofConcreteSlices + 1, nmbofConcreteSlices + 1);

    helpB = [1, ones(1, nmbofConcreteSlices-1)*(1-2*r), 1]';
    B =  [ help1, helpB, help1(end:-1:1) ];
    B = spdiags(B, [-1 0 1], nmbofConcreteSlices + 1, nmbofConcreteSlices + 1);
  
    U(:, i) = A\( B*[SurfConcentr(i-1); U(2:end, i-1)] );
  endfor
  
  endfunction
