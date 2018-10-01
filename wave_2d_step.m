function [v1, s11, s12, v2, s21, s22] = wave_2d_step(v1, s11, s12, v2, s21, s22, dt, rk, dx, c, z, tau, rl, rr, rb, rt, weak, cdiss)

% wave_2d_step takes a time step for solving the 2D wave equation for two elastic blocks coupled
% by a frictional interface. The interface can support a maximum shear stress tau, otherwise slip
% occurs on the interface. Outer boundaries can be absorbing, free surface, or rigid by specifying
% a reflection coefficient (rl, rr, rb, rt for the left, right, top, and bottom boundaries, respectively)
%
% parameters: v1, s11, s12, v2, s21, s22 are the fields for each block at the start of the time step
% dt is time step, rk is a data structure holding the integration information
% dx is the spatial grid spacing, c is the wave speed, z is the shear impedance
% tau is the interface strength, rl, rr, rb, rt are the external boundary conditions
% boundary conditions can be imposed weakly (weak = true) or strongly (weak = false)
% optionally, the method can include artificial dissipation (with coefficient cdiss)
%
% note: can omit the cdiss parameter (cdiss = 0 by default), or cdiss and weak (weak = true by default)
%
% returns updated fields v1, s11, s12, v2, s21, s22

	if nargin == 16
		cdiss = 0;
		weak = 1;
	elseif nargin == 17
		cdiss = 0;
	end

    assert(dt > 0);
    assert(dx > 0);
    assert(c > 0);
    assert(z > 0);
    assert(cdiss >= 0);
    assert(tau > 0);
    assert(rl >= -1 & rl <= 1);
    assert(rr >= -1 & rr <= 1);
    assert(rb >= -1 & rb <= 1);
    assert(rt >= -1 & rt <= 1);
    
    dv1 = zeros(size(v1));
    dv2 = zeros(size(v2));
    ds11 = zeros(size(s11));
    ds12 = zeros(size(s12));
    ds21 = zeros(size(s21));
    ds22 = zeros(size(s22));

    for j=1:rk.nstages

        % calculate finite differences
        dv1 = rk.A(j)*dv1+dt/dx*((calc_diff(s11,1)+calc_diff(s12,2))/z*c+cdiss*(calc_diss(v1,1)+calc_diss(v1,2)));
        ds11 = rk.A(j)*ds11+dt/dx*(calc_diff(v1,1)*z*c+cdiss*(calc_diss(s11,1)+calc_diss(s11,2)));
		ds12 = rk.A(j)*ds12+dt/dx*(calc_diff(v1,2)*z*c+cdiss*(calc_diss(s12,1)+calc_diss(s12,2)));
        dv2 = rk.A(j)*dv2+dt/dx*((calc_diff(s21,1)+calc_diff(s22,2))/z*c+cdiss*(calc_diss(v2,1)+calc_diss(v2,2)));
        ds21 = rk.A(j)*ds21+dt/dx*(calc_diff(v2,1)*z*c+cdiss*(calc_diss(s21,1)+calc_diss(s21,2)));
		ds22 = rk.A(j)*ds22+dt/dx*(calc_diff(v2,2)*z*c+cdiss*(calc_diss(s22,1)+calc_diss(s22,2)));

        % apply outer boundaries if boundary conditions imposed weakly

        if weak
			% left boundary
        	vhat1 = 0.5*(1-rl)*(s11(1,:)/z+v1(1,:));
        	shat1 = 0.5*(1+rl)*(s11(1,:)+z*v1(1,:));
            dv1(1,:) = dv1(1,:)+dt/dx*c*2.*(vhat1-v1(1,:));
            ds11(1,:) = ds11(1,:)+dt/dx*c*2.*(shat1-s11(1,:));
			
			% right boundary
	        vhat2 = 0.5*(rr-1)*(s21(end,:)/z-v2(end,:));
	        shat2 = 0.5*(rr+1)*(s21(end,:)-z*v2(end,:));
            dv2(end,:) = dv2(end,:)+dt/dx*c*2.*(vhat2-v2(end,:));
            ds21(end,:) = ds21(end,:)+dt/dx*c*2.*(shat2-s21(end,:));
			
			% bottom boundary
        	vhat1 = 0.5*(1-rb)*(s12(:,1)/z+v1(:,1));
        	shat1 = 0.5*(1+rb)*(s12(:,1)+z*v1(:,1));
            dv1(:,1) = dv1(:,1)+dt/dx*c*2.*(vhat1-v1(:,1));
            ds12(:,1) = ds12(:,1)+dt/dx*c*2.*(shat1-s12(:,1));
        	vhat2 = 0.5*(1-rb)*(s22(:,1)/z+v2(:,1));
        	shat2 = 0.5*(1+rb)*(s22(:,1)+z*v2(:,1));
            dv2(:,1) = dv2(:,1)+dt/dx*c*2.*(vhat2-v2(:,1));
            ds22(:,1) = ds22(:,1)+dt/dx*c*2.*(shat2-s22(:,1));
			
			% top boundary
	        vhat1 = 0.5*(rt-1)*(s12(:,end)/z-v1(:,end));
	        shat1 = 0.5*(rt+1)*(s12(:,end)-z*v1(:,end));
            dv1(:,end) = dv1(:,end)+dt/dx*c*2.*(vhat1-v1(:,end));
            ds12(:,end) = ds12(:,end)+dt/dx*c*2.*(shat1-s12(:,end));
	        vhat2 = 0.5*(rt-1)*(s22(:,end)/z-v2(:,end));
	        shat2 = 0.5*(rt+1)*(s22(:,end)-z*v2(:,end));
            dv2(:,end) = dv2(:,end)+dt/dx*c*2.*(vhat2-v2(:,end));
            ds22(:,end) = ds22(:,end)+dt/dx*c*2.*(shat2-s22(:,end));

			% apply fault boundary conditions
            phi = (s11(end,:)-z*v1(end,:)+s21(1,:)+z*v2(1,:))/2.;
            shat3 = phi;
            shat3(abs(phi) > tau) = sign(phi(abs(phi) > tau))*tau;

            vhat3 = (shat3-s11(end,:))/z+v1(end,:);
            vhat4 = (-shat3+s21(1,:))/z+v2(1,:);

            dv1(end,:) = dv1(end,:)+dt/dx*c*2.*(vhat3-v1(end,:));
            ds11(end,:) = ds11(end,:)+dt/dx*c*2.*(shat3-s11(end,:));
            dv2(1,:) = dv2(1,:)+dt/dx*c*2.*(vhat4-v2(1,:));
            ds21(1,:) = ds21(1,:)+dt/dx*c*2.*(shat3-s21(1,:));
			
        end

        % update
        v1 = v1+rk.B(j)*dv1;
        s11 = s11+rk.B(j)*ds11;
        s12 = s12+rk.B(j)*ds12;
        v2 = v2+rk.B(j)*dv2;
        s21 = s21+rk.B(j)*ds21;
        s22 = s22+rk.B(j)*ds22;

        % if enforcing boundary conditions strongly, modify boundary data

        if ~weak
			% left boundary
        	vhat1 = 0.5*(1-rl)*(s11(1,:)/z+v1(1,:));
        	shat1 = 0.5*(1+rl)*(s11(1,:)+z*v1(1,:));
            v1(1,:) = vhat1;
            s11(1,:) = shat1;
			
			% right boundary
	        vhat2 = 0.5*(rr-1)*(s21(end,:)/z-v2(end,:));
	        shat2 = 0.5*(rr+1)*(s21(end,:)-z*v2(end,:));
            v2(end,:) = vhat2;
            s21(end,:) = shat2;
			
			% bottom boundary
        	vhat1 = 0.5*(1-rb)*(s12(:,1)/z+v1(:,1));
        	shat1 = 0.5*(1+rb)*(s12(:,1)+z*v1(:,1));
            v1(:,1) = vhat1;
            s12(:,1) = shat1;
        	vhat2 = 0.5*(1-rb)*(s22(:,1)/z+v2(:,1));
        	shat2 = 0.5*(1+rb)*(s22(:,1)+z*v2(:,1));
            v2(:,1) = vhat2;
            s22(:,1) = shat2;
			
			% top boundary
	        vhat1 = 0.5*(rt-1)*(s12(:,end)/z-v1(:,end));
	        shat1 = 0.5*(rt+1)*(s12(:,end)-z*v1(:,end));
            v1(:,end) = vhat1;
            s12(:,end) = shat1;
	        vhat2 = 0.5*(rt-1)*(s22(:,end)/z-v2(:,end));
	        shat2 = 0.5*(rt+1)*(s22(:,end)-z*v2(:,end));
            v2(:,end) = vhat2;
            s22(:,end) = shat2;

			% apply fault boundary conditions
        	phi = (s11(end,:)-z*v1(end,:)+s21(1,:)+z*v2(1,:))/2.;
			shat3 = phi;
			shat3(abs(phi) > tau) = sign(phi(abs(phi) > tau))*tau;
        
			vhat3 = (shat3-s11(end,:))/z+v1(end,:);
        	vhat4 = (-shat3+s21(1,:))/z+v2(1,:);
			
            v1(end,:) = vhat3;
            s11(end,:) = shat3;
            v2(1,:) = vhat4;
            s21(1,:) = shat3;
        end

    end

end

function dx = calc_diff(x,dir)

% calculates finite difference approximation to x along direction dir
% note: does not divide by spatial grid spacing

	assert(dir == 1 | dir == 2);

	if dir == 1
		dx = (x([2:end end],:)-x([1 1:end-1],:))/2.;
    	dx(1,:) = dx(1,:)*2;
    	dx(end,:) = dx(end,:)*2;
	else
		dx = (x(:,[2:end end])-x(:,[1 1:end-1]))/2.;
    	dx(:,1) = dx(:,1)*2;
    	dx(:,end) = dx(:,end)*2;
	end
end

function dx = calc_diss(x,dir)

% calculates finite difference approximation to dissipation operator applied to x along direction dir
% note: does not divide by spatial grid spacing

	assert(dir == 1 | dir == 2);

	if dir == 1
		dx = (x([2:end end],:)-2.*x+x([1 1:end-1],:));
    	dx(1,:) = x(2,:)-x(1,:);
    	dx(end,:) = x(end,:)-x(end-1,:);
	else
		dx = (x(:,[2:end end])-2.*x+x(:,[1 1:end-1]));
    	dx(:,1) = x(:,2)-x(:,1);
    	dx(:,end) = x(:,end)-x(:,end-1);
	end

end