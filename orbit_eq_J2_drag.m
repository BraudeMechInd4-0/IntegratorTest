function [dr] = orbit_eq_J2_drag(t,r,mu,CD,A,m,Re,J2)
%% function [dr] = orbit_eq_J2_drag(t,r,mu,CD,A,m,Re,J2)
% definition of orbit eq with J2 and drag


% d^2r/dt = -mu*r/norm(r);

[n,M] = size(r);
dr = zeros (size(r));


if n == 6

    normr = sqrt(sum(r(1:3,:).^2,1));

    ad = drag_accel(r',CD,A,m,Re)';
    aj2 = -3/2*J2*(mu./(normr.^2)).*(Re./normr).^2.*[(1-5*(r(3,:)./normr).^2).*r(1,:)./normr, (1-5*(r(3,:)./normr).^2).*r(2,:)./normr, (3-5*(r(3,:)./normr).^2).*r(3,:)./normr]';
    
    le =normr.^3;

    dr(1,:) = r(4,:);
    dr(2,:) = r(5,:);
    dr(3,:) = r(6,:);
    dr(4,:) = -mu*r(1,:)./le + ad(1,:)+aj2(1,:);
    dr(5,:) = -mu*r(2,:)./le + ad(2,:)+aj2(2,:);
    dr(6,:) = -mu*r(3,:)./le + ad(3,:)+aj2(3,:);
elseif M == 6

    normr = sqrt(sum(r(:,1:3).^2,2));

    ad = drag_accel(r,CD,A,m,Re);
    aj2 = -3/2*J2*(mu./(normr.^2)).*(Re./normr).^2.*[(1-5*(r(:,3)./normr).^2).*r(:,1)./normr, (1-5*(r(:,3)./normr).^2).*r(:,2)./normr, (3-5*(r(:,3)./normr).^2).*r(:,3)./normr];

    le =normr.^3;

    dr(:,1) = r(:,4);
    dr(:,2) = r(:,5);
    dr(:,3) = r(:,6);
    dr(:,4) = -mu*r(:,1)./le + ad(:,1)+aj2(:,1);
    dr(:,5) = -mu*r(:,2)./le + ad(:,2)+aj2(:,2);
    dr(:,6) = -mu*r(:,3)./le + ad(:,3)+aj2(:,3);
end






