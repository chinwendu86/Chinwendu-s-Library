$title Stochastic Two-stage program

$ontext

Coordinating Pre- and Post-Disaster Investments at Multiple Locations

$offtext

Sets  i 'potential affected locations' / PAL1, PAL2 /
      j 'main distribution centers' / MDC1, MDC2 /
      f ' Budget sensitivity'       /f1*f12/

;

Parameters

alpha(i) 'unit cost of pre-disaster preparedness at PAL' /PAL1 0.1,  PAL2 0.1/

beta(j)  'unit cost of pre-disaster preparedness at MDC' /MDC1 0.05, MDC2 0.1/

gam(i) 'unit cost of relief at PAL'                  /PAL1 0.1, PAL2 0.1/

lamda(i)  'Coeeffcient of the damage function at PAL'   /PAL1 10, PAL2 10/

v(i)      'Coeefficient of the damage function at PAL' /PAL1 10, PAL2 10/

Ecl(i)       'Local preparedness effectiveness at potential affected location'     /PAL1 0.1, PAL2 0.1/

Ecd(j)       'Preparedness effectiveness at main distribution center'              /MDC1 0.1, MDC2 0.1/

Er(i)        'Local relief at potential affected location'                          /PAL1 0.1, PAL2 0.1/

Ecr(i)       'Joint effectiveness of local preparedness and relief at PAL'          /PAL1 0.1, PAL2 0.1/

B_loop(f)   'Sensitivity for Budget'                                               /f1 0.1, f2 0.2, f3 0.3, f4 0.4, f5 0.5, f6 0.6, f7 0.7, f8 0.8, f9 0.9, f10 1.0, f11 1.2, f12 1.4/
;

table Ecdr(i,j) 'Joint effectiveness of the preparedness at MDC and relief at PAL'


         MDC1     MDC2

   PAL1  0.1        0.1

   PAL2  0.1        0.1
;


scalars

P12 'Conditional Probability of disaster magnitude at location 1 given that disaster has occurred at location 2' /0.01/

P1 'Probability of disaster magnitude at location 1' / 0.01/

P2 'Probability of disaster magnitude at location 2'  / 0.01/

m1m 'Disaster magnitude at location 1' / 1 /

m2m 'Disaster magnitude at location 2' / 1/

budget /1/
;

positive variables
                   cd(j) 'pre-disaster preparedness effort main distribution center MDC',
                   r(i)  'post-disaster relief effort at PAL',
                   cl(i) 'pre-disaster preparedness effort at PAL';

binary variable y(i,j)'whether PAL location i seeks post disaster relief from MDC j';

variable z    'loss minimized';

Equation

obj
budgetlimit

;

obj.. z   =e= sum(i,Alpha(i)*cl(i)) + sum(j,Beta(j)*cd(j)) + sum(i,Gam(i)*r(i))+ prod(i,v(i)*(1-exp(-lamda(i)*m1m))*exp(-Ecl(i)*cl(i)-sum(j,Ecd(j)*cd(j)* y(i,j))-Er(i)*r(i)-Ecr(i)*cl(i)*r(i)-sum(j,Ecdr(i,j)*cd(j)*r(i)*y(i,j))))*(P12*P2)

               +sum(i,Gam(i)*r(i))+ prod(i,v(i)*(1-exp(-lamda(i)*m2m))*exp(-Ecl(i)*cl(i)-sum(j,Ecd(j)*cd(j)*y(i,j))-Er(i)*r(i)-Ecr(i)*cl(i)*r(i)-sum(j,Ecdr(i,j)*cd(j)*r(i)* y(i,j))))*(P12*P2)

               + (P1-(P12*P2))*sum(i,Gam(i)*r(i)) + sum(i,Gam(i)*r(i))+ prod(i,v(i)*(1-exp(-lamda(i)*m1m))*exp(-Ecl(i)*cl(i)-sum(j,Ecd(j)*cd(j)* y(i,j))-Er(i)*r(i)-Ecr(i)*cl(i)*r(i)-sum(j,Ecdr(i,j)*cd(j)*r(i)*y(i,j))))*(P1*(1-P12*P2/P1))

               + ((1-P12)*P2)*sum(i,Gam(i)*r(i))+sum(i,Gam(i)*r(i))+ prod(i,v(i)*(1-exp(-lamda(i)*m2m))*exp(-Ecl(i)*cl(i)-sum(j,Ecd(j)*cd(j)* y(i,j))-Er(i)*r(i)-Ecr(i)*cl(i)*r(i)-sum(j,Ecdr(i,j)*cd(j)*r(i)*y(i,j))))*((1-P12)*P2)

               + sum(i,Gam(i)*r(i))*(1-P1-P2+P12*P2) + sum(i,Gam(i)*r(i))*(1-P1-P2+P12*P2);


budgetlimit.. budget=l= sum(i,alpha(i)*cl(i))+ sum(j,beta(j)*cd(j)) + sum(i,gam(i)*cl(i));


Model StochasticTwostage /All/ ;

solve StochasticTwostage minimizing z usin minlp;

*y.l(i,j) $(cd.l(j)=0) = 0

display z.l, cl.l, cd.l, y.l, r.l, p1, p2, budget,m1m, m2m;

$ontext
loop(f,
       budget = B_loop(f);
       solve StochasticTwostage minimizing z usin minlp;
       display budget, P1, P2, P12, m1m, m2m
       display z.l, cd.l, r.l, cl.l, y.l;
);
$offtext
