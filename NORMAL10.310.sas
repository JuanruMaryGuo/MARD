ods rtf close;
ods output close;

%let SampleSize = 5000;
%let mean1 = -0.8607960334;
%let mean2 = -0.2641159284;
%let sig1 = 0.2146751346;
%let sig2 = 0.2946509781;
%let rho = 0.3;
%let number_of_study = 10;
%let number_of_people = 100; 
%let file_name = 'C:\Users\scarlettwxy\Desktop\SIMULATE_NORMAL\5normal10.310.xlsx';



*Generating (p1,p2) pairs;
ods html;
proc iml;

Mean = {&mean1, &mean2};
var1 = &sig1*&sig1;
var2 = &sig2*&sig2;
var = &rho*&sig1*&sig2;
call symputx("var1",var1);
call symputx("var2",var2);
call symputx("var",var);
Cov = {&var1 &var, &var &var2};

call randseed(0);               
list = RandNormal(&SampleSize*&number_of_study, Mean, Cov);
list = cdf('normal',list);

do ij = 1 to &SampleSize;
	do s = 10 to 10*&number_of_study by 10;
		call symputx(cat(cat('p1',s),ij),list[(ij-1)*&number_of_study+s/10,1]);
		call symputx(cat(cat('p2',s),ij),list[(ij-1)*&number_of_study+s/10,2]);
	end;
end;
print(list);

quit;


ods html close;

%macro pcloglog;
 P1 = 1
%do i=-50 %to 50;
 - exp(sig1)/10*(exp(-exp(beta1 + exp(sig1)*&i/10))*pdf('normal', exp(sig1)*&i/10, 0, exp(sig1)))
%end;
;
 P0 = 1
%do i=-50 %to 50;
 - exp(sig0)/10*(exp(-exp(beta0 + exp(sig1)*&i/10))*pdf('normal', exp(sig0)*&i/10, 0, exp(sig0)))
%end;
;
%mend pcloglog;

*Estimation;
%Macro estimate(i);
*Generating dataset according to (p1,p2);
proc iml;
call randseed(&i);
y1 = j(&number_of_study,1);
n1 = j(&number_of_study,1) * &number_of_people;
y0 = j(&number_of_study,1);
n0 = j(&number_of_study,1) * &number_of_people;
id = 1:&number_of_study;

call randgen(y11, "Binomial", &pp11, &number_of_people);call randgen(y01, "Binomial", &pp21, &number_of_people);
call randgen(y12, "Binomial", &pp12, &number_of_people);call randgen(y02, "Binomial", &pp22, &number_of_people);
call randgen(y13, "Binomial", &pp13, &number_of_people);call randgen(y03, "Binomial", &pp23, &number_of_people);
call randgen(y14, "Binomial", &pp14, &number_of_people);call randgen(y04, "Binomial", &pp24, &number_of_people);
call randgen(y15, "Binomial", &pp15, &number_of_people);call randgen(y05, "Binomial", &pp25, &number_of_people);
call randgen(y16, "Binomial", &pp16, &number_of_people);call randgen(y06, "Binomial", &pp26, &number_of_people);
call randgen(y17, "Binomial", &pp17, &number_of_people);call randgen(y07, "Binomial", &pp27, &number_of_people);
call randgen(y18, "Binomial", &pp18, &number_of_people);call randgen(y08, "Binomial", &pp28, &number_of_people);
call randgen(y19, "Binomial", &pp19, &number_of_people);call randgen(y09, "Binomial", &pp29, &number_of_people);
call randgen(y110, "Binomial", &pp110, &number_of_people);call randgen(y010, "Binomial", &pp210, &number_of_people);

y1[1] = y11;y1[2] = y12;y1[3] = y13;y1[4] = y14;y1[5] = y15;y1[6] = y16;y1[7] = y17;y1[8] = y18;y1[9] = y19;y1[10] = y110;
y0[1] = y01;y0[2] = y02;y0[3] = y03;y0[4] = y04;y0[5] = y05;y0[6] = y06;y0[7] = y07;y0[8] = y08;y0[9] = y09;y0[10] = y010;
create rr var {y1 n1 y0 n0 id};
append; 
close rr;
quit;


*Adjust dataset;
data meta1;
 set rr; y=y1; n=n1; x=1; output; y=y0; n=n0; x=0; output; keep y n x id; run;
data meta2; set rr; x1=y1; x2=y0; n2=n0; y=0; run;


*Start estimation with 7 methods;
%put &syserr;
ods output AdditionalEstimates=RD1;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5 ;
 parms fz=3.5330, beta1=-2.2791, beta0=-1.7573, sig1=-4.1988, sig0=-4.3391;
 if x=1 then beta = beta1 + mu1 ;
 if x=0 then beta = beta0 + mu0 ;
 pred=probnorm(beta);
 P1 = probnorm(beta1/sqrt(1+exp(2*sig1)));
 P0 = probnorm(beta0/sqrt(1+exp(2*sig0)));
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 model y~binomial(N, pred);
 random mu1 mu0 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig0), exp(2*sig0)]) subject=id;
 estimate "RD" P1-P0;
run;
ods output;

%put &syserr;
ods output AdditionalEstimates=RD2;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5;
 parms beta1=-2.2791, beta0=-1.7573, sig=-4.2987;
 if x=1 then beta = beta1 + mu ;
 if x=0 then beta = beta0 + mu ;
 pred=probnorm(beta);
 P1 = probnorm(beta1/sqrt(1+exp(2*sig)));
 P0 = probnorm(beta0/sqrt(1+exp(2*sig)));
 model y~binomial(N, pred);
 random mu ~normal(0, exp(2*sig)) subject=id;
 estimate "RD" P1-P0;
run; 
ods output;

%put &syserr;
ods output AdditionalEstimates=RD3;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5 ;
 parms fz1=0.58368, beta11=-0.44688, beta00=-0.31925, sig11=-0.36694, sig00=-0.46919;
 fz = 10*fz1;
 beta1 = 10*beta11;
 beta0 = 10*beta00;
 sig1 = 10*sig11;
 sig0 = 10*sig00;
 if x=1 then beta = beta1 + mu1 ;
 if x=0 then beta = beta0 + mu0 ;
 pred=1/(1+exp(-beta));
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 P1 = 1/(1+exp(-beta1/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig1))**(1/2)));
 P0 = 1/(1+exp(-beta0/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig0))**(1/2)));
 model y~binomial(N, pred);
 random mu1 mu0 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig0), exp(2*sig0)]) subject=id;
 estimate "RD" P1-P0;
run;
ods output;

%put &syserr;
ods output AdditionalEstimates=RD4;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5;
 parms beta1=-4.4688, beta0=-3.1930, sig=-3.4609;
 if x=1 then beta = beta1 + mu ;
 if x=0 then beta = beta0 + mu ;
 sig1=sig; sig0=sig;
 pred=1/(1+exp(-beta));
 P1 = 1/(1+exp(-beta1/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig))**(1/2)));
 P0 = 1/(1+exp(-beta0/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig))**(1/2)));
 model y~binomial(N, pred);
 random mu ~normal(0, exp(2*sig)) subject=id;
 estimate "RD" P1-P0;
 run;
ods output;
/*
%put &syserr;
ods output AdditionalEstimates=RD5;
proc nlmixed data=meta1 fd df=1000 gtol=1e-10;
 parms fz=&fz, beta1=&mean1, beta0=&mean2, sig1=1, sig0=0;
 if x=1 then beta = beta1 + mu1 ;
 if x=0 then beta = beta0 + mu0 ;
 pred=1-exp(-exp((beta)));
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 model y~binomial(N, pred);
 random mu1 mu0 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig0), exp(2*sig0)]) subject=id;
 %pcloglog;
 estimate "RD" P1-P0;
run;
ods output;

%put &syserr;
ods output AdditionalEstimates=RD6;
proc nlmixed data=meta1 fd df=1000 gtol=1e-10;
 parms beta1=&mean1, beta0=&mean2, sig=0;
 if x=1 then beta = beta1 + mu ;
 if x=0 then beta = beta0 + mu ;
 pred=1-exp(-exp((beta)));
 sig1=sig; sig0=sig;
 model y~binomial(N, pred);
 random mu ~normal(0, exp(2*sig)) subject=id;
 %pcloglog;
 estimate "RD" P1-P0;
run;
ods output;
*/
%put &syserr;
ods output AdditionalEstimates=RD7;
proc nlmixed data=meta2 fd=200 df=1000 gtol=1e-5;
 parms a11=0.013040, a22=0.0114467, b11=0.0994078, b22=0.031273, ltw1=0.70496;
 bounds a11>0, a22>0, b11>0, b22>0;
 a1 = a11*10000;
 a2 = a22*100000;
 b1 = b11*100000;
 b2 = b22*1000000;
 ltw = ltw1*100;
 * tw = probnorm(ltw);
 tw = exp(ltw)/(1+exp(ltw));
 lb = (-(a1+b1)*(a2+b2)/max(a1*a2, b1*b2) ) ;
 rb = ( (a1+b1)*(a2+b2)/max(a1*b2, a2*b1) );
 w = tw*(rb-lb)+lb;
 P1=a1/(a1+b1);
 P0=a2/(a2+b2);
 logL1=lgamma(n1+1)-lgamma(x1+1)-lgamma(n1-x1+1)+lgamma(a1+b1)-lgamma(a1)-lgamma(b1)+lgamma(a1+x1)+lgamma(n1+b1-X1)-lgamma(a1+b1+n1);
 logL2=lgamma(n2+1)-lgamma(x2+1)-lgamma(n2-x2+1)+lgamma(a2+b2)-lgamma(a2)-lgamma(b2)+lgamma(a2+x2)+lgamma(n2+b2-X2)-lgamma(a2+b2+n2);
 logL3=log(1+w*(x1-n1*a1/(a1+b1))*(x2-n2*a2/(a2+b2))/(a1+b1+n1)/(a2+b2+n2));
 logL = (logL1+logL2+LogL3); 
 model y~general(logL);
 estimate "RD" P1-P0;
run; 
ods output;

*Save all the results in the Result dataset;
proc sql;
create table result as
select Rd1.Estimate as rd1, Rd2.Estimate as rd2, Rd3.Estimate as rd3, Rd4.Estimate as rd4, Rd7.Estimate as rd7
from Rd1,Rd2,Rd3,Rd4,Rd7
union all
select * from result;
quit;

*Save all the results in the Result dataset;
proc sql;
create table resultlower as
select Rd1.Lower as rd1, Rd2.Lower as rd2, Rd3.Lower as rd3, Rd4.Lower as rd4, Rd7.Lower as rd7
from Rd1,Rd2,Rd3,Rd4,Rd7
union all
select * from resultlower;
quit;

*Save all the results in the Result dataset;
proc sql;
create table resultupper as
select Rd1.Upper as rd1, Rd2.Upper as rd2, Rd3.Upper as rd3, Rd4.Upper as rd4, Rd7.Upper as rd7
from Rd1,Rd2,Rd3,Rd4,Rd7
union all
select * from resultupper;
quit;
%Mend estimate;

*Initialize the result dataset;
proc sql;
create table result
(
rd1 numeric(1,8),
rd2 numeric(1,8),
rd3 numeric(1,8),
rd4 numeric(1,8),
rd7 numeric(1,8)
);
quit;

proc sql;
create table resultlower
(
rd1 numeric(1,8),
rd2 numeric(1,8),
rd3 numeric(1,8),
rd4 numeric(1,8),
rd7 numeric(1,8)
);
quit;

proc sql;
create table resultupper
(
rd1 numeric(1,8),
rd2 numeric(1,8),
rd3 numeric(1,8),
rd4 numeric(1,8),
rd7 numeric(1,8)
);
quit;


%Macro simulate;
%do i=1 %to &SampleSize;
	%put &i;
	%let pp11 = &&p110&i; %let pp21 = &&p210&i; 
	%let pp12 = &&p120&i; %let pp22 = &&p220&i; 
	%let pp13 = &&p130&i; %let pp23 = &&p230&i;
	%let pp14 = &&p140&i; %let pp24 = &&p240&i; 
	%let pp15 = &&p150&i; %let pp25 = &&p250&i; 
	%let pp16 = &&p160&i; %let pp26 = &&p260&i;
	%let pp17 = &&p170&i; %let pp27 = &&p270&i; 
	%let pp18 = &&p180&i; %let pp28 = &&p280&i; 
	%let pp19 = &&p190&i; %let pp29 = &&p290&i;
	%let pp110 = &&p1100&i; %let pp210 = &&p2100&i; 
	%estimate(&i);
%end;
%Mend simulate;

%simulate;

*Save dataset;
proc export data=result outfile=&file_name dbms=xlsx replace; sheet = "result";
run;
proc export data=resultlower outfile=&file_name dbms=xlsx; sheet = "resultlower";
run;
proc export data=resultupper outfile=&file_name dbms=xlsx; sheet = "resultupper";
run;

ods html;
proc means data = result;
run;

proc means data = resultupper missing;
run;
ods html close;
