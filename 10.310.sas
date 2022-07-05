ods rtf close;
ods output close;

%let SampleSize = 5000;
%let alpha1 = 6.863515;
%let alpha2 = 5.597307;
%let beta1 = 27.45406;
%let beta2 = 8.395961;
%let rho = 0.3;
%let number_of_study = 10;
%let number_of_people = 100; 
%let file_name = 'C:\Users\scarlettwxy\Desktop\5beta10.310';



*Generating (p1,p2) pairs;
ods html;
proc iml;
start alpha_function(m,v);
	return( m*(((m*(1-m))/v)-1) );
finish;

start beta_function(m,v);
	return( ((m*(1-m))/v-1)*(1-m) );
finish;

start g_part_function(p1,p2,w,alpha1,alpha2,beta1,beta2);
	mu1 = alpha1/(alpha1+beta1);
	mu2 = alpha2/(alpha2+beta2);
	return((1+w*(p1-mu1)#(p2-mu2)));
finish;

start w_function(rho,alpha1,alpha2,beta1,beta2);
  var1 = alpha1*beta1/(((alpha1+beta1)**2)*(alpha1+beta1+1));
  var2 = alpha2*beta2/(((alpha2+beta2)**2)*(alpha2+beta2+1));
  return(rho/sqrt(var1*var2));
finish;

start M_function(w);
	return(1+w);
finish;

start delete_zeros(nn);
	cols = loc(nn>0);
	D = nn[cols,];
	return(D);
finish;

start bibeta_sampling(n,alpha1,alpha2,beta1,beta2,rho);
	w = w_function(rho,alpha1,alpha2,beta1,beta2);
	M = M_function(w);
	num = 0;
	sampling1 = {};
	sampling2 = {};
	do while(nrow(sampling1)<n);
		p1 = j(n,1);
		p2 = j(n,1);
		call randgen(p1, "Beta", alpha1, beta1);
		call randgen(p2, "Beta", alpha2, beta2);
		g_part = g_part_function(p1,p2,w,alpha1,alpha2,beta1,beta2);
		flag = g_part/M;
		u = j(n,1);
		call randgen(u, "Uniform");
		keep = (u<=flag);
		p1 = p1 # keep;
		p2 = p2 # keep;
		if max(p1) >0 then 
			do;
				p1 = delete_zeros(p1);
				p2 = delete_zeros(p2);
				sampling1 = sampling1//p1;
				sampling2 = sampling2//p2;
		end;
	end;
	return(sampling1||sampling2);
finish;


call randseed(0);
list = bibeta_sampling(&SampleSize*&number_of_study,&alpha1,&alpha2,&beta1,&beta2,&rho);
do ij = 1 to &SampleSize;
	do s = 10 to 10*&number_of_study by 10;
		call symputx(cat(cat('p1',s),ij),list[(ij-1)*&number_of_study+s/10,1]);
		call symputx(cat(cat('p2',s),ij),list[(ij-1)*&number_of_study+s/10,2]);
	end;
end;
print(list);
quit;
ods html close;
/*
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
*/

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
 parms fz1=0.005393, beta11=-0.08204, beta00=-0.02070, sig11=-0.2538, sig00=-0.0878;
 fz = 100*fz1;
 beta1 = 10*beta11;
 beta0 = 10*beta00;
 sig1 = 10*sig11;
 sig0 = 10*sig00;
 if x=1 then beta = beta1 + mu1 ;
 if x=0 then beta = beta0 + mu0 ;
 pred=probnorm(beta);
 P1 = probnorm(beta1/sqrt(1+exp(2*sig1)));
 P0 = probnorm(beta0/sqrt(1+exp(2*sig0)));
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 model y~binomial(N, pred);
 random mu1 mu0 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig0), exp(2*sig0)]) subject=id;
 estimate "RD" P1-P0;
 /*estimate "fz" fz1;
 estimate "beta1" beta11;
 estimate "beta0" beta00;
 estimate "sig1" sig11;
 estimate "sig0" sig00;*/
run;
ods output;

%put &syserr;
ods output AdditionalEstimates=RD2;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5;
 parms beta1=-0.8448, beta0=-0.2605, sig=-10.8173;
 if x=1 then beta = beta1 + mu ;
 if x=0 then beta = beta0 + mu ;
 pred=probnorm(beta);
 P1 = probnorm(beta1/sqrt(1+exp(2*sig)));
 P0 = probnorm(beta0/sqrt(1+exp(2*sig)));
 model y~binomial(N, pred);
 random mu ~normal(0, exp(2*sig)) subject=id;
 estimate "RD" P1-P0;
 /*estimate "beta1" beta1;
 estimate "beta0" beta0;
 estimate "sig" sig;*/
run; 
ods output;

%put &syserr;
ods output AdditionalEstimates=RD3;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5 ;
 parms fz1=0.007477, beta11=-0.01699, beta00=-0.00636, sig11=-0.0099, sig00=-0.00866;
 fz = 100*fz1;
 beta1 = 100*beta11;
 beta0 = 100*beta00;
 sig1 = 100*sig11;
 sig0 = 100*sig00;
 if x=1 then beta = beta1 + mu1 ;
 if x=0 then beta = beta0 + mu0 ;
 pred=1/(1+exp(-beta));
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 P1 = 1/(1+exp(-beta1/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig1))**(1/2)));
 P0 = 1/(1+exp(-beta0/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig0))**(1/2)));
 model y~binomial(N, pred);
 random mu1 mu0 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig0), exp(2*sig0)]) subject=id;
 estimate "RD" P1-P0;
 /*estimate "fz" fz1;
 estimate "beta1" beta11;
 estimate "beta0" beta00;
 estimate "sig1" sig11;
 estimate "sig0" sig00;*/
run;
ods output;

%put &syserr;
ods output AdditionalEstimates=RD4;
proc nlmixed data=meta1 fd df=1000 gtol=1e-5;
 parms beta1=-1.3918, beta0=-0.4169, sig=-9.8826;
 if x=1 then beta = beta1 + mu ;
 if x=0 then beta = beta0 + mu ;
 sig1=sig; sig0=sig;
 pred=1/(1+exp(-beta));
 P1 = 1/(1+exp(-beta1/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig))**(1/2)));
 P0 = 1/(1+exp(-beta0/(1+(16*sqrt(3)/(15*3.14159))**2*exp(2*sig))**(1/2)));
 model y~binomial(N, pred);
 random mu ~normal(0, exp(2*sig)) subject=id;
 estimate "RD" P1-P0;
 /*estimate "beta1" beta1;
 estimate "beta0" beta0;
 estimate "sig" sig;*/
 run;
ods output;

/*
%put &syserr;
ods output AdditionalEstimates=RD5;
proc nlmixed data=meta1 fd df=1000 gtol=1e-10;
 parms fz=&fz, beta1=&mean1, beta0=&mean2, sig1=&sig1, sig0=&sig0;
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
 parms beta1=&mean1, beta0=&mean2, sig=&sig;
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
 parms a11=0.006123, a22=0.000181, b11=0.002463, b22=0.000273, ltw=704.96;
 bounds a11>0, a22>0, b11>0, b22>0;
 * tw = probnorm(ltw);
 a1 = a11*10000;
 a2 = a22*100000;
 b1 = b11*100000;
 b2 = b22*100000;
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
 /*estimate "a11" a11;
 estimate "b11" b11;
 estimate "a22" a22;
 estimate "b22" b22;*/
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
proc means data = resultupper;
run;
ods html close;
