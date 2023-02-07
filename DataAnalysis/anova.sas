
/****************************************************************************************/
/*       Exploratory analysis to compare the inter-group and intra-group  correlation   */
/****************************************************************************************/

PROC IMPORT OUT= WORK.aov 
            DATAFILE= "Source/aov.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data aov_2;
set aov; 
delta_CT_num = input(deltaCT, comma10.);
run;

/****************************************************************************************/
/*         Running proc mixed to fit the random effect model in the supplement          */
/****************************************************************************************/

proc mixed data=aov_2;
class group subject;
model delta_Ct_num=;
random subject group(subject);
run;
quit;

