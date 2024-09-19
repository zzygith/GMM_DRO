* Define sets
Sets
    UH /UH1*UH3/
    LH /LH1*LH3/
    i  /1*10/;

* Define parameters
*Scalar distanceThreshold /1.0/;
Scalar distanceThreshold /10.0/;

Parameter
    wUH(UH) /UH1 0.126, UH2 0.622, UH3 0.252/
    wLH(LH) /LH1 0.126, LH2 0.622, LH3 0.252/;

Table miuUH(UH,i)
         1      2       3       4       5       6       7       8       9       10
    UH1  4.41   8.75    12.98   17.54   20.77   25.52   29.71   34.27   38.40   42.94
    UH2  0.88   1.43    2.28    3.11    3.85    4.42    5.37    5.90    6.97    7.25
    UH3  2.60   5.16    7.69    9.86    12.91   15.14   17.49   20.37   22.33   24.99;


* Define parameter index(i)
Parameter index(i);
index(i) = ord(i);

Table alphaUH(UH,i)
         1      2       3       4       5       6       7       8       9       10
    UH1  1.94   3.31    6.67    6.59    8.26    9.97    11.40   12.94   13.73   20.22
    UH2  1.92   3.63    5.27    7.26    8.79    12.12   11.66   14.56   15.50   15.69
    UH3  1.93   3.88    6.43    7.70    7.65    10.84   13.48   14.14   20.26   21.38;


* Define variables
Variables
    miuLH(LH,i)
    alphaLH(LH,i)
    pi(UH,LH)
    z;

* Set variable bounds
miuLH.lo(LH,i)    = -10000;
miuLH.up(LH,i)    =  10000;

alphaLH.lo(LH,i)  = 0;
alphaLH.up(LH,i)  = 10000;

pi.lo(UH,LH)      = 0;
pi.up(UH,LH)      = 1;

* Initial values for variables
Parameter pi_init(UH,LH) /
    UH1.LH1 0.126,  UH1.LH2 0.0,  UH1.LH3 0.0,
    UH2.LH1 0.0,  UH2.LH2 0.622, UH2.LH3 0.0,
    UH3.LH1 0.0,  UH3.LH2 0.0,  UH3.LH3 0.252 /;

pi.L(UH,LH) = pi_init(UH,LH);

* Establish mapping between UH and LH
Set mapping(UH,LH) /
    UH1.LH1,
    UH2.LH2,
    UH3.LH3 /;

* Initialize miuLH and alphaLH with corresponding miuUH and alphaUH values
Loop((UH,LH)$(mapping(UH,LH)),
    miuLH.L(LH,i)    = miuUH(UH,i);
    alphaLH.L(LH,i)  = alphaUH(UH,i);
);

* Define equations
Equations
    ec1(UH)
    ec2(LH)
    distanceC
    obj1;

* Constraints
ec1(UH).. sum(LH, pi(UH,LH)) =e= wUH(UH);

ec2(LH).. sum(UH, pi(UH,LH)) =e= wLH(LH);

* Distance constraint
distanceC..
    sum((UH,LH),
        (
            sum(i, sqr(miuUH(UH,i) - miuLH(LH,i)))
          + sum(i, alphaUH(UH,i) + alphaLH(LH,i))
          - 2 * sum(i, sqrt(alphaUH(UH,i) * alphaLH(LH,i)))
        ) * pi(UH,LH)
    ) =l= distanceThreshold;

* Objective function
obj1..
*    z =e=
*       sum(i, sqr(miuUH('UH1',i) - miuLH('LH1',i)))
*      + sum(i, alphaUH('UH1',i) + alphaLH('LH1',i))
*      - 2 * sum(i, sqrt(alphaUH('UH1',i) * alphaLH('LH1',i)));
*   z =e=
*       sum(i, sqr(miuUH('UH1',i) - miuLH('LH2',i)))
*     + sum(i, alphaUH('UH1',i) + alphaLH('LH2',i))
*     - 2 * sum(i, sqrt(alphaUH('UH1',i) * alphaLH('LH2',i)));
*    z =e=
*        sum(i, sqr(miuUH('UH1',i) - miuLH('LH3',i)))
*      + sum(i, alphaUH('UH1',i) + alphaLH('LH3',i))
*      - 2 * sum(i, sqrt(alphaUH('UH1',i) * alphaLH('LH3',i)));
*    z =e=
*        sum(i, sqr(miuUH('UH2',i) - miuLH('LH1',i)))
*      + sum(i, alphaUH('UH2',i) + alphaLH('LH1',i))
*      - 2 * sum(i, sqrt(alphaUH('UH2',i) * alphaLH('LH1',i)));
*   z =e=
*       sum(i, sqr(miuUH('UH2',i) - miuLH('LH2',i)))
*     + sum(i, alphaUH('UH2',i) + alphaLH('LH2',i))
*     - 2 * sum(i, sqrt(alphaUH('UH2',i) * alphaLH('LH2',i)));
*    z =e=
*        sum(i, sqr(miuUH('UH2',i) - miuLH('LH3',i)))
*      + sum(i, alphaUH('UH2',i) + alphaLH('LH3',i))
*      - 2 * sum(i, sqrt(alphaUH('UH2',i) * alphaLH('LH3',i)));
*    z =e=
*        sum(i, sqr(miuUH('UH3',i) - miuLH('LH1',i)))
*      + sum(i, alphaUH('UH3',i) + alphaLH('LH1',i))
*      - 2 * sum(i, sqrt(alphaUH('UH3',i) * alphaLH('LH1',i)));
   z =e=
       sum(i, sqr(miuUH('UH3',i) - miuLH('LH2',i)))
     + sum(i, alphaUH('UH3',i) + alphaLH('LH2',i))
     - 2 * sum(i, sqrt(alphaUH('UH3',i) * alphaLH('LH2',i)));
*    z =e=
*        sum(i, sqr(miuUH('UH3',i) - miuLH('LH3',i)))
*      + sum(i, alphaUH('UH3',i) + alphaLH('LH3',i))
*      - 2 * sum(i, sqrt(alphaUH('UH3',i) * alphaLH('LH3',i)));

* Define model
Model mymodel /all/;

* Select solver
*option nlp=ANTIGONE;
Option nlp = baron;

* Solve the model
*Solve mymodel maximizing z using nlp;
Solve mymodel minimizing z using nlp;

* Display results
Display pi.l;
