* Define sets
Sets
    UH /UH1*UH3/
    LH /LH1*LH3/
    i  /1*10/;

* Define parameters
Scalar distanceThreshold /1.0/;

Parameter
    wUH(UH) /UH1 0.6, UH2 0.25, UH3 0.15/
    wLH(LH) /LH1 0.6, LH2 0.25, LH3 0.15/;

Table miuUH(UH,i)
         1      2       3       4       5       6       7       8       9       10
    UH1  0.75   1.5     2.25    3.0     3.75    4.5     5.25    6.0     6.75    7.5
    UH2  2.5    5.0     7.5     10.0    12.5    15.0    17.5    20.0    22.5    25.0
    UH3  4.25   8.5     12.75   17.0    21.25   25.5    29.75   34.0    38.25   42.5;

* Define parameter index(i)
Parameter index(i);
index(i) = ord(i);

*Parameter alphaUH(UH,i);

* Calculate alphaUH(UH,i)
*alphaUH(UH,i) = index(i) * 1.8;
Table alphaUH(UH,i)
         1      2       3       4       5       6       7       8       9       10
    UH1  1.8    3.6     5.4     7.2     9.0     10.8    12.6    14.4    16.2    18.0
    UH2  1.8    3.6     5.4     7.2     9.0     10.8    12.6    14.4    16.2    18.0
    UH3  1.8    3.6     5.4     7.2     9.0     10.8    12.6    14.4    16.2    18.0;

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
    UH1.LH1 0.6,  UH1.LH2 0.0,  UH1.LH3 0.0,
    UH2.LH1 0.0,  UH2.LH2 0.25, UH2.LH3 0.0,
    UH3.LH1 0.0,  UH3.LH2 0.0,  UH3.LH3 0.15 /;

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
*        sum(i, sqr(miuUH('UH1',i) - miuLH('LH1',i)))
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
Solve mymodel maximizing z using nlp;
*Solve mymodel minimizing z using nlp;

* Display results
Display pi.l;
