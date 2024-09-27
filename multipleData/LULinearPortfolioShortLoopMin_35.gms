* Define sets
Sets
    UH /UH1*UH3/
    LH /LH1*LH3/
    i  /1*10/;

* Define parameters
Scalar distanceThreshold /35.0/;

Parameters wUH(UH);
$include wUH.txt

Parameters wLH(LH);
$include wLH.txt

Parameter miuUH(UH,i);
$include miuUH.txt

* Define parameter index(i)
Parameter index(i);
index(i) = ord(i);

Parameter alphaUH(UH,i);
$include alphaUH.txt

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
    UH1.LH1 0.0,  UH1.LH2 0.0,  UH1.LH3 0.0,
    UH2.LH1 0.0,  UH2.LH2 0.0,  UH2.LH3 0.0,
    UH3.LH1 0.0,  UH3.LH2 0.0,  UH3.LH3 0.0 /;

pi.l(UH,LH) = pi_init(UH,LH);

* Establish mapping relationship between UH and LH
Set mapping(UH,LH) /
    UH1.LH1,
    UH2.LH2,
    UH3.LH3 /;

* Initialize miuLH and alphaLH with corresponding miuUH and alphaUH
Loop((UH,LH)$(mapping(UH,LH)),
    miuLH.l(LH,i)    = miuUH(UH,i);
    alphaLH.l(LH,i)  = alphaUH(UH,i);
);

* Define aliases
Alias (UH, UHval), (LH, LHval);

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

* Define parameter
Parameter isActive(UH,LH);

* Define objective function
obj1..
    z =e=
        sum((UH,LH)$(isActive(UH,LH)),
            sum(i, sqr(miuUH(UH,i) - miuLH(LH,i)))
          + sum(i, alphaUH(UH,i) + alphaLH(LH,i))
          - 2 * sum(i, sqrt(alphaUH(UH,i) * alphaLH(LH,i)))
        );

* Define model
Model mymodel / ec1, ec2, distanceC, obj1 /;

* Relative convergence tolerance
mymodel.OptCR = 0.01;
* Absolute convergence tolerance
mymodel.OptCA = 0.1;
* Set maximum solution time (e.g., 10 seconds)
mymodel.ResLim = 10;

* Choose solver
Option nlp = baron;

* Define objective function combination set
Sets
    ObjComb(UHval,LHval);
ObjComb(UHval,LHval) = YES;
* Iterate over all UH and LH combinations

* Define result parameters
Parameter
    z_results(UH,LH)               "Objective function value"
    pi_results(UH,LH,UHval,LHval)  "Value of pi variable";

* Loop through objective function combinations
Loop(ObjComb(UHval,LHval),
* Reset isActive parameter
    isActive(UH,LH) = 0;

* Set the current combination as active
    isActive(UHval,LHval) = 1;

* Reset variable initial values (if necessary)
    Loop(i,
        miuLH.l(LHval,i)    = miuUH(UHval,i);
        alphaLH.l(LHval,i)  = alphaUH(UHval,i);
    );

* Solve the model
    Solve mymodel minimizing z using nlp;
*    Solve mymodel maximizing z using nlp;

* Save results
    z_results(UHval,LHval) = z.l;
    pi_results(UH,LH,UHval,LHval) = pi.l(UH,LH);

    Display 'Objective Value =', z.l;
);

* Display all results
Display 'Final results for saving =', z_results;

file resultsFile /'resultsMIN_35.txt'/;
put resultsFile;
loop(UHval,
	loop(LHval,
		put "z_results(", UHval.tl, ",", LHval.tl, ") = " z_results(UHval,LHval) /;
	);
);
putclose resultsFile;
Display 'Results saved to results.txt';