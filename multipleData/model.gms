Parameter a;
Set i /1*2/;
parameter a(i);
$include a_2.txt

Variable x;
Equation constraint;
constraint.. x =g= a('2');

Model lpModel /constraint/;
Solve lpModel using LP minimizing x;