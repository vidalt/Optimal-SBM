param n;
param m;
param K;
param Cnst;
param w_upper;

set V := {1..n};  #set of nodes
set C := {1..K};  #set of clusters

param A {V,V} >= 0; #adjacency matrix
param k {V} >= 0; #vertex degrees


var z {V,C} binary;  #community assignments
var w {C,C} >=0;     #affinity matrix

minimize L:  - Cnst + 0.5*sum {i in V, j in V, r in C, s in C} (if (A[i,j] > 0) then ( (-A[i,j]*log(w[r,s]) + ((k[i]*k[j])/(2*m))*w[r,s] ) * z[i,r] * z[j,s] ) else ( ((k[i]*k[j])/(2*m))*w[r,s] * z[i,r] * z[j,s]) );

s.t. Constraint {i in V}:
	sum {r in C} (z[i,r]) == 1;

s.t. Wupper {r in C, s in C}:
    w[r,s] <= w_upper;

s.t. Symm_Break_1:
	z[1,1] == 1;

s.t. Symm_Break {r in {2..(K-1)}, j in {r..n}}:
	sum {l in 1..(r-1), i in 2..(j-1)} (z[i,l]) - sum {l in {1..r}} (z[j,l]) <= j - 3;

