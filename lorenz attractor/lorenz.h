#include <cmath>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cerrno>
#include <unistd.h> // for read
#include "gogo.h"

using namespace std;

typedef double     dbl;
typedef const dbl cdbl;

cdbl PGROW = -0.2;
cdbl PSHRINK = -0.25;
cdbl FCOR = 0.066666667;  // 1/15 for 5th order correction
cdbl SAFETY = 0.9;
cdbl ERRCON = 6.0e-4;     // about (4/SAFETY)^(1/PGROW)
cdbl EPS    = 2.0e-5;     // might change this
cdbl INSTEP = 0.005;
cdbl TINY   = 1.0e-30;


// prototypes
void derivs(dbl, vector<dbl>&,  vector<dbl>&);     // user-supplied equations
void rk(int, dbl, dbl, vector<dbl>&,  vector<dbl>&,
        void(*) (dbl, vector<dbl>&, vector<dbl>&));
void rkqc(int, dbl, dbl, vector<dbl>&, vector<dbl>&,
          dbl&, dbl&, void(*) (dbl, vector<dbl>&, vector<dbl>&));
void init_screen(flt, flt, flt, flt);   // sets up tek window
