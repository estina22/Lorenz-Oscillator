// for solving the lorenz equations.  creates a data
// file called "output" that has, on each succeeding
// line, the scaled values of x, y and z.  this file
// is then suitable for viewing with gnuplot (see
// the gnuplot script "lorenz.gp"). the output file
// may also be used as input for the "sonify"
// program to create "output.mup", which in turn may
// be processed with "mup" ("mup -M output.mup") to
// create the MIDI sound file "output.mid", which
// provides an audio representation of the solution
// set.  traditional parameter values are used.
// 
// the lorenz equations are solved using a 4th order
// runge-kutta adaptive step-size routine, patterned
// after a routine described in numerical recipes
// for C (1st edition).  see the files rkqc.c and rk.c
// for details.
//
// created 2/21/06.

#include "lorenz.h"
#include <iomanip>  // for starting point display format

vector<dbl> p;      // global parameter vector
int main(int argc, char** argv)
{
    int nvar = 3;   // number of variables
    dbl t = 0.0;    // time
    dbl h = INSTEP; // initial time step
    dbl hnext;
    dbl hdid;
    vector<dbl> u (nvar);
    vector<dbl> dudt (nvar);
    vector<int> gg (2);
    vector<flt> uu (3);
    p.push_back(10.0);    // p[0]: prandtl number sigma
    p.push_back(25.0);    // p[1]: rayleigh number r
    p.push_back(8.0/3.0); // p[2]: aspect ratio b

    int blocksize = 50;
    string ss;
    cerr << "blocksize [50]: "; // 50 is a good number
    if (getline(cin, ss)){
        if (ss == ""){
            blocksize = 50;
        }
        else {
            blocksize = atoi(ss.c_str());
        }
    }
    cout << "blocksize = " << blocksize << "\n";
    cerr << "\nhit \'x\' or \'q\' to exit\n";
    u[0] = 0.6;            // initial value of x
    u[1] = 0.65;            // initial value of y
    u[2] = 0.7;            // initial value of z

    derivs(t, u, dudt);    // find du/dt at initial point
    uns count = 0;         // to count "for" loop iterations
    char c = 0;            // to hold input char
    dbl delr = 0.4;        // amount by which r (p[1]) is changed
    dbl del_beta = M_PI/16.0;
    dbl beta = M_PI/4.0;   // viewing angle
    dbl xx = u[0]*cos(beta) + u[1]*sin(beta);
    int quit = false;
    init_screen(xx, u[2], p[1], beta);
    M(xx, u[2]);        // starting point
    hotkeys(on);  // sets "read()" wait time to 0.1 sec among other
                  // things.  see the hotkeys.c file for more info
    for (t=0.0 ; ; t += hdid, count++){
        if (!(count%blocksize)){
            read(0, &c, 1);
        }
        switch (c) {
            case 0:
                break;
            case 'x': case 'q':
                quit = true;
                break;
            case 'r': // decrease r by delr
                CLEAR();  // erase screen
                AXES(0.0, 0.0, 0.5, 0.5, 1.0, 1.0);
                p[1] -= delr;
                M(-3.0, 2.0);
                TEXT(0, "beta = ", beta *180.0/M_PI);
                TEXT(0, "r = ", p[1]);
                break;
            case 'R': // increase r by delr
                CLEAR();
                AXES(0.0, 0.0, 0.5, 0.5, 1.0, 1.0);
                p[1] += delr;
                M(-3.0, 2.0);
                TEXT(0, "beta = ", beta *180.0/M_PI);
                TEXT(0, "r = ", p[1]);
                break;
            case 'b': // decrease beta by del_beta
                CLEAR();
                AXES(0.0, 0.0, 0.5, 0.5, 1.0, 1.0);
                beta -= del_beta;
                M(-3.0, 2.0);
                TEXT(0, "beta = ", beta *180.0/M_PI);
                TEXT(0, "r = ", p[1]);
                break;
            case 'B': // increase beta by del_beta
                CLEAR();
                AXES(0.0, 0.0, 0.5, 0.5, 1.0, 1.0);
                beta += del_beta;
                M(-3.0, 2.0);
                TEXT(0, "beta = ", beta *180.0/M_PI);
                TEXT(0, "r = ", p[1]);
                break;
            case 'i': // new initial conditions
                hotkeys(off);
          //      GWIN();
                gg = GIN();
                uu = g2u(gg);
                CLEAR();
                AXES(0.0, 0.0, 0.5, 0.5, 1.0, 1.0);
                M(-3.0, 2.0);
                TEXT(0, "beta = ", beta *180.0/M_PI);
                TEXT(0, "r = ", p[1]);
                u[0] = uu[0]*cos(beta);
                u[1] = uu[0]*sin(beta);
                u[2] = uu[1];
                AWIN();
                cout << "new starting point:\t";
                cout << setprecision(3);
                cout << "x = " << u[0] << '\t';
                cout << "y = " << u[1] << '\t';
                cout << "z = " << u[2] << '\n';
                hotkeys(on);
                GWIN();
                break;
            case 'c': // clear screen and redraw axes
                CLEAR();
                AXES(0.0, 0.0, 0.5, 0.5, 1.0, 1.0);
                M(-3.0, 2.0);
                TEXT(0, "beta = ", beta *180.0/M_PI);
                TEXT(0, "r = ", p[1]);
                u[2] += 0.5; // to get off fixed point if r too small
                break;
            default:
                break;
        } // end of switch
        if (quit) {
            break; // get out of "for" loop
        }
        c = 0;
        xx = u[0]*cos(beta) + u[1]*sin(beta);
        BOX(xx, u[2], 0.3);
//        DOT(xx, u[2]);
        derivs(t, u, dudt);
        rkqc(nvar, h, t, u, dudt, hdid, hnext, derivs);
        h = hnext;
    } // end of "for" loop
    hotkeys(off);
    END();
    cout << "total points plotted = " << count << '\n';
}


