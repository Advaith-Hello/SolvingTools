#ifndef COMPUTING_ODE_UTILS_H
#define COMPUTING_ODE_UTILS_H


typedef double (*UnaryFunc)(double x);
typedef double (*BiVarFunc)(double x, double y);

typedef struct {
    int size;
    double r_start;
    double r_end;
    double data[];
} FunctionTable;

typedef struct {
    double val;
    double loss;
} LossPair;

typedef struct {
    double x;
    double y;
} Coord;

double eval_func(const FunctionTable *table, double x);
double deriv(UnaryFunc f, double x);

LossPair var_solve(
    UnaryFunc lhs, 
    UnaryFunc rhs, 
    double start, 
    double end, 
    int amt);

FunctionTable* ode_solve(
    BiVarFunc f,
    Coord anchor,
    double r_start, 
    double r_end, 
    int r_amt);


#endif // COMPUTING_ODE_UTILS_H
