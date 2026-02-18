#include <stdio.h>
#include <stdlib.h>

#include "ode_utils.h"


double f(const double x, const double y) {
    return x * y;
}

int main(void) {
    FunctionTable* table = ode_solve(f, (Coord){0, 1}, 0, 10, 10000000);
    printf("%f\n", eval_func(table, 5));
    free(table);
    return 0;
}
