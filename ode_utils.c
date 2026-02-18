#include "ode_utils.h"

#include <math.h>
#include <stdlib.h>


double eps = 1e-9;

double eval_func(const FunctionTable *table, const double x) {
    const double step = (table->r_end - table->r_start) / table->size;
    const double pos = (x - table->r_start) / step;
    const int idx = (int)pos;
    const double t = pos - idx;
    return table->data[idx] * (1 - t) + table->data[idx + 1] * t;
}

double deriv(const UnaryFunc f, const double x) {
    return (f(x + eps) - f(x - eps)) / (2 * eps);
}

LossPair var_solve(
    const UnaryFunc lhs,
    const UnaryFunc rhs,
    const double start,
    const double end,
    const int amt) {

    LossPair best = {0, INFINITY};

    const double step = (end - start) / amt;
    for (int i = 0; i <= amt; i++) {
        const double x = start + i * step;
        const double c_loss = fabs(rhs(x) - lhs(x));

        if (c_loss < best.loss) {
            best.val = x;
            best.loss = c_loss;
        }
    }

    return best;
}

FunctionTable* ode_solve(
    const BiVarFunc f,
    const Coord anchor,
    const double r_start,
    const double r_end,
    const int r_amt) {

    FunctionTable* table = malloc(sizeof(FunctionTable) + r_amt * sizeof(double));

    table->size = r_amt;
    table->r_start = r_start;
    table->r_end = r_end;

    const double h = (r_end - r_start) / (double) r_amt;
    const int pos = lround((anchor.x - r_start) / h);

    double x = anchor.x;
    double y = anchor.y;

    table->data[pos] = anchor.y;
    for (int i = pos + 1; i < r_amt; i++) {
        const double k1 = f(x, y);
        const double k2 = f(x + h/2, y + h/2 * k1);
        const double k3 = f(x + h/2, y + h/2 * k2);
        const double k4 = f(x + h, y + h * k3);

        y += h/6 * (k1 + 2*k2 + 2*k3 + k4);
        table->data[i] = y;
        x += h;
    }

    x = anchor.x;
    y = anchor.y;

    for (int i = pos - 1; i >= 0; i--) {
        const double k1 = f(x, y);
        const double k2 = f(x - h/2, y - h/2 * k1);
        const double k3 = f(x - h/2, y - h/2 * k2);
        const double k4 = f(x - h, y - h * k3);

        y -= h/6 * (k1 + 2*k2 + 2*k3 + k4);
        table->data[i] = y;
        x -= h;
    }

    return table;
}
