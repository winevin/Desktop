#include <stdio.h>
#include <stdlib.h>
#include "gauss.h"

enum
{
    MX = 4,
    MY = 9,
};

int main(void) {
    double matrix[MX * MX] = {
        2, 2, -1, 1,
        4, 3, -1, 2,
        8, 5, -3, 4,
        3, 3, -2, 4,
    };
//    double ans[MX] V= {0};
//    Print(matrix, MX, MX + 1);
//    double *m = Gauss(matrix, ans, MX, MY);
//    Print(m, MX, MX);
//    free(m);
    NumOfObus(matrix, MX);
//    Gauss(matrix, ans, MX, MY);
//    Print(ans, MX, MY);
//    Gauss(matrix, ans, MX, MY);
//    Print(ans, MX, 1);
    return 0;
}
