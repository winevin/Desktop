#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum
{
    PERCENT = 100, //100%
    PRECISION = 10000, // Точность округления до 4 знаков
};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }
    double rate = atof(argv[1]);
    for (int i = 2; i < argc; i++) {
        rate += rate * atof(argv[i]) / PERCENT;
        rate = round(rate * PRECISION) / PRECISION;
    }
    printf("%.4lf\n", rate);
    return 0;
}
