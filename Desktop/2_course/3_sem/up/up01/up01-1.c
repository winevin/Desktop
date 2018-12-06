#include <stdio.h>
#include <stdlib.h>

int main(void) {
    double x, y;
    scanf ("%lf%lf", &x, &y);
    printf ("%d\n", x >= 1 && y >= 2 && x <= 7 && y <= 5 && y <= -x + 10);
    return 0;
}
