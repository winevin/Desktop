#include <stdio.h>
#include <stdlib.h>

enum { MAX_MOD = 2000 };

int main(void) {
    int n;
    scanf("%d", &n);
    if (n >= MAX_MOD || n <= 1) {
        return 1;
    }
    unsigned short *mult = calloc(n * n, sizeof(*mult));
    if (mult == NULL) {
        return 2;
    }
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            mult[a * n + (a * b) % n] = b;
        }
    }
    for (int c = 0; c < n; c++, putchar('\n')) {
        for (int a = 1; a < n; a++) {
            printf("%d ", mult[a * n + c]);
        }
    }
    free(mult);
    return 0;
}
