#include <stdio.h>

enum { MAX_LEN_MANTISSA = 23 };

int check(unsigned num) {
    int mantissa = 0;
    while ((num & 1) == 0 && num) {
        num >>= 1;
    }
    while (num) {
        mantissa++;
        num >>= 1;
    }
    return mantissa - 1 <= MAX_LEN_MANTISSA;
}

int main(void) {
    unsigned num;
    while (scanf("%u", &num) == 1) {
        printf("%d\n", check(num));
    }
    return 0;
}
