#include <stdio.h>
#include <stdlib.h>

enum { AMOUNT_BYTE = 8 };

int main(void) {
    unsigned byte[AMOUNT_BYTE];
    while (scanf("%x", byte) == 1) {
        for (int i = 1; i < AMOUNT_BYTE; i++) {
            scanf("%x", byte + i);
        }
        for (int i = 0; i < AMOUNT_BYTE; i++) {
            unsigned t = 0;
            for (int j = 0; j < AMOUNT_BYTE; j++) {
                t |= ((byte[j] >> i) & 1) << j;
            }
            printf("%02X ", t);
        }
    }
    putchar('\n');
    return 0;
}
