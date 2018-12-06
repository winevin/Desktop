#include <stdio.h>

int main(int argc, char *argv[]) {
    const char MASK[] = "rwxrwxrwx";
    for (int i = 1; i < argc; ++i) {
        int bits;
        sscanf(argv[i], "%o", &bits);
        for (int k = sizeof(MASK) - 2; k >= 0; --k) {
            putchar((bits >> k) & 1 ? MASK[sizeof(MASK) - 2 - k] : '-');
        }
        putchar('\n');
    }
    return 0;
}
