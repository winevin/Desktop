#include <stdio.h>
#include <inttypes.h>

int main(void) {
    int n, w;
    scanf ("%o%o", &n, &w);
    int32_t max = ((int32_t) 1) << n;
    for (int32_t i = 0; i < max; ++i) {
        int32_t sign_mag = i < (max >> 1) ? i : -i + (max >> 1);
        printf("|%*" PRIo32 "|%*" PRId32 "|%*" PRId32 "|\n", w, i, w, i, w, sign_mag);
    }
    return 0;
}
