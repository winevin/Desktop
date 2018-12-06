#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

enum { START_YEAR = 1900 };

time_t get_time(FILE *f) {
    struct tm ptime = {0};
    if (fscanf(f, "%d/%d/%d%d:%d:%d", &(ptime.tm_year), &(ptime.tm_mon), &(ptime.tm_mday),
            &(ptime.tm_hour), &(ptime.tm_min), &(ptime.tm_sec)) != 6) {
        return -1;
    }
    ptime.tm_mon -= 1;
    ptime.tm_year -= START_YEAR;
    ptime.tm_isdst = -1;
    return mktime(&ptime);
}

int main(int argc, char *argv[]) {
    FILE *f = fopen(argv[1], "r");
    time_t ttime1, ttime2;
    if (f == NULL || (ttime1 = get_time(f)) == -1) {
        return 1;
    }
    while ((ttime2 = get_time(f)) != -1) {
        printf("%ld\n", ttime2 - ttime1);
        ttime1 = ttime2;
    }
    return 0;
}
