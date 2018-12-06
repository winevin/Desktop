#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

enum
{
    START_YEAR = 1900,
    SUNDAY = 7,
};

int main(int argc, char *argv[]) {
    if (argc < 3) {
        return 1;
    }
    struct tm ptime = {0};
    ptime.tm_year = atoi(argv[1]) - START_YEAR;
    ptime.tm_mday = atoi(argv[2]);
    if (mktime(&ptime) == -1) {
        return 1;
    }
    printf("%d %d %d %d\n", ptime.tm_year + START_YEAR, ptime.tm_mon + 1, ptime.tm_mday,
            ptime.tm_wday ? ptime.tm_wday : SUNDAY);
    return 0;
}
