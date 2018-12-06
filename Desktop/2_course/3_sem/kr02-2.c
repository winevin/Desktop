#include <stdio.h>
#include <time.h>

enum
{
    IDX_YEAR = 1,
    IDX_MDAY,
    IDX_WDAY,
    START_YEAR = 1900,
    HALF_DAY = 12,
    SUNDAY = 7,
    DECEMBER = 12,
};

int WeekDay(int mday, int mon, int year) {
    struct tm t = {
        .tm_mday = mday,
        .tm_mon = mon,
        .tm_year = year - START_YEAR,
        .tm_hour = HALF_DAY
    };
    mktime(&t);
    return t.tm_wday ? t.tm_wday : SUNDAY;
}

int main(int argc, char *argv[]) {
    int year, mday, wday, cnt = 0;
    sscanf(argv[IDX_YEAR], "%d", &year);
    sscanf(argv[IDX_MDAY], "%d", &mday);
    sscanf(argv[IDX_WDAY], "%d", &wday);

    for (int i = 0; i < DECEMBER; ++i) {
        cnt += WeekDay(mday, i, year) == wday;
    }
    printf("%d\n", cnt);
    return 0;
}
