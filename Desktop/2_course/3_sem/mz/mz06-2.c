#include <stdio.h>
#include <time.h>
#include <string.h>

enum
{
    START_YEAR = 1900,
    FIRST_DAY = 7,
    FIRST_MON = 10,
    FIRST_YEAR = 1925,
    SEC_IN_DAY = 86400,
    DAY_IN_MON = 30,
    MON_IN_YEAR = 12,
    SEC_IN_MON = SEC_IN_DAY * DAY_IN_MON,
    SEC_IN_YEAR = SEC_IN_MON * MON_IN_YEAR,
};

int main(void) {
    int year, mon, mday;
    struct tm tt = {0};
    tt.tm_mday = FIRST_DAY;
    tt.tm_year = FIRST_YEAR - START_YEAR;
    tt.tm_mon = FIRST_MON - 1;
    time_t start = mktime(&tt);
    while (scanf("%d%d%d", &year, &mon, &mday) == 3) {
        struct tm t = {0};
        t.tm_mday = mday;
        t.tm_year = year - START_YEAR;
        t.tm_mon = mon - 1;
        unsigned long long sec = (unsigned long long) difftime(mktime(&t), start) + SEC_IN_DAY / 2;
        unsigned long long year1 = sec / SEC_IN_YEAR;
        unsigned long long mon1 = (sec % SEC_IN_YEAR) / SEC_IN_MON;
        unsigned long long day1 = (sec % SEC_IN_MON) / SEC_IN_DAY;
        printf("%llu %llu %llu\n", year1 + 1, mon1 + 1, day1 + 1);
    }
    return 0;
}
