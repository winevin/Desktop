#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>

enum
{
    ONLY_SPACE = 0x0bad1dea,
    NOT_NUM = 0xbedabeda,
    BASIS = 10,
};

char * getline2(FILE *f);

int main(void) {
    char *s;
    unsigned str_ind = 0;
    while ((s = getline2(stdin)) != NULL) {
        char *tmp = s;
        ++str_ind;
        unsigned sum = 0, num_ind = 0;
        for (; *s; ++s) {
            if ((*s >= '0' && *s <= '9') || *s == '-' || *s == '+') {
                ++num_ind;
                char *s_end;
                errno = 0;
                long t = strtol(s, &s_end, BASIS);
                if (!s_end || !isspace((unsigned char) *s_end)) {
                    break;
                } if (errno == ERANGE || t != (int) t) {
                    t = *s == '-' ? -num_ind : num_ind;
                }
                sum = sum + (unsigned) t;
                s = s_end;
            } else if (!isspace((unsigned char) *s)) {
                break;
            }
        }
        if (*s) {
            printf("%d\n", NOT_NUM + str_ind);
        } else {
            printf("%d\n", num_ind ? sum : ONLY_SPACE + str_ind);
        }
        free(tmp);
    }
    return 0;
}
