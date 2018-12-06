#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

enum { BUF_SIZE = 16 };

char * get_long_num(char *buffer) {
    char *buf = buffer;
    char op = '+';
    while (*buf < '0' || *buf > '9') {// start search of first digit
        op = *buf;
        if (*buf++ == '\0') {			// end of buffer. read next path
            buf = buffer;
            int len_buf = read(0, buf, BUF_SIZE);
            if (len_buf <= 0) {
                buf[0] = '\0';
                return NULL;
            }
            buf[len_buf] = '\0';
        }
    }// end search of digit
    int len_num = 2, ind_num = 0;
    char *num = calloc(len_num + 1, sizeof(*num));
    op = op == '-' ? '-' : '+';
    num[ind_num++] = op;
    while (*buf >= '0' && *buf <= '9') {// start get a big number
        if (*buf != '0' || ind_num != 1) {
            num[ind_num++] = *buf;
            if (len_num <= ind_num + 1) {
                len_num *= 2;
                char *tmp = realloc(num, (len_num + 1) * sizeof(*tmp));
                if (tmp == NULL) {
                    free(num);
                    exit(1);
                }
                num = tmp;
            }
        }
        ++buf;
        if (*buf == '\0') { //end of buffer. read next path
            buf = buffer;
            int len_buf = read(0, buf, BUF_SIZE);
            if (len_buf <= 0) {
                buf[0] = '\0';
                break;
            }
            buf[len_buf] = '\0';
        }
    }// end get a big number
    if (ind_num == 1) {
        num[ind_num++] = '0';
    }
    while (*buf) {
        *buffer++ = *buf++;
    }
    *buffer = '\0';
    num[ind_num] = '\0';
    return num; //well done
}

unsigned long long last_bits (char *num) {
    if (num == NULL) {
        return 0;
    }
    unsigned long long ans = 0;
    char op = *num++;
    while (*num) {
        ans *= 10;
        ans += *num - '0';
        ++num;
    }
    if (op == '-') {
        ans = 0 - ans;
    }
    return ans;
}

int main(void) {
    char buf[BUF_SIZE + 1];
    int len_buf = read(0, buf, BUF_SIZE);
    if (len_buf <= 0) {
        printf("0\n");
        return 0;
    }
    buf[len_buf] = '\0';
    char *num = get_long_num(buf);
    unsigned long long ans = last_bits(num);
    free(num);
    while ((num = get_long_num(buf)) != NULL) {
        ans += last_bits(num);
        free(num);
    }
    printf("%lld\n", ans);
    return 0;
}
