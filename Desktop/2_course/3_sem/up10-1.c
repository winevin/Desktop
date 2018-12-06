#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <time.h>

enum { START_YEAR = 1900 };

struct tm * mytime(int fd[]) {
        while(wait(NULL) == -1);
        time_t t;
        read(fd[0], &t, sizeof(t));
        close(fd[1]);
        close(fd[0]);
        return localtime(&t);
}

int main(void) {
    int fd[2];
    pipe(fd);
    if (!fork()) {
        if (!fork()) {
            if (!fork()) {
                time_t t = time(NULL);
                write(fd[1], &t, sizeof(t));
                write(fd[1], &t, sizeof(t));
                write(fd[1], &t, sizeof(t));
                close(fd[0]);
                close(fd[1]);
                exit(0);
            }
            struct tm *tm = mytime(fd);
            printf("D:%02d\n", tm->tm_mday);
            exit(0);
        }
        struct tm *tm = mytime(fd);
        printf("M:%02d\n", tm->tm_mon + 1);
        exit(0);
    }
    struct tm *tm = mytime(fd);
    printf("Y:%04d\n", tm->tm_year + START_YEAR);
    return 0;
}
