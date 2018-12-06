#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

enum
{
    STOP = 6,
    SIG = SIGHUP,
};

int cnt = 0;

void mysig(int sig) {
    signal(SIG, mysig);
    if (cnt >= STOP - 1) {
        exit(0);
    }
    printf("%d\n", cnt++);
    fflush(stdout);
}

int main(void) {
    signal(SIG, mysig);
    printf("%d\n", getpid());
    fflush(stdout);
    while (1) {
        pause();
    }
    exit(1);
}
