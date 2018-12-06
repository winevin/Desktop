#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>

#define TIMER 1

void alrm(int sig) {
    static int cnt = 0;
    system("sudo rfkill block wifi");
    system("sudo rfkill unblock wifi");
    if (sig == SIGALRM) {
        alarm(1200);
    } else {
        if (cnt < 3) {
            ++cnt;
        } else
            signal(SIGINT, SIG_DFL);
    }
}

int main(void) {
    signal(SIGINT, alrm);
#if TIMER
    signal(SIGALRM, alrm);
    alarm(1200);
    while (1)
        sleep(1200);
#else
    alrm(1);
#endif
    return 0;
}
