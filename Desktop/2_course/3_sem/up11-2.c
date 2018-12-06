#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

int flag = SIGINT;

void mysig(int sig) {
    flag = sig;
}

int main(void) {
    int t;
    unsigned sum = 0;
    signal(SIGINT, mysig); 
    signal(SIGQUIT, mysig);
    setbuf(stdout, NULL);
    printf("%d\n", getpid());
    while (scanf("%d", &t) == 1) {
        printf("%d\n", sum = flag == SIGINT ? sum + t : sum * t);
    }
    return 0;
}
