#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

enum { STOP = 4 };

unsigned long long cur = 0;

void mysig(int sig) {
    signal(sig, mysig);
    static int cnt = 0;
    if (sig == SIGINT && ++cnt < STOP) {
        printf("%llu\n", cur);
    } else {
        exit(0);
    }
}

int main(void) {
    signal(SIGINT, mysig);
    signal(SIGTERM, mysig);
    setbuf(stdout, NULL);
    printf("%d\n", getpid());
    unsigned long long low, high;
    scanf("%llu%llu", &low, &high);
    for (unsigned long long j = low; j < high; ++j) {
        int isprime = j > 1;
        for (unsigned long long i = 2; i * i < j + 1 && isprime; ++i) {
            isprime &= j % i != 0; 
        }
        if (isprime) {
            cur = j;
        }
    }
    printf("-1\n");
    return 0;
}
