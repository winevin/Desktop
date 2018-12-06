#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <unistd.h>
#include <sys/wait.h>

int main(void) {
    int a = 0, pid = 0, parent = getpid();
    setbuf(stdin, NULL);
    setbuf(stdout, NULL);
    while (!pid && scanf("%d", &a) == 1) {
        pid = fork();
        if (pid == -1) {
            printf("-1\n");
            return 1;
        }
    }
    int status;
    wait(&status);
    if (WEXITSTATUS(status) == 1) {
        return parent != getpid();
    } else if (pid) {
        printf("%d\n", a);
    }
    return 0;
}
