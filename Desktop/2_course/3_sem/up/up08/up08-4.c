#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

int main(void) {
    int parent = getpid(), pid = 0, i = 0, n;
    setbuf(stdout, NULL);
    scanf("%d", &n);
    printf("%d", ++i);
    while (i < n && !pid) {
        printf(" %d", ++i);
        pid = fork();
    }
    wait(NULL);
    if (parent == getpid()) {
        printf("\n");
    }
    return 0;
}
