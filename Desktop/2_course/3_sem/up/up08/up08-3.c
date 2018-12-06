#include <unistd.h>
#include <stdio.h>
#include <sys/wait.h>

int main(void) {
    int a;
    setbuf(stdin, NULL);
    setbuf(stdout, NULL);
    for (int i = 0; i < 3; ++i) {
        if (!fork()) {
            scanf("%d", &a);
            printf("%d %d\n", i + 1, a * a);
            return 0;
        }
    }
    for (int i = 0; i < 3; ++i) {
        while (wait(NULL) == -1);
    }
    return 0;
}
