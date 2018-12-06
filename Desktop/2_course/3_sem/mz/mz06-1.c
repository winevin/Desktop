#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    off_t size = 0;
    for (int i = 1; i < argc; ++i) {
        struct stat stb;
        if (lstat(argv[i], &stb) == 0 && S_ISREG(stb.st_mode) && stb.st_nlink == 1) {
            size += stb.st_size;
        }
    }
    printf("%lld\n", size);
    return 0;
}
