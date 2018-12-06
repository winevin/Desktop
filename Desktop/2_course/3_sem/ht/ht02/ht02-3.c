#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

enum
{
    LEN_INT = sizeof(unsigned),
    THREE_BYTE = CHAR_BIT * (sizeof(unsigned) - 1),
    ACCESS = S_IRUSR | S_IWUSR
};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }
    int fd = open(argv[1], O_WRONLY | O_TRUNC | O_CREAT, ACCESS);
    if (fd == -1) {
        return 2;
    }
    unsigned t;
    while (scanf("%u", &t) == 1) {
        char s[LEN_INT];
        for (int i = 0; i < LEN_INT; i++) {
            s[i] = t >> (THREE_BYTE - i * CHAR_BIT);
        }
        int len = LEN_INT;
        while (len > 0) {
            int check = write(fd, s + LEN_INT - len, len * sizeof(*s));
            len -= check;
            if (check <= 0) {
                return 1;
            }
        }
    }
    if (close(fd)) {
        return 1;
    }
    return 0;
}
