#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }
    int fd = open(argv[1], O_RDWR);
    if (fd < 0) {
        return 1;
    }
    signed long long t, min;
    signed long long offset_min = 0, offset = sizeof(min);
    if (read(fd, &min, sizeof(min)) != sizeof(min)) {
        return 0;
    }
    while (read(fd, &t, sizeof(t)) == sizeof(t)) {
        if (t < min) {
            min = t;
            offset_min = offset;
        }
        offset += sizeof(t);
    }
    lseek(fd, offset_min, SEEK_SET);
    min = 0 - (unsigned long long) min;
    if (write(fd, &min, sizeof(min)) != sizeof (min)) {
        return 1;
    }
    close(fd);
    return 0;
}
