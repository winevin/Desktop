#include <stdio.h>
#include <stdlib.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

enum
{
    ACCESS = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP,
    FLAGS = O_WRONLY | O_CREAT,
};

unsigned long long hash_function(unsigned long long);

int main(int argc, char *argv[]) {
    int fd = open(argv[1], FLAGS, ACCESS);
    unsigned long long h, count;
    if (fd < 0) {
        return 1;
    }
    sscanf(argv[2], "%llxu", &h);
    sscanf(argv[3], "%llu", &count);
    if (lseek(fd, (count + 1) * sizeof(h), SEEK_END) == -(off_t) 1) {
        return 1;
    }
    for (int i = 0; i < count; ++i) {
        if (lseek(fd, -(off_t) 2 * sizeof(h), SEEK_CUR) == -(off_t) 1) {
            return 1;
        }
        if (write(fd, &h, sizeof(h)) != sizeof(h)) {
            return 2;
        }
        h = hash_function(h);
    }
    if (close(fd)) {
        return 1;
    }
    return 0;
}
