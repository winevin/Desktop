#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>

struct Node
{
    int32_t key;
    int32_t left_idx;
    int32_t right_idx;
};

void trace(int fd, int now) {
    struct Node curr;
    if (lseek(fd, now * sizeof(curr), SEEK_SET) == -(off_t) 1) {
        exit(1);
    }
    if (read(fd, &curr, sizeof(curr)) != sizeof(curr)) {
        exit(0);
    }
    if (curr.right_idx) {
        trace(fd, curr.right_idx);
    }
    printf("%" PRId32 " ", curr.key);
    if (curr.left_idx) {
        trace(fd, curr.left_idx);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) {
        return 0;
    }
    trace(fd, 0);
    putchar('\n');
    if (close(fd)) {
        return 1;
    }
    return 0;
}
