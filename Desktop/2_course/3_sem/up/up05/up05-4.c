#include <sys/stat.h>

struct Task
{
    unsigned uid;
    int gid_count;
    unsigned *gids;
};

enum
{
    FIRST_TRIPLE = 6,
    SECOND_TRIPLE = 3,
};

int myaccess(const struct stat *stb, const struct Task *task, int access) {
    if (task->uid == 0) {
        return 1;
    } else if (task->uid == stb->st_uid) {
        return ((stb->st_mode >> FIRST_TRIPLE) & access) == access;
    } else {
        for (int i = 0; i < task->gid_count; ++i) {
            if (task->gids[i] == stb->st_gid) {
                return ((stb->st_mode >> SECOND_TRIPLE) & access) == access;
            }
        }
        return (stb->st_mode & access) == access;
    }
}
