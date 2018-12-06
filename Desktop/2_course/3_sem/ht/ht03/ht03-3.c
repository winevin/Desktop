#include <stdio.h>
#include <stdlib.h>

enum
{
    A = 1103515245,
    C = 12345,
    MOD = 1ULL << 31,
};

typedef struct RandomOperations
{
    void (*destroy)(void *);
    int (*next)(void *);
} RandomOperations;

typedef struct RandomGenerator
{
    RandomOperations *ops;
    int now;
} RandomGenerator;

void DestroyOperations(RandomGenerator *rr) {
    free(rr->ops);
    free(rr);
}

int NextRand(RandomGenerator *rr) {
    rr->now = (((unsigned long long) rr->now) * A + C) % (unsigned long long) MOD;
    return rr->now;
}

RandomGenerator *random_create(int seed) {
    RandomGenerator *res = malloc(sizeof(*res));
    if (res) {
        res->now = seed;
        res->ops = malloc(sizeof(*(res->ops)));
        if (!res->ops) {
            free(res);
            return NULL;
        }
        res->ops->next = (int(*) (void *)) NextRand;
        res->ops->destroy = (void(*) (void *)) DestroyOperations;
    }
    return res;
}
