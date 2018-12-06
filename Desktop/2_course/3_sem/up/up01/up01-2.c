#include <stdio.h>
#include <ctype.h>

enum
{
    BEGIN_LOWER = 11,
    BEGIN_UPPER = 37,
    END_UPPER = 63,
    THIRD_BIT = 8,
    SECOND_BIT = 4
};

int char_to_ind(int c) {
    if (islower(c)) {
        return c - 'a' + BEGIN_LOWER;
    }
    if (isdigit(c)) {
        return c - '0' + 1;
    }
    if (isupper(c)) {
        return c - 'A' + BEGIN_UPPER;
    }
    return -1;
}

int main(void) {
    int c, ind;
    char disp[END_UPPER + 1] = "@";
    for (int i = 1; i < BEGIN_LOWER; i++) {
        disp[i] = i + '0' - 1;
    }
    for (int i = BEGIN_LOWER; i < BEGIN_UPPER; i++) {
        disp[i] = i - BEGIN_LOWER + 'a';
    }
    for (int i = BEGIN_UPPER; i < END_UPPER; i++) {
        disp[i] = i - BEGIN_UPPER + 'A';
    }
    disp[END_UPPER] = '#';
    while ((c = getchar()) != EOF)  {
        if ((ind = char_to_ind(c)) != -1) {
            putchar(disp[(ind ^ THIRD_BIT) & ~SECOND_BIT]);
        }
    }
    putchar('\n');
    return 0;
}
