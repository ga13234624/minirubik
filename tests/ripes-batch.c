/* Reads one state per line and prints "state|solution" from ripes.c's search,
 * so all 3,674,160 states take one process instead of millions.
 */
#include <string.h>

#define main ripes_main
#include "../ripes.c"
#undef main

int main(void)
{
    char input[64];
    while (fgets(input, sizeof input, stdin)) {
        state_t state;
        input[strcspn(input, "\n")] = '\0';
        if (!parse_state(input, &state)) {
            fprintf(stderr, "invalid state: %s\n", input);
            return 2;
        }
        int length = solve(&state);
        if (length < 0) {
            fprintf(stderr, "no solution: %s\n", input);
            return 1;
        }
        printf("%s|", input);
        for (int i = 0; i < length; ++i)
            printf("%s%s", i ? " " : "", move_names[solution[i]]);
        putchar('\n');
    }
    return fflush(stdout) != 0 || ferror(stdout);
}
