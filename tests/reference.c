/* Prints every one of the 3,674,160 valid states as "state|solution", using
 * solver.c's breadth-first table: the oracle for "make check-all".
 */
#define main solver_main
#include "../solver.c"
#undef main

int main(void)
{
    uint8_t diameter;
    uint8_t *table = build_table(&diameter);
    if (!table) {
        fputs("could not build complete state table\n", stderr);
        return 1;
    }
    for (uint32_t rank = 0; rank < STATES; ++rank) {
        state_t state;
        char digits[2 * CUBIES + 1];
        unrank_state(rank, &state);
        for (uint8_t i = 0; i < CUBIES; ++i) {
            digits[i] = (char) ('1' + state.p[i]);
            digits[i + CUBIES] = (char) ('1' + state.o[i]);
        }
        digits[2 * CUBIES] = '\0';
        printf("%s|", digits);
        const char *separator = "";
        for (uint32_t r = rank; r; r = rank_state(&state)) {
            printf("%s%s", separator, move_names[table[r]]);
            separator = " ";
            state = apply_move(state, table[r]);
        }
        putchar('\n');
    }
    free(table);
    return fflush(stdout) != 0 || ferror(stdout);
}
