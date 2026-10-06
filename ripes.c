#include <stdint.h>

#if defined(__riscv) && !defined(__linux__)
#define BARE_METAL 1
#else
#define BARE_METAL 0
#include <stdio.h>
#endif

/* Ripes passes no command line, so the bare-metal build solves this state.
 * The format is solver.c's: seven cubie digits, then seven twist digits. */
#ifndef SCRAMBLE
#define SCRAMBLE "21345671111111"
#endif

enum {
    CUBIES = 7,
    FACES = 3,
    MOVES = 9,
    MAX_MOVES = 11 /* God's number for the 2x2x2 in the half-turn metric */
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} state_t;

static const char *const move_names[MOVES] = {"R",  "R2", "R'", "B", "B2",
                                              "B'", "D",  "D2", "D'"};
static const uint8_t inverse_move[MOVES] = {2, 1, 0, 5, 4, 3, 8, 7, 6};

/* permutation_move, orientation_move, permutation_depth and
 * orientation_depth: 30,240 + 4,374 + 5,040 + 729 = 40,383 bytes of .rodata. */
#include "ripes_tables.h"

static uint8_t solution[MAX_MOVES];

/* The tables are indexed by solver.c's rank_state, so these two must agree
 * with it. This is its Lehmer rank, adding each inversion's factorial weight
 * instead of multiplying, which keeps mul and __mulsi3 out of RV32I code. */
static uint16_t rank_permutation(const uint8_t p[CUBIES])
{
    static const uint16_t weight[CUBIES - 1] = {720, 120, 24, 6, 2, 1};
    uint16_t rank = 0;
    for (uint8_t i = 0; i + 1 < CUBIES; ++i)
        for (uint8_t j = (uint8_t) (i + 1); j < CUBIES; ++j)
            if (p[j] < p[i])
                rank = (uint16_t) (rank + weight[i]);
    return rank;
}

/* Base-3 rank of the first six twists; the seventh is implied by parity. */
static uint16_t rank_orientation(const uint8_t o[CUBIES])
{
    uint16_t rank = 0;
    for (uint8_t i = 0; i + 1 < CUBIES; ++i)
        rank = (uint16_t) (rank * 3U + o[i]);
    return rank;
}

/* verify.c defines this to count the nodes the search visits. */
#ifndef COUNT_NODE
#define COUNT_NODE() ((void) 0)
#endif

/* Depth-first search for a path of exactly `left` moves from (p, o) to solved,
 * never turning the same face twice in a row, since an optimal path cannot.
 */
static int search(uint16_t p, uint16_t o, uint8_t left, uint8_t last_face)
{
    COUNT_NODE();
    if (left == 0)
        return p == 0 && o == 0;
    for (uint8_t face = 0; face < FACES; ++face) {
        uint16_t next_p = p, next_o = o;
        if (face == last_face)
            continue;
        for (uint8_t turn = 0; turn < 3; ++turn) {
            next_p = permutation_move[next_p][face];
            next_o = orientation_move[next_o][face];
            if (permutation_depth[next_p] >= left ||
                orientation_depth[next_o] >= left)
                continue;
            solution[left - 1] = inverse_move[face * 3U + turn];
            if (search(next_p, next_o, (uint8_t) (left - 1), face))
                return 1;
        }
    }
    return 0;
}

/* Fills solution[] and returns its length, or -1 if no solution exists.
 *
 * solver.c prints, backwards and inverted, the path by which its breadth-first
 * search first reached the state, and that is the lexicographically least
 * shortest path from solved. The same moves lead from the inverse state back
 * to solved, so a depth-first search from the inverse that tries moves in
 * index order meets that path first; storing each move's inverse at index
 * left - 1 reverses it on the way. The output is therefore solver.c's, byte for
 * byte, not merely another optimal solution.
 */
static int solve(const state_t *state)
{
    state_t inverse;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t cubie = state->p[i];
        inverse.p[cubie] = i;
        inverse.o[cubie] = (uint8_t) (state->o[i] ? 3 - state->o[i] : 0);
    }
    uint16_t p = rank_permutation(inverse.p);
    uint16_t o = rank_orientation(inverse.o);
    uint8_t bound = permutation_depth[p] > orientation_depth[o]
                        ? permutation_depth[p]
                        : orientation_depth[o];
    for (; bound <= MAX_MOVES; ++bound)
        if (search(p, o, bound, FACES))
            return bound;
    return -1;
}

/* solver.c's input rules: a permutation of 1-7, then twists 1-3 whose
 * internal values sum to a multiple of three. Reads no further than the first
 * bad character, so a short string is never overrun. */
static int parse_state(const char *input, state_t *state)
{
    uint8_t seen = 0, sum = 0;
    for (uint8_t i = 0; i < 2 * CUBIES; ++i) {
        uint8_t digit = (uint8_t) (input[i] - '1');
        if (digit >= (i < CUBIES ? CUBIES : 3))
            return 0;
        if (i < CUBIES) {
            if (seen >> digit & 1U)
                return 0;
            seen = (uint8_t) (seen | 1U << digit);
            state->p[i] = digit;
        } else {
            sum = (uint8_t) (sum + digit);
            state->o[i - CUBIES] = digit;
        }
    }
    while (sum >= 3)
        sum = (uint8_t) (sum - 3);
    return input[2 * CUBIES] == '\0' && sum == 0;
}

#if BARE_METAL
/* Ripes jumps to the ELF entry with no C runtime behind it, and its command
 * line mode leaves sp at 0. Set the global pointer the linker relaxes against,
 * point sp at a stack of our own, and pass main's result to exit (ecall 93).
 * Nothing reads .bss before writing it, so it is not cleared.
 */
__asm__(
    ".pushsection .text._start, \"ax\", @progbits\n"
    ".globl _start\n"
    "_start:\n"
    ".option push\n"
    ".option norelax\n"
    "    la gp, __global_pointer$\n"
    ".option pop\n"
    "    la sp, stack_top\n"
    "    call main\n"
    "    li a7, 93\n"
    "    ecall\n"
    "1:  j 1b\n"
    ".popsection\n"
    ".pushsection .bss.stack, \"aw\", @nobits\n"
    ".balign 16\n"
    "    .space 2048\n"
    "stack_top:\n"
    ".popsection\n");

/* write(fd, text, length) through Ripes's Linux-numbered system call 64. */
static void emit(uint32_t fd, const char *text, uint32_t length)
{
    register uint32_t a0 __asm__("a0") = fd;
    register const char *a1 __asm__("a1") = text;
    register uint32_t a2 __asm__("a2") = length;
    register uint32_t a7 __asm__("a7") = 64;
    __asm__ volatile("ecall" : "+r"(a0) : "r"(a1), "r"(a2), "r"(a7) : "memory");
}
#else
static void emit(uint32_t fd, const char *text, uint32_t length)
{
    fwrite(text, 1, length, fd == 1 ? stdout : stderr);
}
#endif

#define EMIT(fd, literal) emit(fd, literal, sizeof(literal) - 1)

static int print_solution(const state_t *state)
{
    char line[3 * MAX_MOVES + 1], *end = line;
    int length = solve(state);
    if (length < 0) {
        EMIT(2, "no solution within 11 moves\n");
        return 1;
    }
    for (int i = 0; i < length; ++i) {
        const char *name = move_names[solution[i]];
        if (i > 0)
            *end++ = ' ';
        while (*name)
            *end++ = *name++;
    }
    *end++ = '\n';
    emit(1, line, (uint32_t) (end - line));
    return 0;
}

#if BARE_METAL
int main(void)
{
    state_t state;
    if (!parse_state(SCRAMBLE, &state)) {
        EMIT(2, "SCRAMBLE is not a valid PPPPPPPOOOOOOO state\n");
        return 2;
    }
    return print_solution(&state);
}
#else
int main(int argc, char **argv)
{
    state_t state;
    if (argc != 2 || !parse_state(argv[1], &state)) {
        fprintf(stderr, "usage: %s PPPPPPPOOOOOOO\n",
                argc > 0 && argv[0] ? argv[0] : "ripes");
        return 2;
    }
    int status = print_solution(&state);
    /* As in solver.c, a write error surfaces at the flush. */
    return status ? status : fflush(stdout) != 0 || ferror(stdout);
}
#endif
