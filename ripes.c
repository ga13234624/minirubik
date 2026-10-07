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
static uint32_t rank_permutation(const uint8_t p[CUBIES])
{
    static const uint16_t weight[CUBIES - 1] = {720, 120, 24, 6, 2, 1};
    uint32_t rank = 0;
    for (uint32_t i = 0; i + 1 < CUBIES; ++i)
        for (uint32_t j = i + 1; j < CUBIES; ++j)
            if (p[j] < p[i])
                rank += weight[i];
    return rank;
}

/* Base-3 rank of the first six twists; the seventh is implied by parity. */
static uint32_t rank_orientation(const uint8_t o[CUBIES])
{
    uint32_t rank = 0;
    for (uint32_t i = 0; i + 1 < CUBIES; ++i)
        rank = rank * 3U + o[i];
    return rank;
}

/* verify.c defines this to count the nodes the search visits. */
#ifndef COUNT_NODE
#define COUNT_NODE() ((void) 0)
#endif

/* Depth-first search for a path of exactly `bound` moves from (p, o) to
 * solved, never turning the same face twice in a row, since an optimal path
 * cannot.
 *
 * It does not recurse. The node being expanded lives in locals the compiler
 * keeps in registers: its position (p, o), `left` moves from the goal, the
 * face it is turning, how many quarter turns of that face it has tried and the
 * position (np, no) they reached. Only descending to a child touches memory,
 * saving the parent's p, o, face and turn at level `left`; backing up loads
 * them again, and the parent's (np, no) is the child it is leaving. The nodes
 * are visited in the order a recursive search would visit them, so the first
 * solution found is the same. saved_face[left + 1] is the parent's face, the
 * one this node must not turn, and saved_face[bound + 1] stands in for the
 * root's parent with FACES, which matches none.
 */
static int search(uint32_t p, uint32_t o, uint32_t bound)
{
    uint16_t saved_p[MAX_MOVES + 1], saved_o[MAX_MOVES + 1];
    uint8_t saved_face[MAX_MOVES + 2], saved_turn[MAX_MOVES + 1];
    uint32_t np = p, no = o;
    uint32_t left = bound, face = 0, turn = 0, last_face = FACES;

    COUNT_NODE();
    if (bound == 0)
        return p == 0 && o == 0;
    saved_face[bound + 1] = FACES;
    for (;;) {
        /* This face is done, or was the parent's: try the next, and once all
         * three are done, back up to the parent. */
        if (face == last_face || turn == 3) {
            if (++face == FACES) {
                if (++left > bound)
                    return 0;
                np = p;
                no = o;
                p = saved_p[left];
                o = saved_o[left];
                face = saved_face[left];
                turn = saved_turn[left];
                last_face = saved_face[left + 1];
                continue;
            }
            turn = 0;
            np = p;
            no = o;
            continue;
        }
        np = permutation_move[np][face];
        no = orientation_move[no][face];
        ++turn;
        if (permutation_depth[np] >= left || orientation_depth[no] >= left)
            continue;
        solution[left - 1] = inverse_move[face * 3U + turn - 1U];
        COUNT_NODE();
        if (left == 1) {
            if (np == 0 && no == 0)
                return 1;
            continue;
        }
        saved_p[left] = (uint16_t) p;
        saved_o[left] = (uint16_t) o;
        saved_face[left] = (uint8_t) face;
        saved_turn[left] = (uint8_t) turn;
        --left;
        p = np;
        o = no;
        last_face = face;
        face = 0;
        turn = 0;
    }
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
    for (uint32_t i = 0; i < CUBIES; ++i) {
        uint32_t cubie = state->p[i];
        inverse.p[cubie] = (uint8_t) i;
        inverse.o[cubie] = (uint8_t) (state->o[i] ? 3 - state->o[i] : 0);
    }
    uint32_t p = rank_permutation(inverse.p);
    uint32_t o = rank_orientation(inverse.o);
    uint32_t bound = permutation_depth[p] > orientation_depth[o]
                        ? permutation_depth[p]
                        : orientation_depth[o];
    for (; bound <= MAX_MOVES; ++bound)
        if (search(p, o, bound))
            return (int) bound;
    return -1;
}

/* solver.c's input rules: a permutation of 1-7, then twists 1-3 whose
 * internal values sum to a multiple of three. Reads no further than the first
 * bad character, so a short string is never overrun. */
static int parse_state(const char *input, state_t *state)
{
    uint32_t seen = 0, sum = 0;
    for (uint32_t i = 0; i < 2 * CUBIES; ++i) {
        uint32_t digit = (uint32_t) (input[i] - '1');
        if (digit >= (i < CUBIES ? CUBIES : 3))
            return 0;
        if (i < CUBIES) {
            if (seen >> digit & 1U)
                return 0;
            seen |= 1U << digit;
            state->p[i] = (uint8_t) digit;
        } else {
            sum += digit;
            state->o[i - CUBIES] = (uint8_t) digit;
        }
    }
    while (sum >= 3)
        sum -= 3;
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
